#include "viewer/SimViewer.h"

#include <qopengl.h>
#include <QOpenGLShaderProgram>
#include <QMouseEvent>
#include <QColorDialog>
#include <QOpenGLFramebufferObject>
#include <QGridLayout>

#include <chrono>
#include <iostream>
#include <utility>

#include "contact/Contact.h"
#include "util/MeshUtils.h"
#include "viewer/RigidBodyRenderer.h"
#include "rigidbody/RigidBodySystem.h"
#include "rigidbody/Scenarios.h"

using namespace std;

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

namespace
{
    Eigen::Vector3f sImpulseDirs[8] = {
        Eigen::Vector3f(1, 0, 0),
        Eigen::Vector3f(0, 0, -1),
        Eigen::Vector3f(-1, 0, 0),
        Eigen::Vector3f(0, 0, 1),
        Eigen::Vector3f(0.71f, 0, 0.71f),
        Eigen::Vector3f(0.71f, 0, -0.71f),
        Eigen::Vector3f(-0.71f, 0, 0.71f),
        Eigen::Vector3f(-0.71f, 0, -0.71f)
    };

    static inline QImage createSubImage(const QImage* image, const QRect & rect)
    {
        const size_t offset = rect.x() * image->depth() / 8
                + rect.y() * image->bytesPerLine();
        return QImage(image->bits() + offset, rect.width(), rect.height(),
                      image->bytesPerLine(), image->format());
    }

    static void writeGreyscaleImage(QImage image, int color, const QString& filename)
    {
        for(int i = 0; i < image.width(); ++i)
        {
            for(int j = 0; j < image.height(); ++j)
            {
                QRgb c = image.pixel(i, j);
                switch(color)
                {
                case 0:
                    c = qRgba(4*qRed(c),4*qRed(c),4*qRed(c),255);
                    break;
                case 1:
                    c = qRgba(4*qGreen(c),4*qGreen(c),4*qGreen(c),255);
                    break;
                case 2:
                    c = qRgba(4*qBlue(c),4*qBlue(c),4*qBlue(c),255);
                    break;
                case 3:
                    c = qRgba(4*qAlpha(c),4*qAlpha(c),4*qAlpha(c),255);
                    break;
                }

                image.setPixel(i,j, c);

            }
        }

        image.save(filename, "PNG");
    }

    static RigidBody* findClosestBody(const std::vector<RigidBody*>& bodies, const Eigen::Vector3f& p, float& minDist, Eigen::Vector3f& pout)
    {
        minDist = std::numeric_limits<float>::max();
        RigidBody* closestBody = nullptr;
        for(RigidBody* b : bodies)
        {
            const Eigen::Matrix3f RT = b->R.transpose();
            const Eigen::Vector3f plocal = RT * (p - b->x);

            MeshCollisionData meshData;
            if( b->mesh )
            {
                barycentricClosestPoint(*(b->mesh), plocal, meshData, nullptr);

                if( meshData.distance < minDist )
                {
                    closestBody = b;
                    pout = meshData.p;
                    minDist = meshData.distance;
                }
            }
        }

        return closestBody;
    }

    static bool doPicking(qglviewer::Camera* camera, const QPoint& mouseP, const std::vector<RigidBody*>& bodies, SimViewer::PickingData& pick)
    {
        bool found = false;
        const auto vec = camera->pointUnderPixel(mouseP, found);
        if( found )
        {
            float minDist = std::numeric_limits<float>::max();
            pick.body = findClosestBody(bodies, Eigen::Vector3f(vec[0], vec[1], vec[2]), minDist, pick.plocal);
            if( minDist < 0.1f )
            {
                qglviewer::Vec cameraVec;
                camera->getCameraCoordinatesOf(vec, cameraVec);
                pick.depth = -cameraVec[2];

                qglviewer::Vec src;
                src[0] = mouseP.x();
                src[1] = mouseP.y();
                src[2] = camera->zFar() / (camera->zFar() - camera->zNear()) * (1.0 - camera->zNear() / pick.depth);
                pick.cursor = camera->unprojectedCoordinatesOf(src);
            }
            else
            {
                pick.body = nullptr;
            }

        }
        return found;
    }

}

SimViewer::SimViewer() : QGLViewer(),
    m_gl(),
    m_dt(0.01f), m_subSteps(1),
    m_paused(false), m_stepOnce(false),
    m_rigidBodySystem(), m_rigidBodyRenderer()
{
}

SimViewer::~SimViewer()
{
}

void SimViewer::cleanup()
{
    makeCurrent();
    doneCurrent();
}

void SimViewer::init()
{
    // Setup camera and mouse binding for camera control.
    camera()->setType(qglviewer::Camera::PERSPECTIVE);
    setMouseBinding(Qt::NoModifier, Qt::LeftButton, CAMERA, LOOK_AROUND);
    setMouseBinding(Qt::AltModifier, Qt::LeftButton, CAMERA, ROTATE);
    setMouseBinding(Qt::NoModifier, Qt::MouseButton(Qt::LeftButton + Qt::MidButton), CAMERA, NO_MOUSE_ACTION);
    setMouseBinding(Qt::ControlModifier, Qt::MouseButton(Qt::LeftButton + Qt::MidButton), CAMERA, NO_MOUSE_ACTION);

    // Init OpenGL objects
    connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &SimViewer::cleanup);
    m_gl.initializeOpenGLFunctions();

    qDebug() << "-- Initialized OpenGL ----------";
    qDebug () << (const char*)m_gl.glGetString(GL_VENDOR);
    qDebug () << (const char*)m_gl.glGetString(GL_VERSION);
    qDebug () << (const char*)m_gl.glGetString(GL_RENDERER);
    qDebug() << "--------------------------------";

    // Enable the depth test
    m_gl.glEnable(GL_DEPTH_TEST);
    m_gl.glDepthMask(GL_TRUE);
    m_gl.glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    m_gl.glEnable(GL_MULTISAMPLE);
    m_gl.glEnable(GL_BLEND);
    m_gl.glEnable(GL_POINT_SMOOTH);
    m_gl.glEnable(GL_LINE_SMOOTH);

    m_gl.glClearColor(0, 0, 0, 0);
    m_gl.glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Create the rigid body system and renderer
    m_rigidBodySystem.reset(new RigidBodySystem);
    m_rigidBodyRenderer.reset(new RigidBodyRenderer(&m_gl, m_rigidBodySystem.get()));

    // Set camera parameters to view the whole cloth.
    setSceneRadius(2.0);
    camera()->setPosition(qglviewer::Vec(0,1,1));
    camera()->lookAt(qglviewer::Vec(0,0,0));
    camera()->setZNearCoefficient(0.001);
    camera()->setZClippingCoefficient(10.0);
    showEntireScene();

    // Allow mouse interaction.
    setMouseTracking(true);

    // Animate.
    setAnimationPeriod(10);
    startAnimation();

    // Add pre-step hook.
    m_rigidBodySystem->setPreStepFunc(std::bind(&SimViewer::preStep, this, std::placeholders::_1));

    LightingParams lighting;
    lighting.pos = QVector4D(-20.0f, 50.0f, 0.0f, 1.0f);
    lighting.ambient = QVector3D(0.05f, 0.05f, 0.05f);
    m_rigidBodyRenderer->setLightingParameters( lighting );

    // Setup the mouse spring shader.
    {
        // Init shaders
        m_mouseShader.program = new QOpenGLShaderProgram;
        if (! m_mouseShader.program->addShaderFromSourceFile(QOpenGLShader::Vertex, "glsl/spring.vert")) {
            qDebug() << "Unable to load shader" << endl
                     << "Log file:" << endl;
            qDebug() <<  m_mouseShader.program->log();
        }
        if (! m_mouseShader.program->addShaderFromSourceFile(QOpenGLShader::Fragment, "glsl/spring.frag")) {
            qDebug() << "Unable to load shader" << endl
                     << "Log file:" << endl;
            qDebug() <<  m_mouseShader.program->log();
        }
        m_mouseShader.program->link();
        m_mouseShader.program->bind();

        // The strings "vPosition", "mvMatrix", etc. have to match an attribute name in the vertex shader.
        QString shaderParameter;
        shaderParameter = "vPosition";
        if ((m_mouseShader.vPositionLoc = m_mouseShader.program->attributeLocation(shaderParameter)) < 0)
            qDebug() << "Unable to find shader location for " << shaderParameter;

        shaderParameter = "projMatrix";
        if ((m_mouseShader.projMatrixLoc = m_mouseShader.program->uniformLocation(shaderParameter)) < 0)
            qDebug() << "Unable to find shader location for " << shaderParameter;

        shaderParameter = "mvMatrix";
        if ((m_mouseShader.mvMatrixLoc = m_mouseShader.program->uniformLocation(shaderParameter)) < 0)
            qDebug() << "Unable to find shader location for " << shaderParameter;

        shaderParameter = "Kd";
        if ((m_mouseShader.KdLoc = m_mouseShader.program->uniformLocation(shaderParameter)) < 0)
            qDebug() << "Unable to find shader location for " << shaderParameter;

        m_gl.glGenVertexArrays(1, &m_mouseVAO);
        m_gl.glGenBuffers(1, &m_mouseVBO);

        // Set VAO that binds the shader vertices inputs to the buffer data
        m_gl.glBindVertexArray(m_mouseVAO);
        m_gl.glBindBuffer(GL_ARRAY_BUFFER, m_mouseVBO);

        GLsizei dataSize = 6 * sizeof(GLfloat);

        // Initialize a buffer of versSize+normalsSize
        //
        m_gl.glBufferData(GL_ARRAY_BUFFER, dataSize, NULL, GL_DYNAMIC_DRAW);

        // Setup shader pointers to position and normal data
        //
        m_gl.glVertexAttribPointer(m_mouseShader.vPositionLoc, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));
        m_gl.glEnableVertexAttribArray(m_mouseShader.vPositionLoc);
    }
}

void SimViewer::draw()
{
    makeCurrent();
    glViewport( 0, 0, this->width(), this->height() );

    // Get projection and camera transformations
    QMatrix4x4 projectionMatrix;
    QMatrix4x4 modelViewMatrix;
    camera()->getProjectionMatrix(projectionMatrix);
    camera()->getModelViewMatrix(modelViewMatrix);

    QOpenGLFramebufferObject::bindDefault();
    m_rigidBodyRenderer->draw(projectionMatrix, modelViewMatrix);

    // Do particle picking.
    // Find the closest particle in the cloth and store its initial position.
    //
    if( m_pickFlag != kPickNone )
    {
        doPicking(camera(), QPoint(m_mouseX, m_mouseY), m_rigidBodySystem->getBodies(), m_pickingData);

        if( m_pickFlag == kPickFreeze && m_pickingData.body != nullptr )
        {
            m_pickingData.body->fixed = !m_pickingData.body->fixed;
            m_pickingData.body = nullptr;
        }

        m_pickFlag = kPickNone;
    }

    // Draw the mouse spring if a particle is selected.
    //
    if( m_pickingData.body != nullptr )
    {
        m_mouseShader.program->bind();
        m_mouseShader.program->setUniformValue(m_mouseShader.projMatrixLoc, projectionMatrix);
        m_mouseShader.program->setUniformValue(m_mouseShader.mvMatrixLoc, modelViewMatrix);

        // Set VAO that binds the shader vertices inputs to the buffer data
        m_gl.glBindVertexArray(m_mouseVAO);
        m_gl.glBindBuffer(GL_ARRAY_BUFFER, m_mouseVBO);

        // Copy vertex positions
        const Eigen::Vector3f bodyPos = m_pickingData.body->R * m_pickingData.plocal + m_pickingData.body->x;
        const GLsizei vertsSize = 6 * sizeof(GLfloat);
        const float verts[6] = { bodyPos[0], bodyPos[1], bodyPos[2],
                                 (float)m_pickingData.cursor[0], (float)m_pickingData.cursor[1], (float)m_pickingData.cursor[2] };
        m_gl.glBufferSubData(GL_ARRAY_BUFFER, 0, vertsSize, verts);
        // Colour.
        m_gl.glUniform3f(m_mouseShader.KdLoc, 0.1f, 1.0f, 0.1f);

        m_gl.glDisable(GL_DEPTH_TEST);
        m_gl.glDisable(GL_CULL_FACE);
        m_gl.glLineWidth(8.0f);
        m_gl.glDrawArrays(GL_LINES, 0, 2);
        m_gl.glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        m_gl.glEnable(GL_DEPTH_TEST);
        m_gl.glEnable(GL_CULL_FACE);
    }

    m_gl.glFlush();
}

void SimViewer::animate()
{
    if( !m_paused || m_stepOnce )
    {
        auto dynamics_start = std::chrono::high_resolution_clock::now();

        // Step the simulation.
        // The time step dt is divided by the number of sub-steps.
        //
        const float dt = m_dt / (float)m_subSteps;
        for(int i = 0; i < m_subSteps; ++i)
        {
            m_rigidBodySystem->step(dt);
        }
        auto dynamics_end = std::chrono::high_resolution_clock::now();
        m_dynamicsMicros = std::chrono::duration_cast<std::chrono::microseconds>(dynamics_end-dynamics_start);

        // Clear step-once flag.
        m_stepOnce = false;

        // Update status bar message.
        QString msg = "Frame: " + QString::number(m_rigidBodySystem->getFrameCount());
        const int numContacts = m_rigidBodySystem->getContacts().size();
        msg += ",  numContacts: " + QString::number(numContacts);
        msg += ", dynamics: " + QString::number(m_dynamicsMicros.count()) + " us";
        msg += ", render: " + QString::number(m_renderMicros.count()) + " us";
        emit statusMessageChanged(msg);
    }
}

void SimViewer::mouseMoveEvent(QMouseEvent* e)
{
    // Check if there is a picked body and update the cursor position.
    //
    if ((e->buttons() && Qt::LeftButton) && (e->modifiers() == Qt::ShiftModifier) && (m_pickingData.body != nullptr) )
    {
        qglviewer::Vec src;
        src[0] = e->x();
        src[1] = e->y();
        src[2] = camera()->zFar() / (camera()->zFar() - camera()->zNear()) * (1.0 - camera()->zNear() / m_pickingData.depth);
        m_pickingData.cursor = camera()->unprojectedCoordinatesOf(src);
    }

    // Normal QGLViewer behavior.
    QGLViewer::mouseMoveEvent(e);
}

void SimViewer::mousePressEvent(QMouseEvent *e)
{
    // Handle picking.
    //
    if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ShiftModifier) )
    {
        m_pickFlag = kPickMouseSpring;
        m_mouseX = e->x();
        m_mouseY = e->y();
    }
    else if ((e->button() == Qt::RightButton) && (e->modifiers() == Qt::ShiftModifier) )
    {
        m_pickFlag = kPickFreeze;
        m_mouseX = e->x();
        m_mouseY = e->y();
    }

	// Normal QGLViewer behavior
	QGLViewer::mousePressEvent(e);
}

void SimViewer::mouseReleaseEvent(QMouseEvent* e)
{
    // Release picked body
    if ( m_pickingData.body != nullptr && (e->button() == Qt::LeftButton) )
    {
        m_pickingData.body = nullptr;
    }

	QGLViewer::mouseReleaseEvent(e);
}

void SimViewer::setTimestep(double _dt)
{
    m_dt = _dt;
}

void SimViewer::setSubsteps(int _substeps)
{
    m_subSteps = _substeps;
}

void SimViewer::setContactStiffness(double _stiffness)
{
    m_rigidBodySystem->setContactStiffness(_stiffness);
}

void SimViewer::setContactDamping(double _damping)
{
    m_rigidBodySystem->setContactDamping(_damping);
}

void SimViewer::setPaused(bool _paused)
{
    m_paused = _paused;
}

void SimViewer::stepOnce()
{
    m_stepOnce = true;
}

void SimViewer::setDrawMesh(bool _draw)
{
    m_rigidBodyRenderer->setDrawBodiesEnabled(_draw);
}

void SimViewer::setDrawContacts(bool _draw)
{
    m_rigidBodyRenderer->setDrawContactsEnabled(_draw);
}

void SimViewer::createBoxOnPlane()
{
    Scenarios::createBoxOnPlane(*m_rigidBodySystem);
    m_rigidBodyRenderer->updateMeshVBOs();
}

void SimViewer::createBoxBallStack()
{
    Scenarios::createBoxBallStack(*m_rigidBodySystem);
    m_rigidBodyRenderer->updateMeshVBOs();
}

void SimViewer::createMarbleBox()
{
    Scenarios::createMarbleBox(*m_rigidBodySystem);
    m_rigidBodyRenderer->updateMeshVBOs();
}

void SimViewer::createBunnies()
{
    Scenarios::createBunnies(*m_rigidBodySystem);
    m_rigidBodyRenderer->updateMeshVBOs();
}

void SimViewer::setSolver(eSolverType _solverType)
{
    m_rigidBodySystem->setSolverType(_solverType);
}

void SimViewer::setMaxIterations(int _maxIter)
{
    m_rigidBodySystem->setSolverIterations(_maxIter);
}

void SimViewer::setFrictionCoefficient(double mu)
{
    m_rigidBodySystem->setFrictionCoefficient(mu);
}

void SimViewer::preStep(std::vector<RigidBody*>& _bodies)
{
    // Apply mouse spring forces for picked body
    //
    if( m_pickingData.body != nullptr )
    {
        const float k = 10.0f * m_pickingData.body->mass;
        const float b = 0.5f;
        const Eigen::Vector3f p = (m_pickingData.body->R * m_pickingData.plocal + m_pickingData.body->x);
        const Eigen::Vector3f cursorp = Eigen::Vector3f(m_pickingData.cursor[0], m_pickingData.cursor[1], m_pickingData.cursor[2]);
        const Eigen::Vector3f dx = cursorp - p;
        Eigen::Vector3f v(1, 0, 0);

        const float dx_nrm = dx.norm();
        if( dx_nrm > 1e-3f )
        {
            v = dx / dx_nrm;
        }

        Eigen::Vector3f vel;
        m_pickingData.body->getVelocityAtPos(p, vel);
        m_pickingData.body->addForceAtPos(p, k*dx-b*(vel.dot(v))*v);
    }
}
