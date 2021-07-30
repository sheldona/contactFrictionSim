#pragma once

/**
 * @file SimViewer.h
 *
 * @brief Viewer for a cloth simulation application.
 *
 */

#include <QOpenGLFunctions_4_0_Core>
#include <qglviewer.h>
#include <QPoint>
#include <Eigen/Dense>
#include <memory>

#include "contact/Contact.h"
#include "ShaderVars.h"

QT_FORWARD_DECLARE_CLASS(RigidBodyRenderer)
QT_FORWARD_DECLARE_CLASS(RigidBodySystem)
QT_FORWARD_DECLARE_CLASS(QOpenGLFramebufferObject)

class SimViewer : public QGLViewer
{
    Q_OBJECT

public:
    SimViewer();
    virtual ~SimViewer();

    struct PickingData
    {
        PickingData() : depth(std::numeric_limits<float>::max()), plocal(), body(nullptr) {}
        float depth;
        Eigen::Vector3f plocal;
        RigidBody* body;
        qglviewer::Vec cursor;
    };

public slots:
    void cleanup();
    void setDrawMesh(bool);
    void setDrawContacts(bool);
    void setTimestep(double);
    void setSubsteps(int);
    void setMaxIterations(int);
    void setFrictionCoefficient(double);
    void setContactStiffness(double);
    void setContactDamping(double);
    void setPaused(bool);
    void stepOnce();

    void createBoxOnPlane();
    void createBoxBallStack();
    void createMarbleBox();
    void createBunnies();

signals:

    void statusMessageChanged(const QString&);

protected :
    virtual void animate() override;
    virtual void draw() override;
    virtual void init() override;

    virtual void mouseMoveEvent(QMouseEvent* e) override;
    virtual void mousePressEvent(QMouseEvent *e) override;
    virtual void mouseReleaseEvent(QMouseEvent* e) override;

    QOpenGLFunctions_4_0_Core m_gl;

    void preStep(std::vector<RigidBody*>&);
    void onReset();

private:

    // Simulation parameters
    float m_dt;                         //< Time step parameter.
    int m_subSteps;
    bool m_paused;                      //< Pause the simulation.
    bool m_stepOnce;                    //< Advance the simulation by one frame and then stop.

    // Shader for drawing the mouse spring.
    GLuint m_mouseVAO;
    GLuint m_mouseVBO;
    ShaderVars m_mouseShader;

    // Picking.
    enum ePickType { kPickNone = 0, kPickMouseSpring, kPickFreeze } ;
    ePickType m_pickFlag;               //< Flag indicating if the user tried to pick a particle.
    int m_mouseX, m_mouseY;             //< Mouse x,y position on the screen (in pixels).
    PickingData m_pickingData;          //< Mouse picking data.

    std::chrono::microseconds m_renderMicros;
    std::chrono::microseconds m_dynamicsMicros;

    std::unique_ptr<RigidBodySystem> m_rigidBodySystem;
    std::unique_ptr<RigidBodyRenderer> m_rigidBodyRenderer;
};
