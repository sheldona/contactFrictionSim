#pragma once

/**
 * @file RigidBodyRenderer.h
 *
 * @brief Renderer for system of rigid bodies.
 *
 */

#include "ShaderVars.h"
#include <QVector3D>
#include <QVector4D>

QT_FORWARD_DECLARE_CLASS(QOpenGLFunctions_4_0_Core)
QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)
QT_FORWARD_DECLARE_CLASS(QOpenGLFramebufferObject)
QT_FORWARD_DECLARE_CLASS(QOpenGLTexture)
QT_FORWARD_DECLARE_CLASS(Contact)
QT_FORWARD_DECLARE_CLASS(RigidBodySystem)
QT_FORWARD_DECLARE_CLASS(RigidBody)

struct LightingParams
{
    LightingParams() :
        pos(0,0,0,1), diffuse(0.8f,0.8f,0.8f), specular(0.5f,0.5f,0.5f), ambient(0.1f, 0.1f, 0.1f)
    {

    }

    LightingParams(const QVector4D& _pos, const QVector3D& _diffuse, const QVector3D& _specular, const QVector3D& _ambient) :
        pos(_pos), diffuse(_diffuse), specular(_specular), ambient(_ambient)
    {

    }

    QVector4D pos;      // light position. If pos.w == 0, light is directional, point source otherwise.
    QVector3D diffuse;  // diffuse color
    QVector3D specular; // specular color
    QVector3D ambient;  // ambient color
};

class RigidBodyRenderer
{
public:
    explicit RigidBodyRenderer(QOpenGLFunctions_4_0_Core* _gl, RigidBodySystem* _system);

    virtual ~RigidBodyRenderer();

    virtual void draw(const QMatrix4x4& projectionMatrix, const QMatrix4x4& modelViewMatrix);

    void setDrawBodiesEnabled(bool);

    void setDrawContactsEnabled(bool);

    void updateMeshVBOs();

    void drawDepth(const QMatrix4x4& projectionMatrix, const QMatrix4x4& modelViewMatrix);

    void drawContactCoeffs(Contact* c, int bodyIndex);

    void drawContactWeights(Contact* c, int bodyIndex);

    void setLightingParameters(const LightingParams& lighting) { m_lighting = lighting; }

private:

    void init();
    void drawBodies(const QMatrix4x4& projectionMatrix, const QMatrix4x4& modelViewMatrix);
    void drawContacts(const QMatrix4x4& projectionMatrix, const QMatrix4x4& modelViewMatrix);
    void drawNormals(const QMatrix4x4& projectionMatrix, const QMatrix4x4& modelViewMatrix);

    RigidBodySystem* m_system;

    // GL objects
    GLuint m_meshVAO;
    GLuint m_meshVBO;
    GLuint m_contactVAO;
    GLuint m_contactVBO;

    LightingParams m_lighting;
    ShaderVars m_objectShader;
    ShaderVars m_contactShader;
    QOpenGLFunctions_4_0_Core* m_gl;

    bool m_drawBodies;
    bool m_drawContacts;
    std::vector<GLfloat> contactVerts;
};
