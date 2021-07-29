#pragma once

#include <QOpenGLShaderProgram>

// A simple data structure to store the shader program
// and commonly used shader parameters.
//
struct ShaderVars
{
    QOpenGLShaderProgram* program;
    GLuint projMatrixLoc;
    GLuint mvMatrixLoc;
    GLuint normalMatrixLoc;
    GLuint vPositionLoc;
    GLuint vNormalLoc;
    GLuint vTangentLoc;
    GLuint lPositionLoc;
    GLuint vUvLoc;
    GLuint KdLoc;
    GLuint KsLoc;
    GLuint KnLoc;
    GLuint lKdLoc;
    GLuint lKsLoc;
    GLuint lKaLoc;
    GLuint useTextureLoc;
    GLuint useNormalMapLoc;
    GLuint useLightingLoc;
    GLuint userParams[6];
};
