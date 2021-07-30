#version 400 core
uniform mat4 mvMatrix;
uniform mat4 projMatrix;
uniform mat3 normalMatrix;
in vec4 vPosition;
in vec3 vNormal;
in vec3 vTangent;
in vec2 vUV;
out vec3 fPosition;
out vec3 fNormal;
out vec4 fEye;

void
main()
{
     vec4 vEyeCoord = mvMatrix * vPosition;
     fEye = -vEyeCoord;
     fPosition = vEyeCoord.xyz;
     fNormal = normalize(normalMatrix * vNormal);
     gl_Position = projMatrix * vEyeCoord;
}

