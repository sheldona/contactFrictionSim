#version 430 core
uniform mat4 mvMatrix;
uniform mat4 projMatrix;
uniform mat3 normalMatrix;
in vec4 vPosition;
in vec3 vNormal;
in vec3 vTangent;
in vec2 vUV;
out vec3 fPosition;
out vec3 fNormal;
out vec2 fUV;
out vec4 fEye;
out vec3 fTangent;
out vec3 fBinormal;

void
main()
{
     vec4 vEyeCoord = mvMatrix * vPosition;
     fEye = -vEyeCoord;
     fPosition = vEyeCoord.xyz;
     fNormal = normalize(normalMatrix * vNormal);
     fTangent = mat3(mvMatrix) * vTangent;
     fBinormal = mat3(mvMatrix) * cross(vNormal, vTangent);
     fUV = vUV;
     gl_Position = projMatrix * vEyeCoord;
}

