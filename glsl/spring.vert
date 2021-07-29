#version 430 core
layout (location=0) uniform mat4 mvMatrix;
layout (location=1) uniform mat4 projMatrix;
layout (location=3) in vec4 vPosition;

void
main()
{
     gl_Position = projMatrix * mvMatrix * vPosition;
}

