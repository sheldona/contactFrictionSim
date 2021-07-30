#version 400 core
uniform mat4 mvMatrix;
uniform mat4 projMatrix;
in vec4 vPosition;

void
main()
{
     gl_Position = projMatrix * mvMatrix * vPosition;
}

