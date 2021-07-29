#version 430 core
layout(location=2) uniform vec3 Kd;
out vec4 fColor;

void
main()
{
    fColor = vec4(Kd, 1);
}
