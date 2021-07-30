#version 400 core
uniform vec3 Kd;
out vec4 fColor;

void
main()
{
    fColor = vec4(Kd, 1);
}
