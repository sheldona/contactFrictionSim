#version 430 core
layout(binding=0) uniform sampler2D mapKd;
layout(binding=1) uniform sampler2D mapNormal;
uniform vec3 Kd;
uniform vec3 Ks;
uniform float Kn;
uniform vec4 lPosition;
uniform vec3 lKd;
uniform vec3 lKs;
uniform vec3 lKa;
uniform int useTexture;
uniform int useNormalMap;
uniform int useLighting;
in vec3 fPosition;
in vec3 fNormal;
in vec3 fTangent;
in vec3 fBinormal;
in vec2 fUV;
in vec4 fEye;
out vec4 fColor;

void
main()
{
    vec3 diffuseColor = Kd;
    vec3 specularColor = Kd;
    vec3 lightDiffuseColor = lKd;
    vec3 lightSpecularColor = lKs;
    vec3 lightAmbientColor = lKa;
    vec3 nfNormal = normalize(fNormal);

    if( useTexture > 0)
    {
        diffuseColor = texture(mapKd, fUV).rgb;
    }

    if( useNormalMap > 0 )
    {
        // Build the matrix to transform from XYZ (normal map) space to TBN (tangent) space
        // Each vector fills a column of the matrix
        mat3 tbn = mat3(normalize(fTangent), normalize(fBinormal), normalize(fNormal));
        vec3 normalFromTexture = texture(mapNormal, fUV).rgb * 2.0 - vec3(1.0);
        nfNormal = normalize(tbn * normalFromTexture);
    }

    if( useLighting > 0 )
    {
        // Get lighting vectors
        vec3 l = (lPosition.w == 0) ? normalize(lPosition.xyz) : normalize(lPosition.xyz-fPosition);
        vec3 e = normalize(fEye.xyz);
        vec3 h = normalize(l+e);


        // Compute diffuse component
        float intensity = max(0.0, dot(nfNormal, l));
        vec3 diffuse = diffuseColor * lightDiffuseColor * intensity;

        // Compute specular component
        float spec = max(0.0, dot(h, nfNormal));
        vec3 specular = specularColor * lightSpecularColor * pow(spec, Kn);

        // Ambient component
        vec3 ambient = lightAmbientColor * Kd;

        // Compute final color
        fColor = vec4(ambient + diffuse +  specular, 1);
    }
    else
    {
        fColor = vec4(diffuseColor, 1);
    }
}
