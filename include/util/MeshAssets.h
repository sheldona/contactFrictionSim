#pragma once

#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QImage>
#include <string>
#include <vector>
#include <map>

// Data structure used to store a mesh material.
struct Material
{
    Material() : Kn(0.0f), name(""), img_diffuse(), img_bump(), map_diffuse(nullptr),  map_bump(nullptr)
    {
        memset(Ka, 0, sizeof(Ka));
        memset(Ke, 0, sizeof(Ke));
        memset(Kd, 0, sizeof(Kd));
        memset(Ks, 0, sizeof(Ks));
    }

    float Ka[4];  // Ambient color
    float Ke[4];  // Emissive color
    float Kd[4];  // Diffuse color
    float Ks[4];  // Specular color
    float Kn;     // Specular exponent

    std::string name; // Material's name
    QImage img_diffuse;    // Diffuse image
    QImage img_bump;       // Normal map image
    QOpenGLTexture* map_diffuse; // Diffuse texture
    QOpenGLTexture* map_bump; // Normal map texture
};


// Data structure used to store a vertex.
struct Vertex
{
    Vertex()
    {
        memset(position, 0, sizeof(position));
        memset(normal, 0, sizeof(normal));
        memset(uv, 0, sizeof(uv));
        memset(tangent, 0, sizeof(tangent));
    }

    GLfloat position[3];
    GLfloat normal[3];
    GLfloat tangent[3];
    GLfloat uv[2];
};

// Data structure for a polygon mesh.
struct Mesh
{
    Mesh() : name(""), vertices(), material(), vboOffset(0) { }

    std::string name;
    std::vector<Vertex> vertices;
    Material material;
    GLuint vboOffset;
};


// Cache of mesh data. Each filename is loaded only once.
typedef std::map<unsigned int, Mesh> MeshCache;

// The mesh registry is a static class storing geometry and material information
// for drawing rigid bodies.
//
class MeshAssetRegistry
{
public:

    // Loads a mesh from the OBJ file @a _filename
    // and adds it to the cache.  Returns a pointer to the mesh in the cache.
    static Mesh* loadObj(const std::string& _filename);

    // Clears the mesh cache.
    static void clear();

    // Returns the container of cached meshes.
    static MeshCache& cachedMeshes();

private:

    // The one and only instance of the mesh cache.
    static MeshCache m_meshCache;
};
