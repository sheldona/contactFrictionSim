#pragma once

/**
 * @file OBJLoader.h
 *
 * @brief File loading for OBJ meshes.
 *
 */
#include "MeshAssets.h"

class QOpenGLTexture;

// Class responsible for loading all the meshes included in an OBJ file
class OBJLoader
{
public:
    OBJLoader();
    OBJLoader(const std::string& filename);
    ~OBJLoader();

    bool loadFile(const std::string& filename);
    bool isLoaded() const { return _isLoaded; }
    void unload();

    const std::vector<Mesh>& getMeshes() const { return _meshes; }
    const std::vector<Material>& getMaterials() const { return _materials; }

private:
    void loadMtlFile(const std::string& filename);
    unsigned int findMaterial(const std::string& name);
    unsigned int getMesh(const std::string& name);

    std::vector<Mesh>     _meshes;
    std::vector<Material> _materials;

    bool                  _isLoaded;
};
