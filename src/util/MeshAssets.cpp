#include "util/MeshAssets.h"
#include "util/OBJLoader.h"
#include <cassert>
#include <functional>

std::map<unsigned int, Mesh> MeshAssetRegistry::m_meshCache = std::map<unsigned int, Mesh>();

Mesh* MeshAssetRegistry::loadObj(const std::string& _filename)
{
    unsigned int key = std::hash<std::string>{}(_filename);
    auto cacheItr = m_meshCache.find(key);

    if( cacheItr != m_meshCache.end() )
    {
        return &(cacheItr->second);
    }

    OBJLoader obj(_filename);

    const auto& meshes = obj.getMeshes();
    assert(meshes.size() > 0);

    const auto& materials = obj.getMaterials();
    assert(materials.size() > 0);

    Mesh mesh = meshes[0];
    mesh.material = (materials.size() < 2) ? materials[0] : materials[1]; // 0 = default material, 1 = first material in the file
    m_meshCache[key] = mesh;
    return &(m_meshCache[key]);
}

void MeshAssetRegistry::clear()
{
    m_meshCache.clear();
}

MeshCache& MeshAssetRegistry::cachedMeshes()
{
    return m_meshCache;
}


