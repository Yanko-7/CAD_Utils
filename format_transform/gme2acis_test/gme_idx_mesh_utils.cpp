#include "acis/gme/faceter/gme_idx_mesh_utils.hxx"

INDEXED_MESH* GetIndexedMesh(FACE* f) {
    if(fm_set.find(f) != fm_set.end())
        return fm_set.at(f);
    else
        return nullptr;
}

void attach_indexed_mesh_to_face(FACE* face, INDEXED_MESH* mesh) {
    INDEXED_MESH* old_mesh = GetIndexedMesh(face);
    if(old_mesh != NULL) {
        delete old_mesh;
        fm_set.erase(face);
    }
    fm_set.emplace(face, mesh);
}
