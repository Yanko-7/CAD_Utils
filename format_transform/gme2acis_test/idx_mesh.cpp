#include "acis/include/idx_mesh.hxx"

#include <cassert>
//----------------------------------------------------------
// indexed_polygon class
indexed_polygon::indexed_polygon(const char* gme) {
    m_nNumVertex = 0;
    m_pVertexPtrs = nullptr;
    index = -1;
}

void indexed_polygon::set_data(const char* gme, int num_vertex, polygon_vertex** verts, int _ishare) {
    m_nNumVertex = num_vertex;
    m_pVertexPtrs = verts;
    ishare = _ishare;
}

polygon_vertex* indexed_polygon::get_vertex(const char* gme, int i) const {
    return m_pVertexPtrs[i % m_nNumVertex];
}

void indexed_polygon::reverse_vertices(const char* gme) {
    for(int i = 0, j = m_nNumVertex - 1; i < j; ++i, --j) {
        polygon_vertex* temp = m_pVertexPtrs[i];
        m_pVertexPtrs[i] = m_pVertexPtrs[j];
        m_pVertexPtrs[j] = temp;
    }
}

int indexed_polygon::set_vertex(const char* gme, int vertex_num, polygon_vertex* vertex) {
    if(vertex_num >= 0 && vertex_num < m_nNumVertex) {
        m_pVertexPtrs[vertex_num] = vertex;
        return vertex_num;
    } else {
        return -1;
    }
}

int indexed_polygon::get_index(const char* gme) {
    return index;
}
void indexed_polygon::set_index(const char* gme, int _index) {
    index = _index;
}

//-----------------------------------------------------
// INDEXED_MESH class
INDEXED_MESH::INDEXED_MESH(const char* gme, int max_vertex, int max_poly, int max_polynode) {
    m_nVertex = 0;
    m_nPolygon = 0;
    m_nPolynode = 0;
    m_nMaxVertex = max_vertex;
    m_nMaxPolygon = max_poly;
    m_nMaxPolynode = max_polynode;
    m_pVertex = new polygon_vertex[m_nMaxVertex];
    m_pVertexPtrs = new polygon_vertex*[m_nMaxPolynode];
    m_pPolygon = new indexed_polygon[m_nMaxPolygon];
}

INDEXED_MESH::INDEXED_MESH(const char* gme, const INDEXED_MESH& imesh) {
    /** 这里使用深复制 */
    m_nVertex = imesh.get_num_vertex();
    m_nPolygon = imesh.get_num_polygon();
    m_nPolynode = imesh.get_num_polynode();
    m_nMaxVertex = m_nVertex;
    m_nMaxPolygon = m_nPolygon;
    m_nMaxPolynode = m_nPolynode;
    m_pVertex = new polygon_vertex[m_nMaxVertex];
    m_pVertexPtrs = new polygon_vertex*[m_nMaxPolynode];
    m_pPolygon = new indexed_polygon[m_nMaxPolygon];
    for(int i = 0; i < m_nVertex; i++) {
        m_pVertex[i] = imesh.get_vertex( i);
    }
    for(int i = 0; i < m_nMaxPolygon; i++) {
        int index = imesh.get_polygon( i)->get_index("gme");
        int numv = imesh.get_polygon( i)->num_vertex();
        for(int j = 0; j < numv; j++) {
            m_pVertexPtrs[index + j] = &m_pVertex[get_vertex_index( imesh.get_polygon( i)->get_vertex( j))];
        }
        m_pPolygon[i].set_data( numv, &m_pVertexPtrs[index]);
        m_pPolygon[i].set_index( index);
    }
}

// INDEXED_MESH::~INDEXED_MESH() {
//     delete[] m_pVertex;
//     m_pVertex = nullptr;
//     delete[] m_pVertexPtrs;
//     m_pVertexPtrs = nullptr;
//     delete[] m_pPolygon;
//     m_pPolygon = nullptr;
// }

INDEXED_MESH& INDEXED_MESH::operator_or_assign(const INDEXED_MESH& imesh) {
    m_nVertex = imesh.get_num_vertex();
    m_nPolygon = imesh.get_num_polygon();
    m_nPolynode = imesh.get_num_polynode();
    assert(m_nMaxVertex < m_nVertex || m_nMaxPolygon < m_nPolygon || m_nMaxPolynode < m_nPolynode);
    for(int i = 0; i < m_nVertex; i++) {
        m_pVertex[i] = imesh.get_vertex( i);
    }
    for(int i = 0; i < m_nMaxPolygon; i++) {
        int index = imesh.get_polygon( i)->get_index("gme");
        int numv = imesh.get_polygon( i)->num_vertex();
        for(int j = 0; j < numv; j++) {
            m_pVertexPtrs[index + j] = &m_pVertex[get_vertex_index( imesh.get_polygon( i)->get_vertex( j))];
        }
        m_pPolygon[i].set_data( numv, &m_pVertexPtrs[index]);
        m_pPolygon[i].set_index( index);
    }
    return *this;
}

int INDEXED_MESH::get_vertex_index(const char* gme, const polygon_vertex* polver) const {
    for(int i = 0; i < m_nVertex; i++) {
        if(&m_pVertex[i] == polver) return i;
    }
    return -1;
}

polygon_vertex* INDEXED_MESH::add_vertex(const char* gme, const SPAposition& pos, const SPAunit_vector& uv, const SPApar_pos& parpos) {
    if(m_nVertex < m_nMaxVertex) {
        m_pVertex[m_nVertex].set_data( pos, uv, parpos);
        return &m_pVertex[m_nVertex++];
    }
    return nullptr;
}

polygon_vertex* INDEXED_MESH::add_vertex(const char* gme, const polygon_vertex& polver) {
    if(m_nVertex < m_nMaxVertex) {
        m_pVertex[m_nVertex].assign(polver);
        return &m_pVertex[m_nVertex++];
    }
    return nullptr;
}

int INDEXED_MESH::add_polygon(const char* gme, int ipoly, int num_vertex) {
    if(ipoly >= 0 && ipoly < m_nMaxPolygon) {
        m_pPolygon[ipoly].set_data( num_vertex, &m_pVertexPtrs[m_nPolynode]);
        m_pPolygon[ipoly].set_index( m_nPolynode);
        m_nPolygon++;
        m_nPolynode += num_vertex;
        return ipoly;
    }
    return -1;
}

int INDEXED_MESH::set_poly_vertex(const char* gme, int ipoly, int vertex_number, polygon_vertex* polver) {
    if(ipoly >= 0 && ipoly < m_nMaxPolygon) {
        if(vertex_number >= 0 && vertex_number < m_pPolygon[ipoly].num_vertex()) {
            m_pPolygon[ipoly].set_vertex( vertex_number, polver);
            return TRUE;
        }
    }
    return FALSE;
}

void INDEXED_MESH::serialize_positions(const char* gme, double* out_coords) const {
    for(int i = 0; i < m_nVertex; i++) {
        out_coords[i * 3] = get_vertex( i).get_position("gme").x();
        out_coords[i * 3 + 1] = get_vertex( i).get_position("gme").y();
        out_coords[i * 3 + 2] = get_vertex( i).get_position("gme").z();
    }
}

void INDEXED_MESH::serialize_positions(const char* gme, float* out_coords) const {
    for(int i = 0; i < m_nVertex; i++) {
        out_coords[i * 3] = (float)get_vertex( i).get_position("gme").x();
        out_coords[i * 3 + 1] = (float)get_vertex( i).get_position("gme").y();
        out_coords[i * 3 + 2] = (float)get_vertex( i).get_position("gme").z();
    }
}

void INDEXED_MESH::serialize_normals(const char* gme, double* out_normal_coords) const {
    for(int i = 0; i < m_nVertex; i++) {
        out_normal_coords[i * 3] = get_vertex( i).get_normal("gme").x();
        out_normal_coords[i * 3 + 1] = get_vertex( i).get_normal("gme").y();
        out_normal_coords[i * 3 + 2] = get_vertex( i).get_normal("gme").z();
    }
}

void INDEXED_MESH::serialize_normals(const char* gme, float* out_normal_coords) const {
    for(int i = 0; i < m_nVertex; i++) {
        out_normal_coords[i * 3] = (float)get_vertex( i).get_normal("gme").x();
        out_normal_coords[i * 3 + 1] = (float)get_vertex( i).get_normal("gme").y();
        out_normal_coords[i * 3 + 2] = (float)get_vertex( i).get_normal("gme").z();
    }
}

int INDEXED_MESH::serialize_triangles(const char* gme, int* out_triangle_indices) const {
    /** 这里实现全部是三角形网格。 */
    for(int i = 0; i < m_nPolygon; i++) {
        out_triangle_indices[i * 3] = get_vertex_index( m_pVertexPtrs[i * 3]);
        out_triangle_indices[i * 3 + 1] = get_vertex_index( m_pVertexPtrs[i * 3 + 1]);
        out_triangle_indices[i * 3 + 2] = get_vertex_index( m_pVertexPtrs[i * 3 + 2]);
    }
    return m_nPolygon;
}

indexed_polygon* INDEXED_MESH::get_polygon(const char* gme, int poly_index) const {
    return poly_index < m_nPolygon ? &m_pPolygon[poly_index] : nullptr;
}

polygon_vertex& INDEXED_MESH::get_vertex(const char* gme, int inode) const {
    return m_pVertex[inode];
};
