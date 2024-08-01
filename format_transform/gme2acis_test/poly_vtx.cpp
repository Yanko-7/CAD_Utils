#include "acis/include/poly_vtx.hxx"

polygon_vertex::polygon_vertex(const char* gme) {
    m_Position = SPAposition();
    m_Normal = SPAunit_vector();
    m_SurfaceUV = SPApar_pos();
    m_Color[0] = m_Color[1] = m_Color[2] = -1.0;
}

polygon_vertex::polygon_vertex(const char* gme, const SPAposition& pos, const SPAunit_vector& uv, const SPApar_pos& parpos) {
    m_Position = pos;
    m_Normal = uv;
    m_SurfaceUV = parpos;
    m_Color[0] = m_Color[1] = m_Color[2] = -1.0;
}

polygon_vertex::polygon_vertex(const char* gme, const polygon_vertex& polver) {
    m_Position = polver.m_Position;
    m_Normal = polver.m_Normal;
    m_SurfaceUV = polver.m_SurfaceUV;
    m_Color[0] = polver.m_Color[0];
    m_Color[1] = polver.m_Color[1];
    m_Color[2] = polver.m_Color[2];
}

polygon_vertex& polygon_vertex::assign(const polygon_vertex& other) {
    if(this != &other) {
        m_Position = other.m_Position;
        m_Normal = other.m_Normal;
        m_SurfaceUV = other.m_SurfaceUV;
        m_Color[0] = other.m_Color[0];
        m_Color[1] = other.m_Color[1];
        m_Color[2] = other.m_Color[2];
    }
    return *this;
}

const SPAposition& polygon_vertex::get_position(const char* gme) {
    return m_Position;
}

const SPAunit_vector& polygon_vertex::get_normal(const char* gme) {
    return m_Normal;
}

const SPApar_pos& polygon_vertex::get_uv(const char* gme) {
    return m_SurfaceUV;
}

const double* polygon_vertex::get_color(const char* gme) {
    return m_Color;
}

void polygon_vertex::set_position(const char* gme, const SPAposition& pos) {
    m_Position = pos;
}

void polygon_vertex::set_normal(const char* gme, const SPAunit_vector& norm) {
    m_Normal = norm;
}

void polygon_vertex::set_uv(const char* gmee, const SPApar_pos& uv) {
    m_SurfaceUV = uv;
}

void polygon_vertex::set_uv(const char* gme, double u, double v) {
    m_SurfaceUV.u = u;
    m_SurfaceUV.v = v;
}

void polygon_vertex::set_color(const char* gme, const double* color) {
    if(color != nullptr) {
        m_Color[0] = color[0];
        m_Color[1] = color[1];
        m_Color[2] = color[2];
    }
}

void polygon_vertex::set_data(const char* gme, const SPAposition& pos, const SPAunit_vector& uv, const SPApar_pos& parpos) {
    m_Position = pos;
    m_Normal = uv;
    m_SurfaceUV = parpos;
}

void polygon_vertex::operator*=(const SPAtransf& t) {
    m_Position *= t;
    m_Normal *= t;
}