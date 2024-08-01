#include "acis/include/ptlist.hxx"

#include <iostream>

#include "acis/gme/faceter/gme_facet_utils.hxx"

AF_POINT::AF_POINT(const char* gme, AF_POINT_ID _id, AF_POINT* prev, int sense) {
    id = _id;
    /**
     * sense = 0 = FORWARD, default;
     * sense = 1 = REVERSED.
     */

    if(!sense) {  // FORWARD
        if(prev == (AF_POINT*)nullptr) {
            /**
             * 首节点：加入时prev为空，则前后指针均指向自己
             */
            this->next_point = this;
            this->prev_point = this;
        } else {
            // the 't' will be the head AF_POINT of the ordered AF_POINT linked list
            AF_POINT* t = prev->prev_point;
            // std::cout << t->get_position().x() << " " << t->get_position().y() << " " << t->get_position().z() << std::endl;
            while(t != prev) {
                t = t->prev_point;
            }
            t = t->next_point;  // now this 't' is actually the head AF_POINT
            this->next_point = t;
            t->prev_point = this;
            prev->next_point = this;
            this->prev_point = prev;
#if FACETER_DEBUG_MODE
            std::cout << t->get_position().x() << " " << t->get_position().y() << " " << t->get_position().z() << std::endl;
            std::cout << t->next_point->get_position().x() << " " << t->next_point->get_position().y() << " " << t->next_point->get_position().z() << std::endl;
            std::cout << t->prev_point->get_position().x() << " " << t->prev_point->get_position().y() << " " << t->prev_point->get_position().z() << std::endl;
#endif
        }
    } else {  // REVERSED
        if(prev == (AF_POINT*)nullptr) {
            this->prev_point = this;
            this->next_point = this;
        } else {
            AF_POINT* t = prev;
            while(t != prev) {
                t = prev->next_point;
            }
            t = t->prev_point;
            this->prev_point = t;
            prev->prev_point = this;
            t->next_point = this;
            this->next_point = prev;
        }
    }
}

AF_POINT* AF_POINT::next(const char* gme, int sense) const {
    if(sense) {
        return prev_point;
    } else {
        return next_point;
    }
}

void AF_POINT::set_next(AF_POINT* nxt) {
    this->next_point = nxt;
}

void AF_POINT::set_prev(AF_POINT* pre) {
    this->prev_point = pre;
}

void AF_POINT::insert(AF_POINT* prev, int sense) {
    AF_POINT* next = prev->next(0);

    prev->set_next(this);
    next->set_prev(this);
    this->set_prev(prev);
    this->set_next(next);
}

/** static */
logical AF_POINT::find(const char* gme, ENTITY* e, int sense, AF_POINT*& P0, AF_POINT*& P1) {
    AF_POINT *head = nullptr, *tail = nullptr;
    logical res = FALSE;
    ATTRIB_EYE_POINTLIST_HEADER* attr = ATTRIB_EYE_POINTLIST_HEADER::find( e);
    if(attr != nullptr && attr->get_flag( AF_POINTLIST_DIRTY)) {
        head = attr->get_pointlist("gme");
        tail = head;
        while(tail->next( FORWARD) != head) {
            tail = tail->next( FORWARD);
        }
        if(sense) {
            P1 = head;
            P0 = tail;
        } else {
            P0 = head;
            P1 = tail;
        }
        res = TRUE;
    }
    return res;
}

void AF_POINT::attach(const char* gme, ENTITY* e) {
    ATTRIB_EYE_POINTLIST_HEADER* attr = new ATTRIB_EYE_POINTLIST_HEADER( e);
    AF_POINT* p = this;
    while(p->get_user_id() != 0) {
        // 找到id为0的首节点
        p = p->next( REVERSED);
    }
    attr->replace_pointlist( p);
}

void AF_POINT::set_position(const char* gme,
                            const SPAposition& Xtemp  // Cartesian coordinates
) {
    X = Xtemp;
}

void AF_POINT::set_parameter(const char* gme,
                             const double& _t  // Parametric coordinate
) {
    t = _t;
}

// AF_POINT_LIST::~AF_POINT_LIST() {
//     m_pAfPoint = nullptr;
// }

AF_POINT_LIST::AF_POINT_LIST(const char* gme, AF_POINT* pAfPoint): m_pAfPoint(pAfPoint), m_useCount(0) {
    AF_POINT* p = m_pAfPoint;
    AF_POINT* h = m_pAfPoint;
    while(m_useCount != 0 && p != h) {
        m_useCount++;
        p = p->next( FORWARD);
    }
}

std::unordered_map<::ENTITY*, ATTRIB_EYE_POINTLIST_HEADER*> ATTRIB_EYE_POINTLIST_HEADER::ep_set = {};

ATTRIB_EYE_POINTLIST_HEADER::ATTRIB_EYE_POINTLIST_HEADER(const char* gme, ENTITY* pent) {
    if(find( pent) != nullptr) {
        ep_set.erase(pent);
    }
    ent = pent;
    ep_set.emplace(ent, this);
    flags = 0;
    m_pPointList = nullptr;
}

// ATTRIB_EYE_POINTLIST_HEADER::~ATTRIB_EYE_POINTLIST_HEADER() {
//     delete m_pPointList;
//     m_pPointList = nullptr;
//     ep_set.erase(ent);
// }

AF_POINT* ATTRIB_EYE_POINTLIST_HEADER::get_pointlist(const char* gme) {
    return m_pPointList ? m_pPointList->GetPoint() : nullptr;
}

void ATTRIB_EYE_POINTLIST_HEADER::replace_pointlist(const char* gme, AF_POINT* pNewPoint) {
    if(m_pPointList != nullptr) delete m_pPointList;
    m_pPointList = new AF_POINT_LIST( pNewPoint);
    int v = pNewPoint == nullptr ? 0 : 1;
    set_flag( AF_POINTLIST_DIRTY, v);
}

/** static */
ATTRIB_EYE_POINTLIST_HEADER* ATTRIB_EYE_POINTLIST_HEADER::find(const char* gme, ENTITY* ent) {
    if(ent != nullptr) {
        if(ep_set.find(ent) != ep_set.end()) return ep_set.at(ent);
    }
    return nullptr;
}

unsigned int ATTRIB_EYE_POINTLIST_HEADER::get_flag(const char* gme, unsigned long mask) {
    return flags & mask;
}

void ATTRIB_EYE_POINTLIST_HEADER::set_flag(const char* gme, unsigned long mask, int value) {
    if(value) {
        flags |= mask;
    } else {
        flags &= ~mask;
    }
}