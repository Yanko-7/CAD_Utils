#include "acis/gme/faceter/gme_facet_face_utils.hxx"

#include <vector>

#include "acis/include/intcucu.hxx"
#include "acis/include/intrapi.hxx"
#include "acis/include/position.hxx"
#include "acis/include/ptlist.hxx"

FaceterLoopInfo CheckLoopType(const point_vector& loop, const FaceterFaceInfo& info) {
    // 输入的loop要求最后一个点与第一个点重复
    // 确保输入的loop经过了预处理，是比较正常的，不会在边界处反复横条的loop
    // 上穿下穿次数，其正负暗含了uloop的指向
    int up_cnt = 0;
    // 左穿右穿次数
    int le_cnt = 0;
    bool is_cross_u = false;
    bool is_cross_v = false;
    point curr, next;
    for(size_t i = 0; i < loop.size() - 1; ++i) {
        curr = loop[i];
        next = loop[i + 1];
        CrossType cross_u = CrossBoundary(curr.x, next.x, info.u_tolerance);
        CrossType cross_v = CrossBoundary(curr.y, next.y, info.v_tolerance);
        if(cross_u != CrossType::kNoCross && !info.is_periodic_u || cross_v != CrossType::kNoCross && !info.is_periodic_v) return FaceterLoopInfo(FaceterLoopType::kError, FaceterLoopCrossType::kNoCross);

        switch(cross_u) {
            case CrossType::kCrossLower:
                le_cnt -= 1;
                is_cross_u = true;
                break;
            case CrossType::kCrossUpper:
                le_cnt += 1;
                is_cross_u = true;
                break;
            case CrossType::kNoCross:
                break;
        }
        switch(cross_v) {
            case CrossType::kCrossLower:
                up_cnt -= 1;
                is_cross_v = true;
                break;
            case CrossType::kCrossUpper:
                up_cnt += 1;
                is_cross_v = true;
                break;
            case CrossType::kNoCross:
                break;
        }
    }
    FaceterLoopInfo return_info(FaceterLoopType::kHole, FaceterLoopCrossType::kNoCross);
    if(is_cross_u && is_cross_v) return_info.loop_cross_type = FaceterLoopCrossType::kCrossUV;
    if(is_cross_u) return_info.loop_cross_type = FaceterLoopCrossType::kCrossU;
    if(is_cross_v) return_info.loop_cross_type = FaceterLoopCrossType::kCrossV;

    // 应当是不为偶数
    if(up_cnt != 0 && le_cnt != 0) return_info.loop_type = FaceterLoopType::kError;
    if(up_cnt != 0) {
        return_info.loop_type = FaceterLoopType::kULoop;
        if(up_cnt == 1) return_info.is_up = true;
    }
    if(le_cnt != 0) {
        return_info.loop_type = FaceterLoopType::kVLoop;
        if(le_cnt == 1) return_info.is_up = true;
    }
    return return_info;
}

void PeriodicPointPreprocess(point_vector& loop, const FaceterFaceInfo& info) {
    point curr, next;
    for(size_t i = 0; i < loop.size() - 1; i++) {
        curr = loop[i];
        next = loop[i + 1];
        // if(curr == next) {
        if(curr.x == next.x && curr.y == next.y) {
            loop.erase(loop.begin() + i + 1);
            i--;
            continue;
        }
        PointBoundaryType curr_t = GetPointBoundaryType(curr, info);
        PointBoundaryType next_t = GetPointBoundaryType(next, info);

        if(info.is_periodic_u) {
            CrossType cross_u = CrossBoundary(curr.x, next.x, info.u_tolerance);
            bool curr_on_u = curr_t == PointBoundaryType::kUle || curr_t == PointBoundaryType::kUri || curr_t == PointBoundaryType::klele || curr_t == PointBoundaryType::kleri || curr_t == PointBoundaryType::kriri || curr_t == PointBoundaryType::krile;
            bool next_on_u = next_t == PointBoundaryType::kUle || next_t == PointBoundaryType::kUri || next_t == PointBoundaryType::klele || next_t == PointBoundaryType::kleri || next_t == PointBoundaryType::kriri || next_t == PointBoundaryType::krile;
            // curr不处于U边界，而next处于U边界，则检查next的u值是否正确
            if(!curr_on_u && next_on_u) {
                if(cross_u == CrossType::kCrossUpper) {
                    loop[i + 1].x = info.u_ri;
                } else if(cross_u == CrossType::kCrossLower) {
                    loop[i + 1].x = info.u_le;
                }
            }  // 两者均位于边界
            else if(curr_on_u && next_on_u) {
                loop[i + 1].x = loop[i].x;
                if(equal(curr.y, next.y)) {
                    loop.erase(loop.begin() + i + 1);
                    i--;
                    continue;
                }
            }  // curr处于U边界，next不处于U边界，视情况添加周期点
            else if(curr_on_u && !next_on_u) {
                if(cross_u == CrossType::kCrossUpper) {
                    loop.insert(loop.begin() + i + 1, point(info.u_le, curr.y));
                    continue;
                } else if(cross_u == CrossType::kCrossLower) {
                    loop.insert(loop.begin() + i + 1, point(info.u_ri, curr.y));
                    continue;
                }
            }
        }
        if(info.is_periodic_v) {
            CrossType cross_v = CrossBoundary(curr.y, next.y, info.v_tolerance);
            bool curr_on_v = curr_t == PointBoundaryType::kVle || curr_t == PointBoundaryType::kVri || curr_t == PointBoundaryType::klele || curr_t == PointBoundaryType::krile || curr_t == PointBoundaryType::kriri || curr_t == PointBoundaryType::kleri;
            bool next_on_v = next_t == PointBoundaryType::kVle || next_t == PointBoundaryType::kVri || next_t == PointBoundaryType::klele || next_t == PointBoundaryType::krile || next_t == PointBoundaryType::kriri || next_t == PointBoundaryType::kleri;
            // curr不处于V边界，而next处于V边界，则检查next的v值是否正确
            if(!curr_on_v && next_on_v) {
                if(cross_v == CrossType::kCrossUpper) {
                    loop[i + 1].y = info.v_ri;
                } else if(cross_v == CrossType::kCrossLower) {
                    loop[i + 1].y = info.v_le;
                }
            }  // 两者均位于边界
            else if(curr_on_v && next_on_v) {
                loop[i + 1].y = curr.y;
                if(equal(curr.x, next.x)) {
                    loop.erase(loop.begin() + i + 1);
                    i--;
                    continue;
                }
            }  // curr处于V边界，next不处于V边界，视情况添加周期点
            else if(curr_on_v && !next_on_v) {
                if(cross_v == CrossType::kCrossUpper) {
                    loop.insert(loop.begin() + i + 1, point(curr.x, info.v_le));
                } else if(cross_v == CrossType::kCrossLower) {
                    loop.insert(loop.begin() + i + 1, point(curr.x, info.v_ri));
                }
            }
        }
    }
}

void FPLoopPreprocessor(point_vector& loop, const FaceterFaceInfo& info) {
    if(*loop.begin() != *loop.rbegin()) loop.push_back(*loop.begin());
    // 不处理走回头路的情况，即在参数域周期边界上形成了一小段的缝边loop
    assert(loop.size() > 2);

    point curr, next;
    size_t st = 0;
    curr = loop[st];
    // 向后寻找第一个不在周期边界上的点
    while(st < loop.size() - 1 && (info.is_periodic_u && equal(fabs(curr.x), info.u_ri) || info.is_periodic_v && equal(fabs(curr.y), info.v_ri))) {
        st++;
        curr = loop[st];
    }
    // 边离散保证了不会有重复点
    if(st == loop.size() - 1) assert(false);
    // size_t ed = loop.size() - 1;
    // next = loop[ed];
    //// 向前寻找的第一个不在周期边界上的点
    // while(info.is_periodic_u && equal(fabs(next.x), info.u_ri) || info.is_periodic_v && equal(fabs(next.y), info.v_ri)) {
    //     ed--;
    //     next = loop[ed];
    // }
    // ed += 1;

    // 平移使得整个loop得以顺序处理 todo
    if(st != 0) {
        loop.pop_back();
        loop.insert(loop.end(), loop.begin(), loop.begin() + st);
        loop.erase(loop.begin(), loop.begin() + st);
        loop.push_back(*loop.begin());
    }

    PeriodicPointPreprocess(loop, info);

    // 查找到第一处跨边界的地方设置为loop起始位置
    for(size_t i = 0; i < loop.size() - 1; i++) {
        curr = loop[i];
        next = loop[i + 1];
        CrossType cross_u = CrossBoundary(curr.x, next.x, info.u_tolerance);
        CrossType cross_v = CrossBoundary(curr.y, next.y, info.v_tolerance);
        if(cross_u != CrossType::kNoCross || cross_v != CrossType::kNoCross) {
            st = i + 1;
            break;
        }
    }
    loop.pop_back();
    assert(loop.size() != 0);
    loop.insert(loop.end(), loop.begin(), loop.begin() + st);
    loop.erase(loop.begin(), loop.begin() + st);
    loop.push_back(*loop.begin());
}

std::pair<int, int> GetPeriodicPointsInInterval(double inf, double sup, double pt, double period) {
    /*if(sup <= inf) {
        assert(sup > inf);
    }*/
    int inf_i = ceil((inf - pt) / period);
    int sup_i = floor((sup - pt) / period);
    if(is_equal(pt + inf_i * period, inf)) inf_i += 1;
    if(is_equal(pt + sup_i * period, sup)) sup_i -= 1;
    return std::pair<int, int>(inf_i, sup_i);
}

point_vector SingularityCompletion(double in, double out, SingularityType ty, const FaceterFaceInfo& info) {
    point_vector res;
    std::pair<int, int> interval_i;
    // ty表示奇点在哪条边界上
    switch(ty) {
        case SingularityType::kUle: {
            if(is_equal(in, out)) {
                // 如果两点位置相同 意味着这一段全都需要补或者完全不需要补，并且v方向上为周期性，并且起点位置在v的下界
                if(info.is_periodic_v && is_equal(in, info.v_le)) {
                    // 根据dv获取均匀分段的各个区间，左端点非0
                    interval_i = GetPeriodicPointsInInterval(info.v_le, info.v_ri, info.v_le, info.dv);
                    for(int i = interval_i.first; i <= interval_i.second; i++) {
                        res.push_back(point(info.u_le, info.v_le + i * info.dv));
                    }
                }
            } else if(info.is_periodic_v && is_greater_than(in, out)) {
                // 若起点不在下界，并且跨越边界，则两段拼接起来
                interval_i = GetPeriodicPointsInInterval(in, info.v_ri, info.v_le, info.dv);
                for(int i = interval_i.first; i <= interval_i.second; i++) {
                    res.push_back(point(info.u_le, info.v_le + i * info.dv));
                }
                interval_i = GetPeriodicPointsInInterval(info.v_le, out, info.v_le, info.dv);
                for(int i = interval_i.first; i <= interval_i.second; i++) {
                    res.push_back(point(info.u_le, info.v_le + i * info.dv));
                }
            } else if(is_less_than(in, out)) {
                interval_i = GetPeriodicPointsInInterval(in, out, info.v_le, info.dv);
                for(int i = interval_i.first; i <= interval_i.second; i++) {
                    res.push_back(point(info.u_le, info.v_le + i * info.dv));
                }
            }
            break;
        }
        case SingularityType::kUri:
            if(is_equal(in, out)) {
                if(info.is_periodic_v && is_equal(in, info.v_ri)) {
                    interval_i = GetPeriodicPointsInInterval(info.v_le, info.v_ri, info.v_le, info.dv);
                    for(int i = interval_i.second; i >= interval_i.first; i--) {
                        res.push_back(point(info.u_ri, info.v_le + i * info.dv));
                    }
                }
            } else if(is_greater_than(in, out)) {
                interval_i = GetPeriodicPointsInInterval(out, in, info.v_le, info.dv);
                for(int i = interval_i.second; i >= interval_i.first; i--) {
                    res.push_back(point(info.u_ri, info.v_le + i * info.dv));
                }
            } else if(info.is_periodic_v && is_less_than(in, out)) {
                interval_i = GetPeriodicPointsInInterval(info.v_le, in, info.v_le, info.dv);
                for(int i = interval_i.second; i >= interval_i.first; i--) {
                    res.push_back(point(info.u_ri, info.v_le + i * info.dv));
                }
                interval_i = GetPeriodicPointsInInterval(out, info.v_ri, info.v_le, info.dv);
                for(int i = interval_i.second; i >= interval_i.first; i--) {
                    res.push_back(point(info.u_ri, info.v_le + i * info.dv));
                }
            }
            break;
        case SingularityType::kVle:
            if(is_equal(in, out)) {
                if(info.is_periodic_u && is_equal(in, info.u_ri)) {
                    interval_i = GetPeriodicPointsInInterval(info.u_le, info.u_ri, info.u_le, info.du);
                    for(int i = interval_i.second; i >= interval_i.first; i--) {
                        res.push_back(point(info.u_le + i * info.du, info.v_le));
                    }
                }
            } else if(is_greater_than(in, out)) {
                interval_i = GetPeriodicPointsInInterval(out, in, info.u_le, info.du);
                for(int i = interval_i.second; i >= interval_i.first; i--) {
                    res.push_back(point(info.u_le + i * info.du, info.v_ri));
                }
            } else if(info.is_periodic_u && is_less_than(in, out)) {
                interval_i = GetPeriodicPointsInInterval(info.u_le, in, info.u_le, info.du);
                for(int i = interval_i.second; i >= interval_i.first; i--) {
                    res.push_back(point(info.u_le + i * info.du, info.v_ri));
                }
                interval_i = GetPeriodicPointsInInterval(out, info.u_ri, info.u_le, info.du);
                for(int i = interval_i.second; i >= interval_i.first; i--) {
                    res.push_back(point(info.u_le + i * info.du, info.v_ri));
                }
            }
            break;
        case SingularityType::kVri:
            if(is_equal(in, out)) {
                if(info.is_periodic_u && is_equal(in, info.u_le)) {
                    interval_i = GetPeriodicPointsInInterval(info.u_le, info.u_ri, info.u_le, info.du);
                    for(int i = interval_i.first; i <= interval_i.second; i++) {
                        res.push_back(point(info.u_le + i * info.du, info.v_ri));
                    }
                }
            } else if(info.is_periodic_u && is_greater_than(in, out)) {
                interval_i = GetPeriodicPointsInInterval(in, info.u_ri, info.u_le, info.du);
                for(int i = interval_i.first; i <= interval_i.second; i++) {
                    res.push_back(point(info.u_le + i * info.du, info.v_ri));
                }
                interval_i = GetPeriodicPointsInInterval(info.u_le, out, info.u_le, info.du);
                for(int i = interval_i.first; i <= interval_i.second; i++) {
                    res.push_back(point(info.u_le + i * info.du, info.v_ri));
                }
            } else if(is_less_than(in, out)) {
                interval_i = GetPeriodicPointsInInterval(in, out, info.u_le, info.du);
                for(int i = interval_i.first; i <= interval_i.second; i++) {
                    res.push_back(point(info.u_le + i * info.du, info.v_ri));
                }
            }
            break;
        case SingularityType::kNotSingularity:
        default:
            assert(false);
    }
    return res;
}

Clipper2Lib::PathD ToPathD(const point_vector& loop) {
    Clipper2Lib::PathD path;
    for(const auto& pt: loop) {
        path.push_back(Clipper2Lib::PointD(pt.x, pt.y));
    }
    return path;
}

point_vector ToPointVector(const Clipper2Lib::PathD& path) {
    point_vector loop;
    for(const auto& pt: path) {
        loop.push_back(point(pt.x, pt.y));
    }
    return loop;
}

bool CheckEdgeLoop(const LOOP* loop) {
    COEDGE* coedge = loop->start();
    COEDGE* coedge_nxt = coedge->next();
    if(coedge != coedge_nxt && coedge_nxt->next() == coedge && coedge->edge() == coedge_nxt->edge()) return true;
    return false;
}

FaceterEdgeLoopInfo GetEdgeLoopInfo(const point_vector& loop, const FaceterFaceInfo& info) {
    size_t st = 0;
    if(info.is_periodic_u) {
        while(st < loop.size() && is_equal(loop[st].x, info.u_le)) {
            st++;
        }
        if(st == loop.size()) {
            return FaceterEdgeLoopInfo(FaceterEdgeLoopType::kU, FaceterLoopCrossType::kNoCross);
        }
    }
    if(info.is_periodic_v) {
        st = 0;
        while(st < loop.size() && is_equal(loop[st].y, info.v_le)) {
            st++;
        }
        if(st == loop.size()) {
            return FaceterEdgeLoopInfo(FaceterEdgeLoopType::kV, FaceterLoopCrossType::kNoCross);
        }
    }
    point curr, next;
    bool is_cross_u = false, is_cross_v = false;
    for(size_t i = 0; i < loop.size() - 1; i++) {
        curr = loop[i];
        next = loop[i + 1];
        CrossType cross_u = CrossBoundary(curr.x, next.x, info.u_tolerance);
        CrossType cross_v = CrossBoundary(curr.y, next.y, info.v_tolerance);
        if(cross_u != CrossType::kNoCross && info.is_periodic_u) is_cross_u = true;
        if(cross_v != CrossType::kNoCross && info.is_periodic_v) is_cross_v = true;
    }
    if(is_cross_u && is_cross_v) return FaceterEdgeLoopInfo(FaceterEdgeLoopType::kNotBoundary, FaceterLoopCrossType::kCrossUV);
    if(is_cross_u) return FaceterEdgeLoopInfo(FaceterEdgeLoopType::kNotBoundary, FaceterLoopCrossType::kCrossU);
    if(is_cross_v) return FaceterEdgeLoopInfo(FaceterEdgeLoopType::kNotBoundary, FaceterLoopCrossType::kCrossV);
    return FaceterEdgeLoopInfo(FaceterEdgeLoopType::kNotBoundary, FaceterLoopCrossType::kNoCross);
}

void EdgeLoopProcess(point_vector& loop, const FaceterFaceInfo& info, point_vector& points, edge_vector& edges) {
    FaceterEdgeLoopInfo loop_info = GetEdgeLoopInfo(loop, info);
    if(loop_info.edge_loop_type == FaceterEdgeLoopType::kU) {
        // 添加左边界上的约束
        for(const auto& pt: loop) {
            points.push_back(point(info.u_le, pt.y));
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
        }
        edges.pop_back();
        // 添加右边界上的约束
        for(const auto& pt: loop) {
            points.push_back(point(info.u_ri, pt.y));
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
        }
        edges.pop_back();
    } else if(loop_info.edge_loop_type == FaceterEdgeLoopType::kV) {
        // 添加下边界上的约束
        for(const auto& pt: loop) {
            points.push_back(point(pt.x, info.v_le));
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
        }
        edges.pop_back();
        // 添加上边界上的约束
        for(const auto& pt: loop) {
            points.push_back(point(pt.x, info.v_ri));
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
        }
        edges.pop_back();
    } else {
        if(loop_info.cross_type == FaceterLoopCrossType::kNoCross) {
            for(const auto& point: loop) {
                points.push_back(point);
                edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            }
            edges.pop_back();
        } else {
            size_t st, ed;
            st = ed = 0;
            point curr, next;
            // 在此之前需要进行跨越周期的点添加处理 此处不需要保证首点为跨周期点
            PeriodicPointPreprocess(loop, info);

            // 遍历查找所有跨越边界的段
            st = 0;
            std::vector<point_vector> segments;
            for(size_t i = 0; i < loop.size() - 1; i++) {
                curr = loop[i];
                next = loop[i + 1];
                CrossType cross_v = CrossBoundary(curr.y, next.y, info.v_tolerance);
                CrossType cross_u = CrossBoundary(curr.y, next.y, info.u_tolerance);
                if(cross_v != CrossType::kNoCross || cross_u != CrossType::kNoCross) {
                    ed = i + 1;
                    segments.push_back(point_vector(loop.begin() + st, loop.begin() + ed));
                    st = ed;
                }
            }
            for(const auto& segment: segments) {
                for(const auto& pt: segment) {
                    points.push_back(pt);
                    edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                }
                edges.pop_back();
            }
        }
    }
}

void RefineTriangulation(FACE* f, double nt, double st, CDT::Triangulation<double>* pcdt) {
    bool res = true;

    for(size_t i = 0; i < pcdt->triangles.size(); i++) {
        CDT::Triangle t = pcdt->triangles[i];
        SPApar_pos pl[3];  // 三角形节点的参数值列表
        for(int ind = 0; ind < 3; ++ind) {
            pl[ind].u = (pcdt->vertices)[t.vertices[ind]].x;
            pl[ind].v = (pcdt->vertices)[t.vertices[ind]].y;
        }
        if(checkTriangle(f, nt, st, pl)) {
            continue;
        } else {
            // 边界处的点也直接添加会有水密性问题
            SPApar_pos mid_points[3];  // 三角形各个中点
            for(int ind = 0; ind < 2; ind++) {
                mid_points[ind] = (pl[ind] - pl[ind + 1]) / 2 + pl[ind + 1];
            }
            mid_points[2] = (pl[2] - pl[0]) / 2 + pl[0];
            pcdt->triangles.erase(pcdt->triangles.begin() + i);
            i--;

            // pcdt->triangles.push_back();
        }
    }
    for(CDT::Triangle t: pcdt->triangles) {
        SPApar_pos pl[3];  // 三角形节点的参数值列表
        for(int ind = 0; ind < 3; ++ind) {
            pl[ind].u = (pcdt->vertices)[t.vertices[ind]].x;
            pl[ind].v = (pcdt->vertices)[t.vertices[ind]].y;
        }
        if(checkTriangle(f, nt, st, pl)) {
            continue;
        } else {
            res = false;
            // 添加点.
            // 判断三角形的一边是不是在面的EDGE上
            for(int i = 0; i < 3; ++i) {
                int j = (i + 1) % 3;  // i - j为一条三角形的边
                // ppoints->push_back(CDT::V2d<double>((pl[i].u + pl[j].u) / 2, (pl[i].v + pl[j].v) / 2));
            }
        }
    }
}

small_loop_type judge_type(point_vector small_loop, CDT::V2d<double> le_low, CDT::V2d<double> ri_high) {
    size_t len = small_loop.size();
    CDT::V2d<double> start = small_loop[0];
    CDT::V2d<double> end = small_loop[len - 1];
    bool equalU = fabs(start.x - end.x) < SPAresabs;
    bool equalV = fabs(start.y - end.y) < SPAresabs;
    for(int i = 1; i < small_loop.size(); i++) {
        equalU &= (fabs(small_loop[i].x - small_loop[i - 1].x) < SPAresabs);

        equalV &= (fabs(small_loop[i].y - small_loop[i - 1].y) < SPAresabs);
    }

    bool crossUstart = fabs(start.x - le_low.x) < SPAresabs || fabs(start.x - ri_high.x) < SPAresabs;
    bool crossUend = fabs(end.x - le_low.x) < SPAresabs || fabs(end.x - ri_high.x) < SPAresabs;
    bool crossVstart = fabs(start.y - le_low.y) < SPAresabs || fabs(start.y - ri_high.y) < SPAresabs;
    bool crossVend = fabs(end.y - le_low.y) < SPAresabs || fabs(end.y - ri_high.y) < SPAresabs;
    if(len <= 2 || (equalU && equalV)) {
        return small_loop_type::degenerate;
    } else if(equalU && (crossUstart || crossUend)) {
        return small_loop_type::crossU;
    } else if(equalV && (crossVstart || crossVend)) {
        return small_loop_type::crossV;
    }
    if(crossUstart && crossVend) {
        return small_loop_type::crossUV;
    } else if(crossVstart && crossUend) {
        return small_loop_type::crossVU;
    }
    return small_loop_type::unknown;
}

boundary_edge_loop_type judge_edge_type(point_vector edge_loop, CDT::V2d<double> le_low, CDT::V2d<double> ri_high) {
    assert(edge_loop.size() >= 2);
    boundary_edge_loop_type kind = not_boundary_edge;
    point p1 = edge_loop[0];
    point p2 = edge_loop[1];
    if(equal(p1.x, le_low.x) && equal(p2.x, le_low.x)) {
        kind = be_ule;
    } else if(equal(p1.y, le_low.y) && equal(p2.y, le_low.y)) {
        kind = be_vle;
    } else if(equal(p1.x, ri_high.x) && equal(p2.x, ri_high.x)) {
        kind = be_uri;
    } else if(equal(p1.y, ri_high.y) && equal(p2.y, ri_high.y)) {
        kind = be_vri;
    }

    for(int i = 3; i < edge_loop.size(); i++) {
        switch(kind) {
            case be_vle:
                if(!equal(edge_loop[i].y, le_low.y)) return not_boundary_edge;

                break;
            case be_vri:
                if(!equal(edge_loop[i].y, ri_high.y)) return not_boundary_edge;

                break;
            case be_ule:
                if(!equal(edge_loop[i].x, le_low.x)) return not_boundary_edge;

                break;
            case be_uri:
                if(!equal(edge_loop[i].x, ri_high.x)) return not_boundary_edge;

                break;
            case not_boundary_edge:
                return not_boundary_edge;
                break;
            default:
                break;
        }
    }
    return kind;
}

PointBoundaryType GetPointBoundaryType(const point& a, const FaceterFaceInfo& info) {
    bool equal_ule = is_equal(a.x, info.u_le);
    bool equal_uri = is_equal(a.x, info.u_ri);
    bool equal_vle = is_equal(a.y, info.v_le);
    bool equal_vri = is_equal(a.y, info.v_ri);
    if(equal_ule && equal_vle) return PointBoundaryType::klele;
    if(equal_ule && equal_vri) return PointBoundaryType::kleri;
    if(equal_uri && equal_vle) return PointBoundaryType::krile;
    if(equal_uri && equal_vri) return PointBoundaryType::kriri;
    if(equal_ule) return PointBoundaryType::kUle;
    if(equal_uri) return PointBoundaryType::kUri;
    if(equal_vle) return PointBoundaryType::kVle;
    if(equal_vri) return PointBoundaryType::kVri;
    return PointBoundaryType::kNotBoundary;
    return PointBoundaryType::kNotBoundary;
}

void printLoop(point_vector loop) {
    for(int i = 0; i < loop.size(); i++) {
        printf("%lf %lf\n", loop[i].x, loop[i].y);
    }
}

bool operator==(CDT::V2d<double> a, CDT::V2d<double> b) {
    return equal(a.x, b.x) && equal(a.y, b.y);
}

bool operator!=(CDT::V2d<double> a, CDT::V2d<double> b) {
    return !(a == b);
}

logical get_facet_edge_points_and_params(EDGE* edge, SPAposition*& pos_array, double*& param_array, int& num_pts) {
    AF_POINT* start;
    AF_POINT* end;
    AF_POINT* curr;
    int knt = 0;
    pos_array = nullptr;
    param_array = nullptr;

    // Find the list of AF_POINTS on the given edge, if one exists.
    if(AF_POINT::find("gme", edge, 0, start, end)) {
        // Determine the number of points in the list.
        for(curr = start; curr != end; curr = curr->next("gme", 0)) {
            knt++;
        }
        knt++;

        // Allocate arrays of the proper size.
        pos_array = ACIS_NEW SPAposition[knt];
        param_array = ACIS_NEW double[knt];

        // Populate the arrays.
        int index = 0;
        for(curr = start; index < knt; curr = curr->next("gme", 0)) {
            pos_array[index] = curr->get_position();
            param_array[index] = curr->get_parameter();
            index++;
        }
    }
    num_pts = knt;
    return TRUE;
}

cross_type cross_boundary(double a, double b, double tol) {
    /**
     *
     * 如果a b的差值小于tol，则没有跨越边界
     * 如果大约tol，则进行如下判断
     * 如果a>b则a位于上界附近，跨越上界
     * 如果a<b则a位于下界附近，跨越下界
     */
    if(fabs(a - b) < tol) return no_cross;
    if(a > b) return cross_upper;
    return cross_lower;
}

CrossType CrossBoundary(double a, double b, double tol) {
    if(fabs(a - b) < tol) return CrossType::kNoCross;
    if(a > b) return CrossType::kCrossUpper;
    return CrossType::kCrossLower;
}

void preprocessor(std::vector<CDT::V2d<double>>& periphery, bool close_u, bool close_v, CDT::V2d<double> le_low, CDT::V2d<double> ri_high) {
    CDT::V2d<double> curr, next;
    double tu = (ri_high.x - le_low.x) / 2;
    double tv = (ri_high.y - le_low.y) / 2;
    double Tv = tv;
    assert(periphery.size() >= 2);
    curr = periphery[0];
    next = periphery[1];
    // 判断是否跨u方向
    if(close_u) {
        cross_type u_cross = cross_boundary(curr.x, next.x, tu);
        // 如果起点位于上边界，把其调正到下边界
        if(u_cross == cross_upper && equal(curr.x, ri_high.x)) {
            periphery[0].x = le_low.x;
        } else if(u_cross == cross_lower && equal(curr.x, le_low.x)) {
            periphery[0].x = ri_high.x;
        }
    }
    // 同上
    if(close_v) {
        cross_type v_cross = cross_boundary(curr.y, next.y, tv);
        if(v_cross == cross_upper && equal(curr.y, ri_high.y)) {
            periphery[0].y = le_low.y;
        } else if(v_cross == cross_lower && equal(curr.y, le_low.y)) {
            periphery[0].y = ri_high.y;
        }
    }
    // 依照上面的，按顺序取相邻的点进行调整对其到同一边界上
    for(int i = 0; i < periphery.size() - 1; i++) {
        if(i == 120 || i == 146 || i == 268) {
            int kk = 2;
        }
        curr = periphery[i];
        int nextIndex = (i + 1) % periphery.size();
        next = periphery[nextIndex];
        // 两点相同，则删除该点
        if(curr == next) {
            periphery.erase(periphery.begin() + i);
            i--;
            continue;
        }
        if(close_u) {
            cross_type u_cross = cross_boundary(curr.x, next.x, tu);
            // 在跨越边界的情况下，next点还正好位于下一处边界边上，则改动该点位置到靠近curr的点的边界边上
            if(u_cross == cross_upper && equal(next.x, le_low.x)) {
                periphery[nextIndex].x = ri_high.x;
            } else if(u_cross == cross_lower && equal(next.x, ri_high.x)) {
                periphery[nextIndex].x = le_low.x;
            }
        }

        if(close_v) {
            cross_type v_cross = cross_boundary(curr.y, next.y, tv);
            // 在跨越边界的情况下，next点还正好位于下一处边界边上，则改动该点位置到靠近curr的点的边界边上
            if(v_cross == cross_upper && equal(next.y, le_low.y)) {
                periphery[nextIndex].y = ri_high.y;
            } else if(v_cross == cross_lower && equal(next.y, ri_high.y)) {
                periphery[nextIndex].y = le_low.y;
            }
        }
    }
}

inline double cross(CDT::V2d<double> a, CDT::V2d<double> b) {
    return a.x * b.y - a.y * b.x;
}

// 两条线不退化的情况可以求解，退化的情况有待补充
bool intersect(CDT::V2d<double> p, CDT::V2d<double> pr, CDT::V2d<double> q, CDT::V2d<double> qs, CDT::V2d<double>& intersect_p) {
    CDT::V2d<double> r = {pr.x - p.x, pr.y - p.y};
    CDT::V2d<double> s = {qs.x - q.x, qs.y - q.y};

    /*if(fabs(r.x - 0) < SPAresabs && fabs(r.y - 0) < SPAresabs) {
        if(fabs(s.x - 0) < SPAresabs) {
            r.x += 20 * SPAresabs;
        } else if(fabs(s.y - 0) < SPAresabs) {
            r.y -= 20 * SPAresabs;
        } else {
            r.x += 20 * SPAresabs;
        }
    }*/

    double rs = cross(r, s);
    if(fabs(rs) <= EPSILON) return false;

    CDT::V2d<double> qpSub = {q.x - p.x, q.y - p.y};
    double u = cross(qpSub, r) / rs;
    double t = cross(qpSub, s) / rs;
    if(u < -EPSILON || u > 1.0 + EPSILON || t < -EPSILON || t > 1.0 + EPSILON) return false;
    intersect_p.x = p.x + t * r.x;
    intersect_p.y = p.y + t * r.y;
    return true;
}

double GetSegAngle(double dt, double at, double radius) {
    auto val = 1.0 - dt / radius;
    val = std::min(1.0, std::max(-1.0, val));
    auto angle = 2.0 * acos(val);
    angle = std::min(at * M_PI / 180.0, angle);
    return angle;
}

void PostProcessTriangle(FACE* f, double nt, double st, CDT::Triangulation<double>* pcdt) {
    for(int i = 0; i < pcdt->triangles.size(); i++) {
        const CDT::Triangle& t = pcdt->triangles[i];
        SPApar_pos pl[3];  // 三角形节点的参数值列表
        for(int j = 0; j < 3; j++) {
            pl[j].u = (pcdt->vertices)[t.vertices[j]].x;
            pl[j].v = (pcdt->vertices)[t.vertices[j]].y;
        }
        for(auto& neibor_id: t.neighbors) {
            if(neibor_id > pcdt->triangles.size()) continue;
            const CDT::Triangle& neibor = pcdt->triangles[neibor_id];
            SPApar_pos npl[3];
            for(int j = 0; j < 3; j++) {
                npl[j].u = (pcdt->vertices)[neibor.vertices[j]].x;
                npl[j].v = (pcdt->vertices)[neibor.vertices[j]].y;
            }
            // 获取邻居三角形上与t对着的的顶点
            const SPApar_pos& nopsv = npl[CDT::opposedVertexInd(neibor.neighbors, i)];

            // 获取t上与邻居三角形对着的顶点
            CDT::Index id = CDT::opposedVertexInd(t.neighbors, neibor_id);
            const SPApar_pos& opsv = pl[id];
            // 获取两个公共点
            const SPApar_pos& v1 = pl[CDT::cw(id)];
            const SPApar_pos& v2 = pl[CDT::ccw(id)];

            // 获取参数域两个对角线上的
            const SPAposition& pos_center = f->geometry()->equation().eval_position({(v1.u + v2.u) / 2, (v1.v + v2.v) / 2});
            const SPAposition& pos_center_main = f->geometry()->equation().eval_position({(nopsv.u + opsv.u) / 2, (nopsv.v + opsv.v) / 2});
            // 判断两个三角形参数域能否翻转
            if(pos_center != pos_center_main) {
                continue;
            }
            Vector3 center(pos_center.x(), pos_center.y(), pos_center.z());
            const SPAposition& pos_v1 = f->geometry()->equation().eval_position(v1);
            const SPAposition& pos_v2 = f->geometry()->equation().eval_position(v2);
            const SPAposition& pos_opsv = f->geometry()->equation().eval_position(opsv);
            const SPAposition& pos_nopsv = f->geometry()->equation().eval_position(nopsv);
            // 当前为副对角线状态
            const Vector3& sub_diag_center = {(pos_v1.x() + pos_v2.x()) / 2, (pos_v1.y() + pos_v2.y()) / 2, (pos_v1.z() + pos_v2.z()) / 2};

            const Vector3& diag_center = {(pos_opsv.x() + pos_nopsv.x()) / 2, (pos_opsv.y() + pos_nopsv.y()) / 2, (pos_opsv.z() + pos_nopsv.z()) / 2};
            // 换对角线
            if(center.distance(sub_diag_center) > center.distance(diag_center)) {
                pcdt->flipEdge(i, neibor_id);
                break;
            };
        }
    }
}
bool checkTriangle(FACE* f, double nt, double st, SPApar_pos* p) {
    std::vector<SPAposition> posl(5, SPAposition(0, 0, 0));
    // 三个顶点的坐标数组，3 * 3，i为顶点下标
    for(int i = 0; i < 3; ++i) {
        posl[i] = f->geometry()->equation().eval_position(p[i]);
        posl[3].set_x(posl[3].x() + posl[i].x() / 3.0);
        posl[3].set_y(posl[3].y() + posl[i].y() / 3.0);
        posl[3].set_z(posl[3].z() + posl[i].z() / 3.0);
    }
    // 重心的坐标 posl[3]在三角形上，posl[4]在曲面上。
    // posl[3] = (posl[0] + posl[1] + posl[2]) * (1.0 / 3);  // 在三角形上的几何坐标

    SPApar_pos centerPar(0, 0);
    for(int j = 0; j < 3; ++j) {
        centerPar.u += p[j].u / 3.0;  // 在曲面上的参数坐标
        centerPar.v += p[j].v / 3.0;
    }
    // 判断参数是否过近？可能需要优化
    double deltauv = SPAresabs;

    SPAvector v0 = posl[1] - posl[0];
    SPAvector v1 = posl[2] - posl[1];
    SPAvector v2 = posl[0] - posl[2];
    bool b0 = v0.len() < deltauv;
    bool b1 = v1.len() < deltauv;
    bool b2 = v2.len() < deltauv;
    if(b0 && b1 && b2) return true;

    posl[4] = f->geometry()->equation().eval_position(centerPar);  // 参数面的重心的几何坐标

    // 判断距离容差 / 表面容差
    if(st < (posl[4] - posl[3]).len()) {
        return false;
    }

    // 判断角度容差/法向容差
    SPAunit_vector nor = f->geometry()->equation().eval_normal(centerPar);

    if(b0 || b1 || b2) return true;

    v0 = normalise(v0);
    // v1 /= v1.len();
    v2 = normalise(v2);
    // if(f->sense()) nor = -nor;
    // norl[0] = temp1.Cross(temp2);  // 三角形法向

    // norl[1] = Vector3(nor.x(), nor.y(), nor.z());  // 曲面重心处法向
    // norl[1] = norl[1] * (1.0 / norl[1].length());
    SPAvector normal_approx = v0 * (-v2);
    if(normal_approx.len() < deltauv) return true;
    if(f->sense()) nor = -nor;
    normal_approx = normalise(normal_approx);
    double dot = normal_approx % nor;
    if(dot > 1.0)
        dot = 1.0 - SPAresabs / 2;
    else if(dot < -1.0)
        dot = -1.0 + SPAresabs / 2;
    double theta = 180.0 * acos(dot) / M_PI;
    if(theta > nt && theta < 180 - nt) {
        return false;
    }
    /*double sintheta = norl[0].Cross(norl[1]).length() / norl[0].length() / norl[1].length();
    if(nt > 0 && sin(nt * M_PI / 180.0) < sintheta) {
        return false;
    }*/
    return true;
}
bool checkTriangle_(FACE* f, double nt, double st, SPApar_pos* p) {
    std::vector<Vector3> posl(5, Vector3(0, 0, 0));
    // 三个顶点的坐标数组，3 * 3，i为顶点下标
    for(int i = 0; i < 3; ++i) {
        SPAposition pos = f->geometry()->equation().eval_position(p[i]);
        posl[i][0] = pos.x();
        posl[i][1] = pos.y();
        posl[i][2] = pos.z();
    }
    if(posl[0] == posl[1] || posl[1] == posl[2] || posl[2] == posl[0]) return true;
    // 重心的坐标 posl[3]在三角形上，posl[4]在曲面上。
    posl[3] = (posl[0] + posl[1] + posl[2]) * (1.0 / 3);  // 在三角形上的几何坐标
    SPApar_pos centerPar(0, 0);
    for(int j = 0; j < 3; ++j) {
        centerPar.u += p[j].u / 3.0;  // 在曲面上的参数坐标
        centerPar.v += p[j].v / 3.0;
    }
    // 判断参数是否过近？可能需要优化
    /*double deltauv = SPAresabs;
    if(posl[0].distance(posl[1]) < deltauv && posl[1].distance(posl[2]) < deltauv && posl[2].distance(posl[0]) < deltauv) {
        return true;
    }
    if((posl[1] - posl[0]).Cross(posl[2] - posl[0]).length() < deltauv) {
        return true;
    }*/

    SPAposition centerPos = f->geometry()->equation().eval_position(centerPar);  // 参数面的重心的几何坐标
    posl[4][0] = centerPos.x();
    posl[4][1] = centerPos.y();
    posl[4][2] = centerPos.z();
    // 判断距离容差 / 表面容差
    if(st > 0 && st < posl[3].distance(posl[4])) {
#if FACETER_DEBUG_MODE
        printf("! st = %f < %f = distance\n", st, posl[3].distance(posl[4]));
        printf("Triangle = (%f, %f, %f)-(%f, %f, %f)-(%f, %f, %f)\n", posl[0][0], posl[0][1], posl[0][2], posl[1][0], posl[1][1], posl[1][2], posl[2][0], posl[2][1], posl[2][2]);
#endif
        return false;
    }

    // 判断角度容差/法向容差
    std::vector<Vector3> norl(3, Vector3(0, 0, 0));
    SPAunit_vector nor = f->geometry()->equation().eval_normal(centerPar);
    norl[0] = (posl[1] - posl[0]).Cross(posl[2] - posl[0]);  // 三角形法向
    norl[1] = Vector3(nor.x(), nor.y(), nor.z());            // 曲面重心处法向
    double one = norl[0].length();
    double two = norl[1].length();
    double sintheta = norl[0].Cross(norl[1]).length() / norl[0].length() / norl[1].length();
    if(nt > 0 && sin(nt * M_PI / 180.0) < sintheta) {
#if FACETER_DEBUG_MODE
        printf("normal is in (%f, %f)\n", centerPar.u, centerPar.v);
        printf("Triangle = (%f, %f)-(%f, %f)-(%f, %f)\n", p[0].u, p[0].v, p[1].u, p[1].v, p[2].u, p[2].v);
#endif
        return false;
    }
    return true;
}
bool fixTriangle(FACE* f, double nt, double st, CDT::Triangulation<double>* pcdt, std::vector<CDT::V2d<double>>* ppoints) {
    bool res = true;

    /**
     * @todo: 合适的递归方式
     * 如果使用flag判断的效率比较低，会重复判断
     *
     * 先尝试使用只加点的方式，让CDT再进行重新连接。
     */
    for(CDT::Triangle t: pcdt->triangles) {
        SPApar_pos pl[3];  // 三角形节点的参数值列表
        for(int ind = 0; ind < 3; ++ind) {
            pl[ind].u = (pcdt->vertices)[t.vertices[ind]].x;
            pl[ind].v = (pcdt->vertices)[t.vertices[ind]].y;
        }
        if(checkTriangle(f, nt, st, pl)) {
            continue;
        } else {
            res = false;
            // 添加点.
            // 判断三角形的一边是不是在面的EDGE上
            SPApar_pos mid(0.0, 0.0);
            for(int i = 0; i < 3; ++i) {
                mid.u += pl[i].u / 3.0;
                mid.v += pl[i].v / 3.0;
            }
            ppoints->push_back(CDT::V2d<double>(mid.u, mid.v));
        }
    }
    return res;
}

void WritePointsToCSV(const std::string& filename, const point_vector& points) {
    std::ofstream file(filename);

    if(!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    file << "index,x,y\n";
    file << std::setiosflags(std::ios::fixed) << std::setprecision(8);
    for(size_t i = 0; i < points.size(); i++) {
        file << i << "," << points[i].x << "," << points[i].y << "\n";
    }
    file.close();
}

void WriteTrianglesToCSV(const std::string& filename, const CDT::TriangleVec& triangles) {
    std::ofstream file(filename);

    if(!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    file << "index,v0,v1,v2,n0,n1,n2\n";
    for(size_t i = 0; i < triangles.size(); i++) {
        file << i << "," << triangles[i].vertices[0] << "," << triangles[i].vertices[1] << "," << triangles[i].vertices[2] << "," << triangles[i].neighbors[0] << "," << triangles[i].neighbors[1] << "," << triangles[i].neighbors[2] << "\n";
    }
    file.close();
}

void SaveFaceToSAT(const char* filename, FACE* face) {
    FileInfo fileinfo;
    fileinfo.set_units(1.0);
    fileinfo.set_product_id("Example Application");
    api_set_file_info((FileIdent | FileUnits), fileinfo);
    api_set_int_option("sequence_save_files", 1);

    ENTITY_LIST entity_list_save;
    entity_list_save.add(face);

    FILE* file_ptr = fopen(filename, "w");
    if(file_ptr == nullptr) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
    }

    // 设置保存选项
    AcisOptions* ao = NULL;    // 使用默认选项
    logical text_mode = TRUE;  // 以文本模式保存

    outcome result = api_save_entity_list(file_ptr, text_mode, entity_list_save, ao);
    if(!result.ok()) {
        std::cerr << "Failed to save entity list to file." << std::endl;
        fclose(file_ptr);
        return;
    }
    fclose(file_ptr);
}

void RemoveDuilplicatePointsAndEdges(point_vector& points, edge_vector& edges) {
    int idp = 0;                                          // 用于给不重复的点分配新的索引
    std::unordered_map<CDT::V2d<double>, int> mp_points;  // 点到新索引的映射
    std::vector<int> points_mapping(points.size(), 0);    // 原始点索引到新点索引的映射
    CDT::EdgeUSet mp_edge;                                // 用于存储不重复的边
    point_vector final_points;                            // 不重复的点集
    edge_vector final_edges;                              // 不重复的边集

    for(int i = 0; i < points.size(); i++) {
        if(mp_points.find(points[i]) == mp_points.end()) {  // 如果点不在映射中，则添加
            mp_points[points[i]] = idp++;
            final_points.push_back(points[i]);  // 添加到不重复的点集中
        }
        points_mapping[i] = mp_points[points[i]];  // 更新点的映射
    }
    for(int i = 0; i < edges.size(); i++) {
        auto [x, y] = edges[i].verts();                            // 获取边的两个顶点
        CDT::Edge new_edge(points_mapping[x], points_mapping[y]);  // 根据点的新索引创建新边
        if(mp_edge.find(new_edge) == mp_edge.end()) {              // 如果边不在集合中，则添加
            mp_edge.insert(new_edge);
            final_edges.push_back(new_edge);  // 添加到不重复的边集中
        }
    }
    points = std::move(final_points);  // 更新点集
    edges = std::move(final_edges);    // 更新边集
    return;
}

/// DEBUG
void DeBugOutputUVGraph::init(FACE* f) {
    SPAinterval u_interval = f->geometry()->equation().param_range_u();
    SPAinterval v_interval = f->geometry()->equation().param_range_v();
    b_u = std::make_pair(u_interval.start_pt(), u_interval.start_pt() + u_interval.length());
    b_v = std::make_pair(v_interval.start_pt(), v_interval.start_pt() + v_interval.length());
    dis_points.clear();
    triangulation_edges.clear();
    constraint_edges.clear();
    boundary_indices.clear();
    sample_indices.clear();
}
void DeBugOutputUVGraph::AddTriangles(point_vector points, CDT::EdgeUSet edges_set) {
    for(auto edge: edges_set) {
        dis_points.emplace_back(points[edge.v1()].x, points[edge.v1()].y);
        dis_points.emplace_back(points[edge.v2()].x, points[edge.v2()].y);
        triangulation_edges.emplace_back(dis_points.size() - 1, dis_points.size() - 2);
    }
}
void DeBugOutputUVGraph::AddconstraintEdge(std::pair<double, double> pointa, std::pair<double, double> pointb) {
    dis_points.push_back(pointa);
    dis_points.push_back(pointb);
    constraint_edges.emplace_back(dis_points.size() - 1, dis_points.size() - 2);
}
void DeBugOutputUVGraph::AddBoundaryPoint(std::pair<double, double> point) {
    dis_points.push_back(point);
    boundary_indices.push_back(dis_points.size());
}
void DeBugOutputUVGraph::AddSamplePoint(std::pair<double, double> point) {
    dis_points.push_back(point);
    sample_indices.push_back(dis_points.size());
}
void DeBugOutputUVGraph::output() {
    std::ofstream outfile("data.txt");
    outfile << "Uinterval:\n";
    outfile << b_u.first << " " << b_u.second << "\n";
    outfile << "Vinterval:\n";
    outfile << b_v.first << " " << b_v.second << "\n";

    // 输出二维点集数组
    outfile << "Points:\n";
    for(const auto& p: dis_points) {
        outfile << p.first << " " << p.second << "\n";
    }

    // 输出边界点索引数组
    outfile << "Boundary Indices:\n";
    for(int idx: boundary_indices) {
        outfile << idx << "\n";
    }

    // 输出边界点索引数组
    outfile << "Sample Indices:\n";
    for(int idx: sample_indices) {
        outfile << idx << "\n";
    }

    // 输出三角化后形成的边集合
    outfile << "Triangulation Edges:\n";
    for(const auto& e: triangulation_edges) {
        outfile << e.first << " " << e.second << "\n";
    }

    // 输出限制边集合
    outfile << "Constraint Edges:\n";
    for(const auto& e: constraint_edges) {
        outfile << e.first << " " << e.second << "\n";
    }

    outfile.close();
}