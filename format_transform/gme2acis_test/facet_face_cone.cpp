#include "acis/gme/faceter/gme_facet_face_cone.hxx"

#include <vector>

#include "acis/gme/faceter/gme_facet_face_torus.hxx"
#include "acis/gme/faceter/gme_facet_face_utils.hxx"
#include "acis/include/acistol.hxx"
#include "acis/include/add_pcu.hxx"
#include "acis/include/api_sample_faces.hxx"
#include "acis/include/condef.hxx"
#include "acis/include/cone.hxx"
#include "acis/include/intrapi.hxx"
#include "acis/include/pcurve.hxx"
#include "acis/include/transf.hxx"

#ifdef USE_ACIS

#else
void cone_faceter::init() {
    // sphere is only closed on v , not on u
    // mainly [-pi/2,pi/2]*[-pi, pi)
    // 获取锥面
#    ifdef FACETER_USE_ACIS
    localCone = (cone*)&f->geometry()->equation_for_update();
    SPAinterval u_interval = f->geometry()->equation().param_range_u();
    SPAinterval v_interval = f->geometry()->equation().param_range_v();
#    else
    localCone = (cone*)&f->geometry()->equation_for_update();
    // 获取锥面在u上的范围
    SPAinterval u_interval = f->geometry()->equation().param_range_u();
    // 获取锥面在v上的范围
    SPAinterval v_interval = f->geometry()->equation().param_range_v();
#    endif
    // 根据u参数的范围
    switch(u_interval.type()) {
        case interval_finite:  // interval为finite时
            // u的左边界是起始点
            u_le = u_interval.start_pt();
            // u的右边界是结束点
            u_ri = u_interval.end_pt();
            break;
        case interval_finite_below:  // interval类型是finite_below时
            // u的左边界是起始点
            u_le = u_interval.start_pt();
            // u的右边界取左边界的负数
            u_ri = -u_le;
            // 如果u_le被映射到锥面是一个奇异点，需要指出
#    ifdef FACETER_USE_ACIS
            if(localCone->singular_u(u_le)) {
#    else
            if(localCone->singular_u(u_le)) {
#    endif
                singular_ule = true;
            }
            break;
        case interval_finite_above:  // interval类型是finite_abow时
            // u的右边界时结束点
            u_ri = u_interval.end_pt();
            // u的左边界取右边界的负值
            u_le = -u_ri;
            // 如果u_ri被映射到锥面是一个奇异点，需要指出
#    ifdef FACETER_USE_ACIS
            if(localCone->singular_u(u_ri)) {
#    else
            if(localCone->singular_u(u_ri)) {
#    endif
                singular_uri = true;
            }
            break;
        case interval_infinite:
        case interval_unknown:
            // 暂时通过改变系数的方案解决
            u_ri = localCone->u_param_scale;
            u_le = -u_ri;
            // is_infinite_u = true;
    }
    // 由于锥面在参数v上存在周期性，无需对其进行范围化约束
    // v的左边界是起始点
    v_le = v_interval.start_pt();
    // v的右边界是结束点
    v_ri = v_interval.end_pt();
    // 默认u,v的容差
    u_tolerance = 1e15;
    v_tolerance = 1e15;
#    ifdef FACETER_USE_ACIS
    if(f->geometry()->equation().closed_u()) u_tolerance = localCone->param_range_u().length() / 2;  // 判断跨边界所用值
    if(f->geometry()->equation().closed_v()) v_tolerance = localCone->param_range_v().length() / 2;
#    else
    if(f->geometry()->equation().closed_u()) u_tolerance = localCone->param_range_u().length() / 2;  // 判断跨边界所用值
    if(f->geometry()->equation().closed_v()) v_tolerance = localCone->param_range_v().length() / 2;
#    endif
    assert(equal(u_tolerance, 1e15));

    le_low.x = u_le;  // 左下角点
    le_low.y = v_le;
    ri_high.x = u_ri;
    ri_high.y = v_ri;
    if(localCone->cosine_angle < 0 && f->sense() == FORWARD || localCone->cosine_angle > 0 && f->sense() == REVERSED) normalR = true;
}

bool cone_faceter::getLoops() {
    // in cone, v_speration is not valid.
    bool has_periphery_loop = false;
    SPApar_pos param;  // param是一个参数域上的点

    int max_u_loop_size = 0;
#    ifdef FACETER_USE_ACIS
    for(LOOP* lp = f->loop(); lp; lp = lp->next()) {
        std::vector<CDT::V2d<double>> tempPS;
        loop_type lptype;
        api_loop_type(lp, lptype);  // 获取loop的类型
#    else
    for(LOOP* lp = f->loop(); lp; lp = lp->next()) {
        std::vector<CDT::V2d<double>> tempPS;  // tempPS用于存放二维参数域上的点
        loop_type lptype;                      // 用于记录loop的类型
        api_loop_type(lp, lptype);         // 获取loop的类型
#    endif
        unsigned ic = 0;               // coegde的数量？
        bool has_singularity = false;  // loop是否有奇异点
        bool has_null_edge = false;    // loop是否有null类型的边。
        // 获取组成该loop的所有的点
        // 根据先验知识需要对获取到的点集进行额外处理
        for(COEDGE* coedge = lp->start(); coedge != lp->start() || ic == 0; coedge = coedge->next(), ic++) {
            SPAposition* polyline;                 // 边离散化后3维空间上的点
            double* params;                        // 边离散化后2维参数空间上的点，（好像没什么用，只是get_facet_edge_points_and_params函数需要一个参数）
            REVBIT r = coedge->sense();            // 判断edge的指向与coedge指向是否相反，关系到点添加顺序问题
            std::vector<CDT::V2d<double>> tempps;  // 用于存放二维参数域上的点，最后需要合并到tempPS
            int nP = 0;                            // 用于记录有多少个点
            get_facet_edge_points_and_params(coedge->edge(), polyline, params, nP);
#    ifdef FACETER_USE_ACIS
            // api_get_facet_edge_points();
#    else
            // api_get_facet_edge_points();
#    endif
            // 可能造成内存泄漏；暂不解耦
            // sg_add_pcurve_to_coedge(coedge);
            // PCURVE* PC = coedge->geometry();
            // 该边为null edge 有可能为奇点 有可能为非流形点
            if(coedge->edge()->geometry() == nullptr) {
#    ifdef FACETER_USE_ACIS
                param = f->geometry()->equation().param(polyline[0]);
#    else
                // 因为只有一个点polyline[0]
                param = f->geometry()->equation().param(polyline[0]);
                // 将这条边变成(u,le_low.y)和(u,ri_high.y)两个点
#    endif
                tempps.push_back(point(param.u, le_low.y));
                tempps.push_back(point(param.u, ri_high.y));
                // 表明存在null边
                has_null_edge = true;
            } else {  // 正常边进行离散化
                // pcurve pc = PC->equation();
                for(int i = 0; i < nP; i++) {
                    // param = pc.eval_position(polyline[i], params[i], 1);
#    ifdef FACETER_USE_ACIS
                    SPApar_pos param1 = f->geometry()->equation().param(polyline[i]);
#    else
                    // 获得边上的离散化点
                    SPApar_pos param1 = f->geometry()->equation().param(polyline[i]);
                    // 转化至二维参数域上的点
#    endif

                    CDT::V2d<double> temp = {param1.u, param1.v};
                    // 向点集内插入点
                    tempps.push_back(temp);
                }
            }
            // if(r) std::reverse(tempps.begin(), tempps.end());

            // 提前去重
            for(size_t i = 0; i < tempps.size() - 1; i++) {
                if(tempps[i] == tempps[i + 1]) {
                    tempps.erase(tempps.begin() + i + 1);
                    i--;
                }
            }
            // 隐式的null边，点很多，但是这些点都是重复的。
            if(tempps.size() == 1) {
                // 获得点的参数u
                double tmp_u = tempps[0].x;
                // 删除这个点
                tempps.pop_back();
                // 将这条边变成(u,le_low.y)和(u,ri_high.y)两个点
                tempps.push_back(point(tmp_u, le_low.y));
                tempps.push_back(point(tmp_u, ri_high.y));
                // 表明存在null边
                has_null_edge = true;
            }

            // 奇点处理,将奇点附近的点在参数域上拉到对应位置
            if(!localCone->_IsCylinder && !has_null_edge) {
                // 获取起始点
                auto start = tempps.begin();
                // 获取起始点后的一个点
                auto start_next = tempps.begin() + 1;
                // 获取结束点
                auto end = tempps.rbegin();
                // 获取结束点前的一个点
                auto end_prev = tempps.rbegin() + 1;
                if(singular_ule) {                  // 如果奇点在ule上
                    if(is_equal(start->x, u_le)) {  // 判断开始点是否是奇点
                        assert(!is_equal(start_next->x, u_le));
                        // 用起始点的后一个点的v值来表示这个奇点处的v值(由于满足容差，所以可以用这种方法)
                        start->y = start_next->y;
                        // 标记存在奇点
                        has_singularity = true;
                    }
                    if(is_equal(end->x, u_le)) {  // 判断结束点是否是奇点
                        assert(!is_equal(end_prev->x, u_le));
                        // 用结束点前一个点的v值来表示这个奇点处的v值
                        end->y = end_prev->y;
                        // 标记该边存在奇点
                        has_singularity = true;
                    }
                } else if(singular_uri) {           // 如果奇点在uri处
                    if(is_equal(start->x, u_ri)) {  // 判断起始点是否是奇点
                        assert(!is_equal(start_next->x, u_ri));
                        // 用起始点的前一个点的v值来表示这个奇点处的v值
                        start->y = start_next->y;
                        // 标记该边存在奇点
                        has_singularity = true;
                        // 找到下一个不为奇点的点，找不到则报错
                        // 待推广到其余情况
                        /* int i = 1;
                         for(; i < tempps.size(); i++) {
                             if(equal(tempps[i].x, u_ri)) {
                                 continue;
                             }
                             break;
                         }
                         assert(i != tempps.size());
                         start->y = start_next->y;*/
                    }
                    if(is_equal(end->x, u_ri)) {  // 判断结束点是否是奇点
                        assert(!is_equal(end_prev->x, u_ri));
                        // 用结束点后一个点的v值来表示这个奇点处的v值
                        end->y = end_prev->y;
                        // 标记该边存在奇点
                        has_singularity = true;
                    }
                }
            }
            // 如果EDGE方向和COEDGE的方向相反，将tempps内的点换个方向
            if(r) std::reverse(tempps.begin(), tempps.end());
            // 将tempps插入到最终结果tempPS内
            if(tempPS.size() != 0 && *tempPS.rbegin() == *tempps.begin()) {
                tempPS.pop_back();
            }
            tempPS.insert(tempPS.end(), tempps.begin(), tempps.end());
        }
        // 根据normalR来判断是否进行倒转
        if(normalR) {
            std::reverse(tempPS.begin(), tempPS.end());
        }

        // 经过奇点的Loop，loop本身不为奇点
        if(has_singularity) {
#    ifdef FACETER_USE_ACIS
            double radius = localCone->base.GetMajorAxisLength();
#    else
            double radius = localCone->base.GetMajorAxisLength();  // 圆锥的半径长度
#    endif
            double seg_angle = GetSegAngle(this->st / 2.0, this->nt / 2.0, radius);  // 分段最小的角度值
            // 去重
            /*for(size_t i = 0; i < tempPS.size() - 1; i++) {
                if(tempPS[i] == tempPS[i + 1]) {
                    tempPS.erase(tempPS.begin() + i + 1);
                    i--;
                }
            }*/
            // 保证性质
            // loop是闭合的
            if(*tempPS.begin() != *tempPS.rbegin()) {
                tempPS.push_back(tempPS[0]);
            }

            for(size_t i = 0; i < tempPS.size() - 1; i++) {
                auto curr = tempPS[i];      // 当前的点
                auto next = tempPS[i + 1];  // 后一个点

                if(singular_uri && is_equal(curr.x, u_ri) && is_equal(next.x, u_ri)) {  // uri是奇点，且curr和next是奇点
                    double sup, inf;                                                    // sup是上界   inf是下界
                    int j_inf, j_sup;                                                   // j_inf是下界索引 j_sup是上界索引
                    point_vector addeds;                                                // 需要添加的点集
                    // 根据curr.y和next.y的关系，确定sup和inf
                    if(is_less_than(curr.y, next.y)) {  // 如果next.y-curr.y小于SPAresmch（容差），添加[v_le,curr.y]和[next.y,v_ri]上的点到addeds
                        sup = curr.y, inf = v_le;
                        j_inf = 0, j_sup = floor((sup - v_le) / seg_angle);
                        // 从上界j_sup到下界j_inf，按照seg_angle递减添加点到addeds
                        for(int j = j_sup; j >= j_inf; j--) {
                            addeds.push_back(point(u_ri, v_le + j * seg_angle));
                        }
                        // 更新下界j_inf，从next.y到v_ri，按照seg_angle递减添加点到addeds
                        j_inf = ceil((next.y - v_le) / seg_angle);
                        j_sup = (v_ri - v_le) / seg_angle;
                        for(int j = j_sup; j >= j_inf; j--) {
                            addeds.push_back(point(u_ri, v_le + j * seg_angle));
                        }
                        // 将addeds中的点插入到tempPS的i+1位置，并更新i
                        tempPS.insert(tempPS.begin() + i + 1, addeds.begin(), addeds.end());
                        i += addeds.size();  // 更新i，跳过刚插入的点
                    } else {                 // next.y-curr.y大于SPAresmch(容差)，添加[curr.y]next,y]上的点到addeds
                        sup = curr.y, inf = next.y;
                        // 根据curr.y和next.y计算j_inf和j_sup
                        j_inf = ceil((inf - v_le) / seg_angle), j_sup = floor((sup - v_le) / seg_angle);
                        // 从上界j_sup到下界j_inf，按照seg_angle递减添加点到addeds
                        for(int j = j_sup; j >= j_inf; j--) {
                            addeds.push_back(point(u_ri, v_le + j * seg_angle));
                        }
                        // 将addeds中的点插入到tempPS的i+1位置，并更新i
                        tempPS.insert(tempPS.begin() + i + 1, addeds.begin(), addeds.end());
                        i += addeds.size();  // 更新i，跳过刚插入的点
                    }
                }
                if(singular_ule && is_equal(curr.x, u_le) && is_equal(next.x, u_le)) {  // ule是奇点，且curr和next是奇点
                    double sup, inf;                                                    // sup是上界   inf是下界
                    int j_inf, j_sup;                                                   // j_inf是下界索引 j_sup是上界索引
                    point_vector addeds;                                                // 需要添加的点集
                    // 根据curr.y和next.y的关系，确定sup和inf
                    if(is_greater_than(curr.y, next.y)) {  // 如果next.y-curr.y小于SPAresmch，添加[v_le,curr.y]和[next.y,v_ri]上的点到addeds
                        sup = v_ri, inf = curr.y;
                        j_inf = ceil((curr.y - v_le) / seg_angle), j_sup = (v_ri - v_le) / seg_angle;
                        // 从上界j_sup到下界j_inf，按照seg_angle递减添加点到addeds
                        for(int j = j_inf; j <= j_sup; j++) {
                            addeds.push_back(point(u_le, v_le + j * seg_angle));
                        }
                        // 更新下界j_inf，从next.y到v_ri，按照seg_angle递减添加点到addeds
                        j_inf = 0, j_sup = floor((next.y - v_le) / seg_angle);
                        for(int j = j_inf; j <= j_sup; j++) {
                            addeds.push_back(point(u_le, v_le + j * seg_angle));
                        }
                        // 将addeds中的点插入到tempPS的i+1位置，并更新i
                        tempPS.insert(tempPS.begin() + i + 1, addeds.begin(), addeds.end());
                        i += addeds.size();  // 更新i，跳过刚插入的点
                    } else {                 // next.y-curr.y大于SPAresmch，添加[curr.y]next,y]上的点到addeds
                        sup = next.y, inf = curr.y;
                        // 根据curr.y和next.y计算j_inf和j_sup
                        j_inf = ceil((inf - v_le) / seg_angle), j_sup = floor((sup - v_le) / seg_angle);
                        // 从上界j_sup到下界j_inf，按照seg_angle递减添加点到addeds
                        for(int j = j_inf; j <= j_sup; j++) {
                            addeds.push_back(point(u_le, v_le + j * seg_angle));
                        }
                        // 将addeds中的点插入到tempPS的i+1位置，并更新i
                        tempPS.insert(tempPS.begin() + i + 1, addeds.begin(), addeds.end());
                        i += addeds.size();  // 更新i，跳过刚插入的点
                    }
                }
            }
        }

        // nulledge的loop，loop本身为奇点
        // 将这个nulledge的loop从一个点，变成[v_le,v_ri]的离散线段
        if(has_null_edge) {
#    ifdef FACETER_USE_ACIS
            double radius = localCone->base.GetMajorAxisLength();
#    else
            double radius = localCone->base.GetMajorAxisLength();  // 圆锥面的半径
#    endif
            double seg_angle = GetSegAngle(this->st / 2.0, this->nt / 2.0, radius);  // 分段最小的角度值
            point_vector addeds;                                                     // 需要增加的点集
            int j_sup = (v_ri - v_le) / seg_angle;                                   // 索引
            if(is_equal(tempPS[0].y, v_ri)) {
                for(int i = 1; i < j_sup; i++) {
                    addeds.push_back(point(u_ri, v_ri - i * seg_angle));
                }
            } else if(is_equal(tempPS[0].y, v_le)) {
                for(int i = 1; i < j_sup; i++) {
                    addeds.push_back(point(u_le, v_le + i * seg_angle));
                }
            }
            /*else {
                assert(false);
            }*/
            tempPS.insert(tempPS.begin() + 1, addeds.begin(), addeds.end());
        }
        // 如果是含有奇点的uloop或者periphery，

        primary_loops.push_back(tempPS);
        lptypes.push_back(lptype);
        switch(lptype) {
            case loop_type::loop_hole: {
                hole_loops.push_back(tempPS);
                break;
            }
            case loop_type::loop_periphery: {
                // printf("periphery\n");
                periphery = tempPS;
                has_periphery_loop = true;
                break;
            }
            case loop_type::loop_u_separation: {
                // printf("useperation \n");
                u_loops.push_back(tempPS);
                if(tempPS.size() > max_u_loop_size) max_u_loop_size = tempPS.size();
                break;
            }
            case loop_type::loop_v_separation: {
                v_loops.push_back(tempPS);
                // printLoop(tempPS);
                // printf("\n");
                break;
            }
            case loop_type::loop_uv_separation: {
                uv_loops.push_back(tempPS);
                error("loop uv seperation");
                break;
            }
            default: {
                error("unknown loop");
            }
        }
    }

    // if(is_infinite_u) {
    double u_max = u_ri, u_min = u_le;
    for(const auto& loop: primary_loops) {
        for(const auto& pt: loop) {
            if(pt.x > u_max) u_max = pt.x;
            if(pt.x < u_min) u_min = pt.x;
        }
    }
    u_ri = u_max + 2 * SPAresabs;
    u_le = u_min - 2 * SPAresabs;
    //}

    std::vector<size_t> edge_loop_index;
    size_t j = 0;
    // 进行特殊处理
#    ifdef FACETER_USE_ACIS
    for(LOOP* lp = f->loop(); lp; lp = lp->next(), j++) {
#    else
    for(LOOP* lp = f->loop(); lp; lp = lp->next(), j++) {
#    endif
        bool is_edge_loop = CheckEdgeLoop(lp);
        // 若为退化边Loop的情况
        if(is_edge_loop) {
            COEDGE* coedge = lp->start();
            SPAposition* polyline;
            double* params;
            point_vector tempps;
            int nP = 0;
            get_facet_edge_points_and_params(coedge->edge(), polyline, params, nP);
            for(int i = 0; i < nP; i++) {
#    ifdef FACETER_USE_ACIS
                param = f->geometry()->equation().param(polyline[i]);
                tempps.push_back(point(param.u, param.v));
            }
            FaceterFaceInfo face_info(u_ri, u_le, v_ri, v_le, u_tolerance, v_tolerance, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v());
#    else
                param = f->geometry()->equation().param(polyline[i]);
                tempps.push_back(point(param.u, param.v));
            }
            FaceterFaceInfo face_info(u_ri, u_le, v_ri, v_le, u_tolerance, v_tolerance, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v());
#    endif
            EdgeLoopProcess(tempps, face_info, points, edges);
            edge_loop_index.push_back(j);
        }
    }
    for(int i = edge_loop_index.size() - 1; i >= 0; i--) {
        primary_loops.erase(primary_loops.begin() + edge_loop_index[i]);
        lptypes.erase(lptypes.begin() + edge_loop_index[i]);
    }

    // writeLoopsToFile(primary_loops, "coneLoops1.txt");

    return has_periphery_loop;
}

void cone_faceter::SplitBoundarySeg() {
    // 输入的loop首尾点重复
    // loop的起始点不能为边界点开头
    for(size_t i = 0; i < primary_loops.size(); i++) {
        size_t len = primary_loops[i].size();
        for(size_t j = 0; j < len; j++) {
            if(equal(fabs(primary_loops[i][j].y), v_ri)) {
                size_t tmp_j = j + 1;
                while(equal(fabs(primary_loops[i][tmp_j % len].y), v_ri)) {
                    tmp_j++;
                }
                if(tmp_j > j + 2) {
                    boundary_segments.push_back(point_vector(primary_loops[i].begin() + j, primary_loops[i].begin() + tmp_j));
                    primary_loops[i].erase(primary_loops[i].begin() + j + 1, primary_loops[i].begin() + tmp_j);
                    j += 1;
                }
            }
        }
    }
}

void cone_faceter::CopySplitLoopSeg(const point_vector& loop) {
    // 输入的loop首尾点重复
    // loop的起始点以边界点开头
    // 输入的loop必须是跨越一次边界的
    // cone类必须是跨越v缝边的
    size_t st, ed;
    st = ed = 0;
    point curr, next;
    // 向后找到跨越边界处
    // for(size_t i = 0; i < loop.size() - 1; i++) {
    //    curr = loop[i];
    //    next = loop[i + 1];
    //    CrossType cross_v = CrossBoundary(curr.y, next.y, v_tolerance);
    //    if(cross_v != CrossType::kNoCross) {
    //        ed = i + 1;
    //    }
    //}
    //// 向前找到跨越边界处
    // for(size_t i = loop.size() - 1; i >= 1; i--) {
    //     curr = loop[i - 1];
    //     next = loop[i];
    //     CrossType cross_v = CrossBoundary(curr.y, next.y, v_tolerance);
    //     if(cross_v != CrossType::kNoCross) {
    //         st = i;
    //     }
    // }
    // point_vector tmp_segment(loop.begin() + st, loop.end());
    // tmp_segment.insert(tmp_segment.end(), loop.begin(), loop.begin() + ed);
    // boundary_segments.push_back(tmp_segment);

    // 遍历查找所有跨越边界的段
    st = 0;
    for(size_t i = 0; i < loop.size() - 1; i++) {
        curr = loop[i];
        next = loop[i + 1];
        CrossType cross_v = CrossBoundary(curr.y, next.y, v_tolerance);
        if(cross_v != CrossType::kNoCross) {
            ed = i + 1;
            split_loop_segs.push_back(point_vector(loop.begin() + st, loop.begin() + ed));
            st = ed;
        }
    }

    std::vector<point_vector> tmp_upper, tmp_lower;
    // 将被分割的段复制三份
    for(size_t i = 0; i < split_loop_segs.size(); i++) {
        bool copy_upper = true;
        bool copy_lower = true;
        if(is_equal(split_loop_segs[i].begin()->y, v_ri) && is_equal(split_loop_segs[i].rbegin()->y, v_ri)) copy_upper = false;
        if(is_equal(split_loop_segs[i].begin()->y, v_le) && is_equal(split_loop_segs[i].rbegin()->y, v_le)) copy_lower = false;
        point_vector upper;
        point_vector lower;
        for(auto& pt: split_loop_segs[i]) {
            if(copy_upper) upper.push_back(point(pt.x, pt.y + v_ri - v_le));
            if(copy_lower) lower.push_back(point(pt.x, pt.y - v_ri + v_le));
        }
        if(copy_upper) tmp_upper.push_back(upper);
        if(copy_lower) tmp_lower.push_back(lower);
    }
    split_loop_segs.insert(split_loop_segs.begin(), tmp_lower.begin(), tmp_lower.end());
    split_loop_segs.insert(split_loop_segs.end(), tmp_upper.begin(), tmp_upper.end());
}

void cone_faceter::ConnectSegments() {
    // 将boundary_segments里的所有段按序连接
    // 简单粗暴平方算法
    std::vector<bool> loop_seg_used(split_loop_segs.size(), false);
    int not_use_count = split_loop_segs.size();

    point_vector curr_loop;

    // 运行时间过长报错
    int iter_count = 0;
    while(not_use_count > 0) {
        if(curr_loop.size() == 0) {
            for(size_t i = 0; i < split_loop_segs.size(); i++) {
                if(!loop_seg_used[i]) {
                    curr_loop.insert(curr_loop.end(), split_loop_segs[i].begin(), split_loop_segs[i].end());
                    loop_seg_used[i] = true;
                    not_use_count--;
                    break;
                }
            }
            continue;
        }
        size_t count = 0;
        for(size_t i = 0; i < split_loop_segs.size(); i++, count++) {
            if(!loop_seg_used[i]) {
                if(*curr_loop.rbegin() == *split_loop_segs[i].begin()) {
                    curr_loop.insert(curr_loop.end(), split_loop_segs[i].begin() + 1, split_loop_segs[i].end());
                    loop_seg_used[i] = true;
                    not_use_count--;

                    if(*curr_loop.begin() == *curr_loop.rbegin()) {
                        connect_loops.push_back(curr_loop);
                        curr_loop.clear();
                    }
                    break;
                } else if(*curr_loop.begin() == *split_loop_segs[i].rbegin()) {
                    curr_loop.insert(curr_loop.begin(), split_loop_segs[i].begin(), split_loop_segs[i].end() - 1);
                    loop_seg_used[i] = true;
                    not_use_count--;

                    if(*curr_loop.begin() == *curr_loop.rbegin()) {
                        connect_loops.push_back(curr_loop);
                        curr_loop.clear();
                    }
                    break;
                }
            }
        }
        // 没有匹配的段，就直接加入connect_loops
        if(count == split_loop_segs.size()) {
            connect_loops.push_back(curr_loop);
            curr_loop.clear();
        }
        iter_count++;
        if(iter_count > 10000) assert(false);
    }
    if(curr_loop.size() != 0) {
        connect_loops.push_back(curr_loop);
    }
}

void cone_faceter::Clipper() {
    Clipper2Lib::FillRule fr = Clipper2Lib::FillRule::EvenOdd;
    // 构造cone类需要的周期边界，一个u向上宽度略大于原区间的矩形
    Clipper2Lib::PathsD boundary_loops;
    Clipper2Lib::PathD boundary_loop;

    boundary_loop.push_back(Clipper2Lib::PointD(u_le - 2 * SPAresabs, v_le));
    boundary_loop.push_back(Clipper2Lib::PointD(u_le - 2 * SPAresabs, v_ri));
    boundary_loop.push_back(Clipper2Lib::PointD(u_ri + 2 * SPAresabs, v_ri));
    boundary_loop.push_back(Clipper2Lib::PointD(u_ri + 2 * SPAresabs, v_le));

    boundary_loops.push_back(boundary_loop);

    Clipper2Lib::PathsD results;
    for(auto& loop: clipper_input_loops) {
        std::reverse(loop.begin(), loop.end());
    }

    bool has_uloop = false;
    int index0 = -1, index1 = -2;

    for(size_t i = 0; i < clipper_input_loops.size(); i++) {
        if(input_loops_info[i].loop_type == FaceterLoopType::kULoop) {
            if(index0 == -1) {
                index0 = i;
                has_uloop = true;
            } else {
                index1 = i;
                break;
            }
        }
    }
    if(has_uloop) {
        Clipper2Lib::PathsD u0, u1;
        u0.push_back(clipper_input_loops[index0]);
        u1.push_back(clipper_input_loops[index1]);
        results = Clipper2Lib::Intersect(u0, u1, fr, 7);

        clipper_input_loops.erase(clipper_input_loops.begin() + index1);
        clipper_input_loops.erase(clipper_input_loops.begin() + index0);
        clipper_input_loops.insert(clipper_input_loops.begin(), results.begin(), results.end());
    }

    // 除了Uloop的情况，剩下的loop均不是相交的
    results = Clipper2Lib::Intersect(clipper_input_loops, boundary_loops, fr, 7);
    for(const auto& path: results) {
        CDT_input.push_back(ToPointVector(path));
    }
    // writeLoopsToFile(CDT_input, "coneLoops2.txt");
    assert(CDT_input.size() != 0);
}

void cone_faceter::NineSheetAlgo() {
#    ifdef FACETER_USE_ACIS
    FaceterFaceInfo face_info(u_ri, u_le, v_ri, v_le, u_tolerance, v_tolerance, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v());
#    else
    FaceterFaceInfo face_info(u_ri, u_le, v_ri, v_le, u_tolerance, v_tolerance, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v());  // 整个面的信息
#    endif
    FaceterLoopInfo loop_info;  // 存放loop的信息
    for(size_t i = 0; i < primary_loops.size(); i++) {
        FPLoopPreprocessor(primary_loops[i], face_info);                   // 对FPloop进行预处理
        loop_info = CheckLoopType(primary_loops[i], face_info);            // 提取出相关的FPloop信息
        if(loop_info.loop_cross_type == FaceterLoopCrossType::kNoCross) {  // 判断fploop是否跨界，未跨界就无需处理，直接添加进return_fploop即可
            clipper_input_loops.push_back(ToPathD(primary_loops[i]));      // 将primary_loops添加进输入队列
            input_loops_info.push_back(loop_info);                         // loop_info添加进输入队列
            continue;
        }
        CopySplitLoopSeg(primary_loops[i]);                   // 对跨越边界的loop进行处理，复制loop段
        ConnectSegments();                                    // 对复制后的loop段进行连接
        if(loop_info.loop_type == FaceterLoopType::kULoop) {  // 由于v上周期性，u上没有周期性，所以只需要考虑跨越边界u的情况
            if(loop_info.is_up) {
                for(size_t j = 0; j < connect_loops.size(); j++) {
                    if(*connect_loops[j].begin() != *connect_loops[j].rbegin()) {
                        // connect_loops[j].pop_back();
                        double tmp_vri = connect_loops[j].rbegin()->y;
                        double tmp_vle = connect_loops[j].begin()->y;
                        connect_loops[j].push_back(point(u_ri + 2 * SPAresabs, tmp_vri));
                        connect_loops[j].push_back(point(u_ri + 2 * SPAresabs, tmp_vle));
                        connect_loops[j].push_back(*connect_loops[j].begin());
                    }
                }
            } else {
                for(size_t j = 0; j < connect_loops.size(); j++) {
                    if(*connect_loops[j].begin() != *connect_loops[j].rbegin()) {
                        // connect_loops[j].pop_back();
                        double tmp_vle = connect_loops[j].rbegin()->y;
                        double tmp_vri = connect_loops[j].begin()->y;
                        connect_loops[j].push_back(point(u_le - 2 * SPAresabs, tmp_vle));
                        connect_loops[j].push_back(point(u_le - 2 * SPAresabs, tmp_vri));
                        connect_loops[j].push_back(*connect_loops[j].begin());
                    }
                }
            }
        }
        for(size_t j = 0; j < connect_loops.size(); j++) {
            clipper_input_loops.push_back(ToPathD(connect_loops[j]));
            input_loops_info.push_back(loop_info);
        }
        split_loop_segs.clear();
        connect_loops.clear();
    }
    Clipper();
}

void cone_faceter::peripheryProcess() {
#    ifdef FACETER_USE_ACIS
    bool closed_u = f->geometry()->equation().closed_u();
    bool closed_v = f->geometry()->equation().closed_v();
#    else
    bool closed_u = f->geometry()->equation().closed_u();
    bool closed_v = f->geometry()->equation().closed_v();
#    endif

    point curr, next;

    preprocessor(periphery, closed_u, closed_v, le_low, ri_high);
    // printLoop(periphery);
    std::vector<point_vector> small_loops;
    boundary_loop_split(periphery, le_low, ri_high, u_tolerance, v_tolerance, small_loops);

    for(int j = 0; j < small_loops.size(); j++) {
        point_vector loop = small_loops[j];
        size_t start = points.size();
        bool boundary_point_end = false;  // 边界点结尾，
        DiscreteLoop periphery_loop;
        periphery_loop.is_periphery = true;
        for(int p = 0; p < loop.size(); p++) {
            points.push_back(loop[p]);
            curr = loop[p];
            periphery_loop.points.push_back(loop[p]);
            next = loop[(p + 1) % loop.size()];
            if(next == curr) {  // 重复点需要丢弃
                points.pop_back();
                continue;
            }

            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            // cone不需要面内采点
            // 如果两点位于边界，需要添加在面内的点
            // if((fabs(curr.x - u_ri) < EPSILON && fabs(next.x - u_ri) < EPSILON) || (fabs(curr.x - u_le) < EPSILON && fabs(next.x - u_le) < EPSILON)) {
            //     if(fabs(curr.y - next.y) < dv) continue;
            //     printf("%lf %lf %lf %lf\n", curr.x, curr.y, next.x, next.y);
            //     if(curr.y < next.y) {
            //         for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low < next.y - EPSILON; low += dv) {
            //             if(low < curr.y) continue;
            //             CDT::V2d<double> temp = {curr.x, low};
            //             points.push_back(temp);
            //             edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            //         }
            //     } else {
            //         for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low > next.y + EPSILON; low -= dv) {
            //             if(low > curr.y) continue;
            //             CDT::V2d<double> temp = {curr.x, low};
            //             points.push_back(temp);
            //             edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            //         }
            //     }
            // } else if((fabs(curr.y - v_ri) < EPSILON && fabs(next.y - v_ri) < EPSILON) || (fabs(curr.y - v_le) < EPSILON && fabs(next.y - v_le) < EPSILON)) {
            //     if(fabs(curr.x - next.x) < du) continue;
            //     if(curr.x < next.x) {
            //         for(double low = get_cloest_boundary_point(u_le, du, curr.x); low < next.x; low += du) {
            //             if(low < curr.x) continue;
            //             CDT::V2d<double> temp = {low, curr.y};
            //             points.push_back(temp);
            //             edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            //         }
            //     } else {
            //         for(double low = get_cloest_boundary_point(u_le, du, curr.x); low > next.x; low -= du) {
            //             if(low > curr.x) continue;
            //             CDT::V2d<double> temp = {low, curr.y};
            //             points.push_back(temp);
            //             edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            //         }
            //     }
            // }
        }
        edges.pop_back();
        edges.push_back(CDT::Edge(points.size() - 1, start));
        periphery_loops.push_back(periphery_loop);
    }
}

bool cone_faceter::USeperationProcess() {
    // 不存在uloop仅有一个的情况
    SPApar_pos param;
    point curr, next;
    point le_high(le_low.x, ri_high.y);
    point ri_low(ri_high.x, le_low.y);

    std::vector<std::vector<std::pair<unsigned, double>>> u_points;  // 记录跨域点的索引，及其y值，并且索引0记录y最小值，索引1记录y最大值

    // 如果其中一个loop只有一个点，认为其为孤立点，即不存在loop；
    int loops_size = u_loops.size();
    // 孤立点的u值
    double u_isolated_value = 0;
    for(unsigned j = 0; j < u_loops.size(); j++) {
        if(u_loops[j].size() == 2) {
            loops_size -= 1;
            u_isolated_value = u_loops[j][0].x;
            continue;
        }
        std::vector<CDT::V2d<double>> u_loop = u_loops[j];

        std::vector<std::pair<unsigned, double>> u_point;  // 记录当前vloop的跨域点

        // 进行前处理，确保位于参数域边界的边不会反复横跳
#    ifdef FACETER_USE_ACIS
        preprocessor(u_loop, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v(), le_low, ri_high);
#    else
        preprocessor(u_loop, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v(), le_low, ri_high);
#    endif
        // printLoop(v_loop);
        // printf("\n");

        // 遍历寻找跨域点
        for(int i = 0; i < u_loop.size(); i++) {
            curr = u_loop[i];
            next = u_loop[(i + 1) % u_loop.size()];

            cross_type cross_v = cross_boundary(curr.y, next.y, v_tolerance);

            if(cross_v == cross_upper) {
                // 指向上
                CDT::V2d<double> temp = {0, 0};
                bool has_intersect = intersect(next, curr, le_low, le_high, temp);
                if(!has_intersect) {
                    temp.x = curr.x;
                    temp.y = v_ri;
                }

                u_loop.insert(u_loop.begin() + i + 1, temp);

                // 仅记录最大和最小的即可
                if(u_point.size() == 0) {
                    u_point.push_back(std::make_pair(i + 1, temp.x));
                    u_point.push_back(std::make_pair(i + 1, temp.x));
                } else {
                    if(temp.y < u_point[0].second) {
                        u_point[0].first = i + 1;
                        u_point[0].second = temp.x;
                    } else if(temp.y > u_point[1].second) {
                        u_point[1].first = i + 1;
                        u_point[1].second = temp.x;
                    }
                }

                temp.y = v_le;
                u_loop.insert(u_loop.begin() + i + 2, temp);
                i += 2;
            } else if(cross_v == cross_lower) {
                // 指向下
                CDT::V2d<double> temp = {0, 0};
                bool has_intersect = intersect(next, curr, le_low, le_high, temp);
                if(!has_intersect) {
                    temp.x = curr.x;
                    temp.y = v_le;
                }

                u_loop.insert(u_loop.begin() + i + 1, temp);

                // 仅记录最大和最小的即可
                if(u_point.size() == 0) {
                    u_point.push_back(std::make_pair(i + 1, temp.x));
                    u_point.push_back(std::make_pair(i + 1, temp.x));
                } else {
                    if(temp.y < u_point[0].second) {
                        u_point[0].first = i + 1;
                        u_point[0].second = temp.x;
                    } else if(temp.y > u_point[1].second) {
                        u_point[1].first = i + 1;
                        u_point[1].second = temp.x;
                    }
                }

                temp.y = v_ri;
                u_loop.insert(u_loop.begin() + i + 2, temp);
                i += 2;
            }
        }

        u_loops[j] = u_loop;
        u_points.push_back(u_point);
    }

    // 仅有一个uloop的情况
    if(loops_size == 1) {
        if(u_loops[0].size() == 2) {
            std::reverse(u_loops.begin(), u_loops.end());
        }
        periphery = u_loops[0];
        if(u_loops[0][u_points[0][0].first].y > u_loops[0][u_points[0][0].first + 1].y) {
            point_vector temp_vector2;
            point_vector temp_vector;
            // for(double x = get_cloest_boundary_point(u_le, du, u_loops[0][u_points[0][0].first].x); x < u_ri; x += du) {
            //     if(x < u_loops[0][u_points[0][0].first].x) continue;
            //     temp_vector.push_back(point(x, v_ri));
            //     temp_vector2.push_back(point(x, v_le));
            // }
            temp_vector.push_back(point(u_isolated_value, v_ri));
            for(int i = u_loops[0].size() - 1; i > 0; --i) {
                temp_vector.push_back(point(u_isolated_value, i * (v_ri - v_le) / (u_loops[0].size() - 1) + v_le));
            }
            temp_vector.push_back(point(u_isolated_value, v_le));
            temp_vector.insert(temp_vector.end(), temp_vector2.rbegin(), temp_vector2.rend());
            periphery.insert(periphery.begin() + u_points[0][1].first + 1, temp_vector.begin(), temp_vector.end());
        } else {
            point_vector temp_vector2;
            point_vector temp_vector;
            // for(double x = get_cloest_boundary_point(u_le, du, u_loops[0][u_points[0][0].first].x); x < u_ri; x += du) {
            //     if(x < u_loops[0][u_points[0][0].first].x) continue;
            //     temp_vector2.push_back(point(x, v_ri));
            //     temp_vector.push_back(point(x, v_le));
            // }
            temp_vector.push_back(point(u_isolated_value, v_le));
            for(int i = 1; i < u_loops[0].size() - 1; ++i) {
                temp_vector.push_back(point(u_isolated_value, i * (v_ri - v_le) / (u_loops[0].size() - 1) + v_le));
            }
            temp_vector.push_back(point(u_isolated_value, v_ri));
            temp_vector.insert(temp_vector.end(), temp_vector2.rbegin(), temp_vector2.rend());
            periphery.insert(periphery.begin() + u_points[0][0].first + 1, temp_vector.begin(), temp_vector.end());
        }
        return true;
    }

    // 如果第一个uloop不是指向上面的，则翻转之
    if(u_loops[0][u_points[0][0].first].y < u_loops[0][u_points[0][0].first + 1].y) {
        std::reverse(u_loops.begin(), u_loops.end());
        std::reverse(u_points.begin(), u_points.end());
    }

    // 有两个uloop的情况，由于其在u向上不跨越边界，故其必定顺序连接两处跨域点即可
    if(u_points[0][1].second < u_points[1][0].second) {
        periphery = u_loops[0];
        point_vector temp_vector;
        // 不需要额外添加点以防止后续处理认为其跨越边界
        // cross_type left_right_cross = cross_boundary(u_points[0][1].second, u_points[1][0].second, v_tolerance);
        // if(left_right_cross != no_cross) temp_vector.push_back(point(u_le, (u_points[0][1].second + u_points[1][0].second) / 2));
        temp_vector.insert(temp_vector.end(), u_loops[1].begin() + u_points[1][0].first + 1, u_loops[1].end());
        temp_vector.insert(temp_vector.end(), u_loops[1].begin(), u_loops[1].begin() + u_points[1][0].first + 1);
        // if(left_right_cross != no_cross) temp_vector.push_back(point(u_ri, (u_points[0][1].second + u_points[1][0].second) / 2));
        periphery.insert(periphery.begin() + u_points[0][1].first + 1, temp_vector.begin(), temp_vector.end());
    }
    return true;
}

void cone_faceter::holeProcess() {
    for(int i = 0; i < hole_loops.size(); i++) {
        std::vector<point_vector> small_loops;
        boundary_loop_split(hole_loops[i], le_low, ri_high, u_tolerance, v_tolerance, small_loops);
        for(int j = 0; j < small_loops.size(); j++) {
            if(small_loops[j].size() > 2) {
                DiscreteLoop temp;
                temp.is_periphery = false;
                temp.points = small_loops[j];
                splited_hole_loops.push_back(temp);
            }
            point_vector loop = small_loops[j];
            size_t start = points.size();
            for(int p = 0; p < loop.size(); p++) {
                points.push_back(loop[p]);
                edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            }
            edges.pop_back();
            edges.push_back(CDT::Edge(points.size() - 1, start));
        }
    }
}

void cone_faceter::attachMesh() {
    double expansion_cofficient = 1e-3;
    int point_cnt = 0;
    // check 若为退化face则直接return

    for(const auto& loop: CDT_input) {
        point_cnt = points.size();
        for(size_t i = 0; i < loop.size(); i++) {
            points.push_back(loop[i]);
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
        }
        edges.pop_back();
        edges.push_back(CDT::Edge(points.size() - 1, point_cnt));
    }

    for(int i = 0; i < points.size(); i++) {
        points[i].x *= expansion_cofficient;
    }
    // for(auto edge: edges) {
    //     point curr = points[edge.v1()];
    //     point next = points[edge.v2()];
    //     printf("%lf %lf %lf %lf\n", curr.x, curr.y, next.x, next.y);
    // }
    // printf("\n");
    REVBIT r = f->sense();
    // SPApar_pos param;

    // 面上采点
    // assert(spline_faceter::spline_map_mesh.find(f) != spline_faceter::spline_map_mesh.end());
    // my_mesh init_mesh = spline_faceter::spline_map_mesh[f];

    // for(int i = 0; i < init_mesh.points.size(); i++) {
    //     point_face_containment judge;
    //     param.u = init_mesh.points[i].first.x;
    //     param.v = init_mesh.points[i].first.y;
    //     SPAposition pos = surf->eval_position(param);
    //     api_point_in_face(pos, f, SPAtransf(), judge, param);  // 调用该接口即可判断该点是否位于该face上
    //     if(judge == point_inside_face) {
    //         points.push_back(init_mesh.points[i].first);
    //     }
    // }

    // point_vector my_vertices;
    // std::vector<std::array<int, 3>> my_triangles;
    //// std::vector<std::array<unsigned, 3>> my_triangles;
    ////  meshGenerate(my_vertices, my_triangles);

    //// 此处去重可以在getloops时进行
    // for (int i = 0; i < loops.size(); i++) {
    //     preprocessor(loops[i], false, false, le_low, ri_high);
    //     if (*loops[i].begin() != *loops[i].rbegin()) loops[i].push_back(loops[i][0]);
    // }

    //// 判断是否存在非流形点
    // std::vector<std::vector<BoundaryPoint>> non_manifold_points(potential_non_manifold_points.size(), std::vector<BoundaryPoint>());
    // for (int i = 0; i < loops.size(); i++) {
    //     for (int j = 0; j < loops[i].size() - 1; j++) {
    //         for (int k = 0; k < potential_non_manifold_points.size(); k++) {
    //             if (loops[i][j] == potential_non_manifold_points[k]) {
    //                 non_manifold_points[k].push_back(BoundaryPoint(i, j, (j - 1 + loops[i].size()) % loops[i].size(), j+1));
    //             }
    //         }
    //     }
    // }

    // point_vector changed_points;
    // for (int i = 0; i < potential_non_manifold_points.size(); i++) {
    //     if (non_manifold_points[i].size() >= 2) {
    //         assert(non_manifold_points[i].size() == 2);
    //         BoundaryPoint curr = non_manifold_points[i][0];
    //         BoundaryPoint next = non_manifold_points[i][1];
    //         point original, a0_t, a1_t, b0_t, b1_t;
    //         original = loops[curr.loop_i][curr.i];
    //         a0_t = loops[curr.loop_i][curr.prev];
    //         a1_t = loops[curr.loop_i][curr.next];
    //         b0_t = loops[next.loop_i][next.prev];
    //         b1_t = loops[next.loop_i][next.next];

    // SPApar_pos nm_point(original.x, original.y),
    //     a0(a0_t.x, a0_t.y),
    //     a1(a1_t.x, a1_t.y),
    //     b0(b0_t.x, b0_t.y),
    //     b1(b1_t.x, b1_t.y);

    // SPApar_dir v0 ( a0 - nm_point), v1(a1 - nm_point), u0(b0-nm_point), u1(b1-nm_point);
    // SPApar_vec res0, res1;
    // res0 = (v0 + v1) / 2;
    // res1 = (u0 + u1) / 2;
    // if (greater(res1.dv, res0.dv)) {
    //     loops[curr.loop_i][curr.i].y -= 2*SPAresabs * SPAresabs;
    //     changed_points.push_back(loops[curr.loop_i][curr.i]);
    // }
    // else {
    //     loops[next.loop_i][next.i].y -= 2 * SPAresabs * SPAresabs;
    //     changed_points.push_back(loops[next.loop_i][next.i]);
    // }
    // }
    // }

    //// 进行Seidel算法离散
    // SeidelTriangulation seidelTriangulation(loops);
    // seidelTriangulation.triangulation(my_vertices, my_triangles);

    // for (int i = 0; i < my_vertices.size(); i++) {
    //     for (int j = 0; j < changed_points.size(); j++) {
    //         if (equal(my_vertices[i], changed_points[j], 3 * SPAresabs * SPAresabs)) {
    //             my_vertices[i].y += 2 * SPAresabs * SPAresabs;
    //         }
    //     }
    // }

    // writeTrianglesToFile(trianglesInfo(my_vertices, my_triangles), "triangles.txt");

    // INDEXED_MESH* mesh = new INDEXED_MESH( my_vertices.size(), my_triangles.size(), my_triangles.size() * 3);
    // for(int i = 0; i < my_vertices.size(); i++) {
    //     CDT::V2d<double> pp = my_vertices[i];
    //     SPApar_pos param(pp.x, pp.y);
    //     SPAposition pos = f->geometry()->equation().eval_position(param);
    //     SPAunit_vector normal = f->geometry()->equation().eval_normal(param);
    //     /*if(!localCone->_IsCylinder) {
    //         SPAposition apex = localCone->get_apex();
    //         if(apex == pos) {
    //             SPApar_pos param_tmp(param.u - SPAresabs, param.v);
    //             if(equal(param.u, u_le)) param_tmp.u += 2 * SPAresabs;
    //             normal = f->geometry()->equation().eval_outdir(param_tmp);
    //             printf("%lf %lf %lf\n", normal.x(), normal.y(), normal.z());
    //         }
    //     }*/
    //     if(r) normal = -normal;
    //     mesh->add_vertex( pos, normal, param);
    // }

    // for(int i = 0; i < my_triangles.size(); i++) {
    //     int p1 = my_triangles[i][0];
    //     int p2 = my_triangles[i][1];
    //     int p3 = my_triangles[i][2];
    //     if (!localCone->_IsCylinder) {
    //         SPAposition apex = localCone->get_apex();
    //         double sign = 1.0;
    //         if (singular_uri)
    //             sign *= -1.0;
    //         if (mesh->get_vertex(p1).get_position() == apex ) {
    //             double tmp_u = mesh->get_vertex(p1).get_uv().u;
    //             double tmp_v1 = mesh->get_vertex(p2).get_uv().v;
    //             double tmp_v2 = mesh->get_vertex(p3).get_uv().v;
    //             SPApar_pos tmp_p(tmp_u + sign * SPAresabs, (tmp_v1 + tmp_v2) / 2);
    //             SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
    //             if (f->sense()) tmp_n = -tmp_n;
    //             // printf("%lf %lf\t%lf %lf %lf\n", tmp_p.u, tmp_p.v, tmp_n.x(), tmp_n.y(), tmp_n.z());
    //             // mesh->get_vertex(p1).set_uv(tmp_p);
    //             mesh->get_vertex(p1).set_normal(tmp_n);
    //         }
    //         else if (mesh->get_vertex(p2).get_position() == apex ) {
    //             double tmp_u = mesh->get_vertex(p2).get_uv().u;
    //             double tmp_v1 = mesh->get_vertex(p1).get_uv().v;
    //             double tmp_v2 = mesh->get_vertex(p3).get_uv().v;
    //             SPApar_pos tmp_p(tmp_u + sign * SPAresabs, (tmp_v1 + tmp_v2) / 2);
    //             SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
    //             if (f->sense()) tmp_n = -tmp_n;
    //            // printf("%lf %lf\t%lf %lf %lf\n", tmp_p.u, tmp_p.v, tmp_n.x(), tmp_n.y(), tmp_n.z());
    //             // mesh->get_vertex(p2).set_uv(tmp_p);
    //             mesh->get_vertex(p2).set_normal(tmp_n);
    //         }
    //         else if (mesh->get_vertex(p3).get_position() == apex) {
    //             double tmp_u = mesh->get_vertex(p3).get_uv().u;
    //             double tmp_v1 = mesh->get_vertex(p1).get_uv().v;
    //             double tmp_v2 = mesh->get_vertex(p2).get_uv().v;
    //             SPApar_pos tmp_p(tmp_u + sign * SPAresabs, (tmp_v1 + tmp_v2) / 2);
    //             SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
    //             if (f->sense()) tmp_n = -tmp_n;
    //             // printf("%lf %lf\t%lf %lf %lf\n",tmp_p.u, tmp_p.v, tmp_n.x(), tmp_n.y(), tmp_n.z());
    //             // mesh->get_vertex(p3).set_uv(tmp_p);
    //             mesh->get_vertex(p3).set_normal(tmp_n);
    //         }
    //     }

    // mesh->add_polygon( i, 3);
    // indexed_polygon* poly0 = mesh->get_polygon( i);
    // poly0->set_vertex( 0, &mesh->get_vertex( p1));
    // poly0->set_vertex( 1, &mesh->get_vertex( p2));
    // poly0->set_vertex( 2, &mesh->get_vertex( p3));
    // }
    // attach_indexed_mesh_to_face(f, mesh);
    // return;

    // 如需，则取消注释
    CDT::Triangulation<double>* pcdt = nullptr;  // CDT进行三角剖分
                                                 // int _count = 0;
                                                 // do {
    if(pcdt) delete pcdt;
    pcdt = new CDT::Triangulation<double>;
    CDT::RemoveDuplicatesAndRemapEdges(points, edges);

    for(int i = 0; i < points.size(); i++) {
        points[i].x /= expansion_cofficient;
    }
    // writeLoopsToFile(points, edges, "coneLoops.txt");
    for(int i = 0; i < points.size(); i++) {
        points[i].x *= expansion_cofficient;
    }
    pcdt->insertVertices(points);
    pcdt->insertEdges(edges);
    pcdt->eraseOuterTrianglesAndHoles();
    //} while(_count++ < 6 && !fixTriangle(f, nt, st, pcdt, &points));

    CDT::TriangleVec triangles = pcdt->triangles;
    std::vector<CDT::V2d<double>> vertice = pcdt->vertices;
    INDEXED_MESH* mesh = new INDEXED_MESH( vertice.size(), pcdt->triangles.size(), pcdt->triangles.size() * 3);

    for(int i = 0; i < vertice.size(); i++) {
        CDT::V2d<double> pp = vertice[i];
        SPApar_pos param(pp.x / expansion_cofficient, pp.y);
#    ifdef FACETER_USE_ACIS
        SPAposition pos = f->geometry()->equation().eval_position(param);
        SPAunit_vector normal = f->geometry()->equation().eval_outdir(param);
#    else
        SPAposition pos = f->geometry()->equation().eval_position(param);
        SPAunit_vector normal = f->geometry()->equation().eval_outdir(param);
#    endif
        if(r) normal = -normal;
        mesh->add_vertex( pos, normal, param);
    }

    for(int i = 0; i < triangles.size(); i++) {
        int p1 = triangles[i].vertices[0];
        int p2 = triangles[i].vertices[1];
        int p3 = triangles[i].vertices[2];
        if(!localCone->_IsCylinder) {
#    ifdef FACETER_USE_ACIS
            SPAposition apex = localCone->get_apex();
#    else
            SPAposition apex = localCone->get_apex();
#    endif
            double sign = 1.0;
            if(singular_uri) sign *= -1.0;
            if(mesh->get_vertex(p1).get_position() == apex) {
                double tmp_u = mesh->get_vertex(p1).get_uv().u;
                double tmp_v1 = mesh->get_vertex(p2).get_uv().v;
                double tmp_v2 = mesh->get_vertex(p3).get_uv().v;
                SPApar_pos tmp_p(tmp_u + sign * SPAresabs, (tmp_v1 + tmp_v2) / 2);
#    ifdef FACETER_USE_ACIS
                SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
#    else
                SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
#    endif
                if(f->sense()) tmp_n = -tmp_n;
                // printf("%lf %lf\t%lf %lf %lf\n", tmp_p.u, tmp_p.v, tmp_n.x(), tmp_n.y(), tmp_n.z());
                // mesh->get_vertex(p1).set_uv(tmp_p);
                mesh->get_vertex(p1).set_normal(tmp_n);
            } else if(mesh->get_vertex(p2).get_position() == apex) {
                double tmp_u = mesh->get_vertex(p2).get_uv().u;
                double tmp_v1 = mesh->get_vertex(p1).get_uv().v;
                double tmp_v2 = mesh->get_vertex(p3).get_uv().v;
                SPApar_pos tmp_p(tmp_u + sign * SPAresabs, (tmp_v1 + tmp_v2) / 2);
#    ifdef FACETER_USE_ACIS
                SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
#    else
                SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
#    endif
                if(f->sense()) tmp_n = -tmp_n;
                // printf("%lf %lf\t%lf %lf %lf\n", tmp_p.u, tmp_p.v, tmp_n.x(), tmp_n.y(), tmp_n.z());
                // mesh->get_vertex(p2).set_uv(tmp_p);
                mesh->get_vertex(p2).set_normal(tmp_n);
            } else if(mesh->get_vertex(p3).get_position() == apex) {
                double tmp_u = mesh->get_vertex(p3).get_uv().u;
                double tmp_v1 = mesh->get_vertex(p1).get_uv().v;
                double tmp_v2 = mesh->get_vertex(p2).get_uv().v;
                SPApar_pos tmp_p(tmp_u + sign * SPAresabs, (tmp_v1 + tmp_v2) / 2);
#    ifdef FACETER_USE_ACIS
                SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
#    else
                SPAunit_vector tmp_n = f->geometry()->equation().eval_normal(tmp_p);
#    endif
                if(f->sense()) tmp_n = -tmp_n;
                // printf("%lf %lf\t%lf %lf %lf\n",tmp_p.u, tmp_p.v, tmp_n.x(), tmp_n.y(), tmp_n.z());
                // mesh->get_vertex(p3).set_uv(tmp_p);
                mesh->get_vertex(p3).set_normal(tmp_n);
            }
        }
        mesh->add_polygon( i, 3);
        indexed_polygon* poly0 = mesh->get_polygon( i);
        poly0->set_vertex( 0, &mesh->get_vertex( p1));
        poly0->set_vertex( 1, &mesh->get_vertex( p2));
        poly0->set_vertex( 2, &mesh->get_vertex( p3));
    }
    attach_indexed_mesh_to_face(f, mesh);
}

void cone_faceter::PointFaceProcess() {
}

void cone_faceter::ReunitLoops(bool has_periphery) {
    std::vector<BoundaryPoint> boundary_index;
    int i = 0;
    for(auto& boundary_pts: boundary_points) {
        bool no_cross = true;
        if(primary_loops[i].size() <= 2) {
            i++;
            continue;
        }

        for(int k = 0; k < boundary_pts.size(); ++k) {
            boundary_index.push_back(BoundaryPoint(i, boundary_pts[k], boundary_pts[(k - 1 + boundary_pts.size()) % boundary_pts.size()], boundary_pts[(k + 1) % boundary_pts.size()]));
            no_cross = false;
        }
        if(no_cross) loops.push_back(primary_loops[i]);
        i++;
    }
    if(boundary_index.size() == 0) return;

    std::sort(boundary_index.begin(), boundary_index.end(), [this](const BoundaryPoint& a, const BoundaryPoint& b) { return primary_loops[a.loop_i][a.i].x < primary_loops[b.loop_i][b.i].x; });

    std::stack<point_vector> up_loops;
    std::stack<point_vector> low_loops;
    std::stack<BoundaryPoint> up_starts;
    std::stack<BoundaryPoint> low_starts;
    std::stack<BoundaryPoint> up_prevs;
    std::stack<BoundaryPoint> low_prevs;
    point_vector up_loop;
    point_vector low_loop;

    i = 0;
    BoundaryPoint temp = boundary_index[i];
    int sign = -1;
    BoundaryPoint up_start = temp;
    BoundaryPoint low_start = temp;
    BoundaryPoint up_prev = temp;
    BoundaryPoint low_prev = temp;
    up_loop.push_back(primary_loops[temp.loop_i][(temp.i + sign + primary_loops[temp.loop_i].size()) % primary_loops[temp.loop_i].size()]);
    low_loop.push_back(primary_loops[temp.loop_i][temp.i]);

    for(i = 1; i < boundary_index.size(); ++i) {
        sign *= -1;
        temp = boundary_index[i];
        if(sign == 1) {
            up_loop.push_back(primary_loops[temp.loop_i][(temp.i + sign) % primary_loops[temp.loop_i].size()]);
            // 构成了一个完整的loop，并且是从左往右连
            // todo 理清各个i的关系
            if(temp.next == up_start.i && temp.loop_i == up_start.loop_i) {
                // temp.i + 2 -> up_start.i
                if((temp.i + 2) % primary_loops[temp.loop_i].size() < up_start.i) {
                    up_loop.insert(up_loop.end(), primary_loops[temp.loop_i].begin() + (temp.i + 2) % primary_loops[temp.loop_i].size(), primary_loops[temp.loop_i].begin() + up_start.i);
                } else {
                    // assert(temp.i + 2 < primary_loops[temp.loop_i].size());
                    up_loop.insert(up_loop.end(), primary_loops[temp.loop_i].begin() + (temp.i + 2) % primary_loops[temp.loop_i].size(), primary_loops[temp.loop_i].end());
                    up_loop.insert(up_loop.end(), primary_loops[temp.loop_i].begin(), primary_loops[temp.loop_i].begin() + up_start.i);
                }
                loops.push_back(up_loop);
                up_loop.clear();

                if(!up_loops.empty()) {
                    up_loop = up_loops.top();
                    up_loops.pop();
                    up_start = up_starts.top();
                    up_starts.pop();
                    up_prev = up_prevs.top();
                    up_prevs.pop();
                }
            }
            // 无法构成完整的loop，必须再考虑下一个点
            else {
                up_prev = temp;
            }

            low_loop.insert(low_loop.begin(), primary_loops[temp.loop_i][temp.i]);
            // 底部的loop从右往左连
            if(temp.prev == low_start.i && temp.loop_i == low_start.loop_i) {
                // low_start.i -> temp.i
                if(low_start.i < temp.i)
                    low_loop.insert(low_loop.begin(), primary_loops[temp.loop_i].begin() + low_start.i, primary_loops[temp.loop_i].begin() + temp.i);
                else {
                    low_loop.insert(low_loop.begin(), primary_loops[temp.loop_i].begin(), primary_loops[temp.loop_i].begin() + temp.i);
                    low_loop.insert(low_loop.begin(), primary_loops[temp.loop_i].begin() + low_start.i, primary_loops[temp.loop_i].end());
                }
                loops.push_back(low_loop);
                low_loop.clear();

                if(!low_loops.empty()) {
                    low_loop = low_loops.top();
                    low_loops.pop();
                    low_start = low_starts.top();
                    low_starts.pop();
                    low_prev = low_prevs.top();
                    low_prevs.pop();
                }
            } else {
                low_prev = temp;
            }
        } else {
            if(up_loop.size() == 0) {
                up_loop.push_back(primary_loops[temp.loop_i][(temp.i + sign + primary_loops[temp.loop_i].size()) % primary_loops[temp.loop_i].size()]);
                up_start = temp;
            } else {
                if(up_prev.loop_i == temp.loop_i && temp.prev == up_prev.i) {
                    // up_prev.i + 2 -> temp.i
                    if((up_prev.i + 2) % primary_loops[temp.loop_i].size() < temp.i)
                        up_loop.insert(up_loop.end(), primary_loops[temp.loop_i].begin() + (up_prev.i + 2) % primary_loops[temp.loop_i].size(), primary_loops[temp.loop_i].begin() + temp.i);
                    else {
                        up_loop.insert(up_loop.end(), primary_loops[temp.loop_i].begin() + (up_prev.i + 2) % primary_loops[temp.loop_i].size(), primary_loops[temp.loop_i].end());
                        up_loop.insert(up_loop.end(), primary_loops[temp.loop_i].begin(), primary_loops[temp.loop_i].begin() + temp.i);
                    }
                } else {
                    up_loops.push(up_loop);
                    up_starts.push(up_start);
                    up_prevs.push(up_prev);
                    up_loop.clear();
                    up_loop.push_back(primary_loops[temp.loop_i][(temp.i + sign + primary_loops[temp.loop_i].size()) % primary_loops[temp.loop_i].size()]);
                    up_start = temp;
                }
            }
            up_prev = temp;

            if(low_loop.size() == 0) {
                low_loop.push_back(primary_loops[temp.loop_i][temp.i]);
                low_start = temp;
            } else {
                if(low_prev.loop_i == temp.loop_i && temp.next == low_prev.i) {
                    // temp.i -> low_prev.i
                    if(temp.i < low_prev.i)
                        low_loop.insert(low_loop.begin(), primary_loops[temp.loop_i].begin() + temp.i, primary_loops[temp.loop_i].begin() + low_prev.i);
                    else {
                        low_loop.insert(low_loop.begin(), primary_loops[temp.loop_i].begin(), primary_loops[temp.loop_i].begin() + low_prev.i);
                        low_loop.insert(low_loop.begin(), primary_loops[temp.loop_i].begin() + temp.i, primary_loops[temp.loop_i].end());
                    }
                } else {
                    low_loops.push(low_loop);
                    low_starts.push(low_start);
                    low_prevs.push(low_prev);
                    low_loop.clear();
                    low_loop.push_back(primary_loops[temp.loop_i][temp.i]);
                    low_start = temp;
                }
            }
            low_prev = temp;
        }
    }

    // 最外部剩余的up_loop和low_loop相连
    if(up_loop.size() > 0) {
        if(!(up_prev.next == low_prev.i && up_prev.loop_i == low_prev.loop_i)) {
            assert(up_prev.next == low_prev.i && up_prev.loop_i == low_prev.loop_i);
        }
        if(!(low_start.next == up_start.i && low_start.loop_i == up_start.loop_i)) {
            assert(low_start.next == up_start.i && low_start.loop_i == up_start.loop_i);
        }
        if((up_prev.i + 2) % primary_loops[up_prev.loop_i].size() < low_prev.i)
            up_loop.insert(up_loop.end(), primary_loops[up_prev.loop_i].begin() + (up_prev.i + 2) % primary_loops[up_prev.loop_i].size(), primary_loops[low_prev.loop_i].begin() + low_prev.i);
        else {
            up_loop.insert(up_loop.end(), primary_loops[up_prev.loop_i].begin() + (up_prev.i + 2) % primary_loops[up_prev.loop_i].size(), primary_loops[low_prev.loop_i].end());
            up_loop.insert(up_loop.end(), primary_loops[up_prev.loop_i].begin(), primary_loops[low_prev.loop_i].begin() + low_prev.i);
        }
        up_loop.insert(up_loop.end(), low_loop.begin(), low_loop.end());
        if((low_start.i + 1) % primary_loops[low_start.loop_i].size() < up_start.i)
            up_loop.insert(up_loop.end(), primary_loops[low_start.loop_i].begin() + (low_start.i + 1) % primary_loops[low_start.loop_i].size(), primary_loops[up_start.loop_i].begin() + up_start.i);
        else {
            up_loop.insert(up_loop.end(), primary_loops[low_start.loop_i].begin() + low_start.i + 1, primary_loops[up_start.loop_i].end());
            up_loop.insert(up_loop.end(), primary_loops[low_start.loop_i].begin(), primary_loops[up_start.loop_i].begin() + up_start.i);
        }
        loops.push_back(up_loop);
    }
}

point_boundary_type cone_faceter::is_on_boundary(point a) {
    if(equal(a.y, v_ri))
        return b_vri;
    else if(equal(a.y, v_le))
        return b_vle;
    // else if (equal(a.x, u_le))
    //   return b_ule;
    // else if (equal(a.x, u_ri))
    //  return b_uri;
    else
        return not_boundary;
}

outcome cone_faceter::facet() {
    init();

    bool has_periphery_loop = getLoops();

    NineSheetAlgo();

    attachMesh();
    return outcome();
}

outcome gme_facet_face_cone(FACE* f, double nt, double st) {
    // api_split_periodic_faces(f);
    cone_faceter faceter(f, nt, st);
    return faceter.facet();
}

#endif