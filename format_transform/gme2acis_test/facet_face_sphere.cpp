#include "acis/gme/faceter/gme_facet_face_sphere.hxx"

#include <iostream>
#include <vector>

#include "acis/gme/faceter/gme_facet_face_torus.hxx"
#include "acis/include/add_pcu.hxx"
#include "acis/include/api_sample_faces.hxx"
#include "acis/include/cucuint.hxx"
#include "acis/include/getbox.hxx"
#include "acis/include/intrapi.hxx"
#include "acis/include/pcurve.hxx"
#include "acis/include/sphere.hxx"
#include "acis/include/strdef.hxx"
#include "acis/include/transf.hxx"

struct scan_inter {
    double u;
    int low_inter;
    int up_inter;
};

void sphere_faceter::init() {
    // sphere is only closed on v , not on u
    // mainly [-pi/2,pi/2]*[-pi, pi)
    localSphere = (sphere*)&f->geometry()->equation_for_update();
    SPAinterval u_interval = f->geometry()->equation().param_range_u();
    SPAinterval v_interval = f->geometry()->equation().param_range_v();
    u_le = u_interval.start_pt();  // u的左边界
    u_ri = u_interval.end_pt();    // u的右边界
    v_le = v_interval.start_pt();  // v的左边界
    v_ri = v_interval.end_pt();    // v的右边界

    u_tolerance = 1e15;  // u上容差
    v_tolerance = 1e15;  // v上容差
    if(f->geometry()->equation().closed_u()) u_tolerance = localSphere->param_range_u().length() / 2;  // 判断跨边界所用值
    if(f->geometry()->equation().closed_v()) v_tolerance = localSphere->param_range_v().length() / 2;
    assert(equal(u_tolerance, 1e15));

    le_low.x = u_le;  // 左下角点
    le_low.y = v_le;
    ri_high.x = u_ri;  // 右上角点
    ri_high.y = v_ri;

    SPAunit_vector zeroNormal = f->geometry()->equation().eval_normal(SPApar_pos(0, 0));
    if(f->sense()) zeroNormal = -zeroNormal;
    if(zeroNormal != localSphere->uv_oridir) normalR = REVERSED;
}

void sphere_faceter::decideUVlen() {
    // 根据容差和SPHERE的几何特性确定ulen和vlen
    ulen = 2, vlen = 4;
    double r = localSphere->radius;  // 球面的半径
    // 粗调，找到满足容差的ulen和vlen的范围在2的多少次方之间
    while(true) {
        ulen *= 2, vlen *= 2;  // 每次循环的progress

        double theta = M_PI / 2 / ulen;  // 参数域上的值，用于表示划分三角形
        // 因为ulen和vlen划分的参数域只能表示一个矩形区域，用p0和p1分成两个三角形
        SPApar_pos* tmp_p0 = new SPApar_pos[3];
        SPApar_pos* tmp_p1 = new SPApar_pos[3];
        // p0的三个点
        tmp_p0[0] = SPApar_pos(0, theta);
        tmp_p0[1] = SPApar_pos(0, -theta);
        tmp_p0[2] = SPApar_pos(2 * theta, theta);
        // p1的三个点
        tmp_p1[0] = SPApar_pos(2 * theta, -theta);
        tmp_p1[1] = SPApar_pos(0, -theta);
        tmp_p1[2] = SPApar_pos(2 * theta, theta);
        if(checkTriangle(f, nt, st, tmp_p0) && checkTriangle(f, nt, st, tmp_p1)) {  // 检查三角形p0和p1是否满足容差要求
            // 如果满足要求，则说明满足容差要求的参数域范围在[2^(n-1),2^n]之间，可以进行下一步细分
            break;
        }
    }
    // 退回到不满足要求的最大参数范围，寻找最合适的参数范围。
    // 个人见解：可以类比成二分查找。
    ulen /= 2, vlen /= 2;
    int bias_u = ulen / 2, bias_v = vlen / 2;  // bias每次步进的大小
    while(bias_u >= 1 && bias_v >= 1) {
        ulen += bias_u, vlen += bias_v;  // 以ulen+bias_u和vlen+bias_v的参数域范围
        bias_u /= 2, bias_v /= 2;        // 每次循环步长缩小到原来的1/2

        double theta = M_PI / 2 / ulen;
        SPApar_pos* tmp_p0 = new SPApar_pos[3];
        SPApar_pos* tmp_p1 = new SPApar_pos[3];
        if(ulen % 2 == 0) {  // 应该可以理解是bias_u>1的时候
            // 生成p0和p1两个三角形
            tmp_p0[0] = SPApar_pos(0, theta);
            tmp_p0[1] = SPApar_pos(0, -theta);
            tmp_p0[2] = SPApar_pos(2 * theta, theta);
            tmp_p1[0] = SPApar_pos(2 * theta, -theta);
            tmp_p1[1] = SPApar_pos(0, -theta);
            tmp_p1[2] = SPApar_pos(2 * theta, theta);
        } else {
            tmp_p0[0] = SPApar_pos(theta / 2, theta / 2);
            tmp_p0[1] = SPApar_pos(-theta / 2, theta / 2);
            tmp_p0[2] = SPApar_pos(-theta / 2, -theta / 2);
            tmp_p1[0] = SPApar_pos(theta / 2, theta / 2);
            tmp_p1[1] = SPApar_pos(theta / 2, -theta / 2);
            tmp_p1[2] = SPApar_pos(-theta / 2, -theta / 2);
        }
        // 如果if成立，说明在[ulen,ulen+bias_u]和[vlen,vlen+bias_v]之间还有更合适的取值
        // 如果if不成立，说明在ulen+bias_u和vlen+bias_v的取值不合适，需要从[ulen+bias_u,ulen+2*bias_u]和[vlen+bias_v,ulen+2*bias_u]的范围里继续寻找
        if(checkTriangle(f, nt, st, tmp_p0) && checkTriangle(f, nt, st, tmp_p1)) {
            ulen -= 2 * bias_u, vlen -= 2 * bias_v;
            continue;
        }
    }
    ulen += 1, vlen += 2;
    du = (u_ri - u_le) / ulen;
    dv = (v_ri - v_le) / vlen;
}

// store the points of each crossing-point pairs, with the first point's index
// static std::vector<std::tuple<int, point, point>> cross_pairs;

bool sphere_faceter::getLoops() {
    // in sphere, v_speration is not valid.
    bool has_periphery_loop = false;
    SPApar_pos param;
    for(LOOP* lp = f->loop(); lp; lp = lp->next()) {
        std::vector<std::tuple<int, point, point>> cross_pairs;
        std::vector<CDT::V2d<double>> tempPS;
        loop_type lptype;
        api_loop_type(lp, lptype);  // 获取loop的类型
        unsigned ic = 0;
        bool periphery_has_singularity = false;  // 这个periphery loop是否有奇点
        // 获取组成该loop的所有的点
        // 根据先验知识需要对获取到的点集进行额外处理
        for(COEDGE* coedge = lp->start(); coedge != lp->start() || ic == 0; coedge = coedge->next(), ic++) {
            SPAposition* polyline;                                                   // 离散化后三维空间上的点
            double* params;                                                          // 用于计算的初始迭代点数组
            REVBIT r = coedge->sense();                                              // 判断edge的指向与coedge指向是否相反，关系到点添加顺序问题
            std::vector<CDT::V2d<double>> tempps;                                    // 一条边离散化后的点集
            int nP = 0;                                                              // 点的数量
            get_facet_edge_points_and_params(coedge->edge(), polyline, params, nP);  // 获得边离散化后的点和参数
                                                                                     // api_get_facet_edge_points();
                                                                                     // 可能造成内存泄漏；
            sg_add_pcurve_to_coedge(coedge);
            PCURVE* PC = coedge->geometry();
            // 该边为null edge
            if(PC == nullptr) {
                param = f->geometry()->equation().param(polyline[0]);
                tempps.push_back(point(param.u, le_low.y));
                tempps.push_back(point(param.u, ri_high.y));
                if(equal(param.u, ri_high.x)) {  // 如果是在右边界上，边的方向是从(param.u, ri_high.y)到(param.u, le_low.y)
                    std::reverse(tempps.begin(), tempps.end());
                }
            } else {  // 这条边是一条正常的边
                pcurve pc = PC->equation();
                for(int i = 0; i < nP; i++) {
                    param = pc.eval_position(polyline[i], params[i], 1);
                    SPApar_pos param1 = f->geometry()->equation().param(polyline[i]);

                    // printf("param: %lf \t pc:%lf %lf\tf: %lf %lf\n", params[i], param.u, param.v, param1.u, param1.v);
                    // if(tori->degenerate() && lptype == loop_type::loop_u_separation) param = coedge->geometry()->equation().eval_position(polyline[i], 1);  // @todo: 为该接口第二个参数寻找一个合适的参数par
                    // printf("%lf %lf %lf\t%lf %lf\t%lf\n", polyline[i].x(), polyline[i].y(), polyline[i].z(), param.u, param.v, params[i]);
                    if(f->geometry()->equation().closed_u() && equal(param1.u, u_le) && equal(param.u, u_ri)) {
                        param1.u = u_ri;
                    }
                    if(f->geometry()->equation().closed_v() && equal(param1.v, v_le) && equal(param.v, v_ri)) {
                        param1.v = v_ri;
                    }
                    CDT::V2d<double> temp = {param1.u, param1.v};  // 处理完后，待放入tempps的点
                    tempps.push_back(temp);
                }
            }
            if(r) std::reverse(tempps.begin(), tempps.end());           // 如果coedge和edge的指向相反，则需要反转插入顺序
            tempPS.insert(tempPS.end(), tempps.begin(), tempps.end());  // 向tempPS插入这个loop的点
        }

        if(normalR) {
            std::reverse(tempPS.begin(), tempPS.end());
        }

        // delete duplicate points, so that the same process won't run again in prepocessor
        for(int i = 0; i < tempPS.size() - 1; i++) {
            auto prev = tempPS[i];
            auto next = tempPS[i + 1];
            if(prev == next) tempPS.erase(tempPS.begin() + i);
        }

        // the following code deals with the singularities
        SPAposition sphere_centre = localSphere->centre;               // 球面的球心
        SPAunit_vector sphere_pole_direction = localSphere->pole_dir;  // 球面的极坐标方向
        double sphere_radius = localSphere->radius;                    // 球面半径

        // singularities
        SPAposition upper_singularity = sphere_centre + sphere_radius * sphere_pole_direction;  // 球面上方奇点
        SPApar_pos upp_pp = localSphere->param(upper_singularity);
        SPAposition lower_singularity = sphere_centre - sphere_radius * sphere_pole_direction;
        SPApar_pos lpp_pp = localSphere->param(lower_singularity);

        for(int i = 0; i < tempPS.size();) {
            point pt = tempPS[i];  // 待处理点
            // SPAposition pt_pos = f->geometry()->equation().eval_position(SPApar_pos(pt.x, pt.y));
            // if(pt_pos == upper_singularity || pt_pos == lower_singularity) {
            if(fabs(pt.x - upp_pp.u) <= SPAresabs || fabs(pt.x - lpp_pp.u) <= SPAresabs) {                                     // 如果点pt的参数u接近某个奇点的参数u，则说明点pt映射到三维空间是奇点
                SPAposition matched_singularity = fabs(pt.x - upp_pp.u) <= SPAresabs ? upper_singularity : lower_singularity;  // 匹配的奇点在三维空间的坐标

                if(i == 0) {  // the first point is singularity
                    int j = i + 1;
                    // 找到tempPS中的奇点索引范围[i,j]
                    while(j < tempPS.size()) {
                        SPAposition following_pt = f->geometry()->equation().eval_position(SPApar_pos(tempPS[j].x, tempPS[j].y));
                        // 如果匹配，那就说明tempPS[j]也是奇点，继续寻找下一个。如果不匹配，就结束
                        if(following_pt == matched_singularity) {
                            j++;
                        } else {
                            break;
                        }
                    }
                    // 如果tempPS[j]和tempPS[j]这两个点的参数v的插值大于容差，说明跨越了边界
                    if(fabs(tempPS[j].y - tempPS[i].y) > v_tolerance) {
                        cross_pairs.push_back({i, tempPS[i], tempPS[j]});
                    }
                    // delete every duplicate points in bewteen that are singularities when mapped back to 3D space
                    tempPS.erase(tempPS.begin() + 1, tempPS.begin() + j);
                    j = i + 1;
                    // only adjust the param position of the first point
                    tempPS[i].y = tempPS[j].y;

                    point_vector addtional_singularities;
                    if(lptype == loop_periphery) {
                        if(fabs(tempPS[i].x - u_ri) <= SPAresabs && tempPS[i].y < tempPS[j].y) {  // 因为是periphery，右手面内，所以点的方向是在u_le是向上(v轴正方向)，在u_ri是向下(u轴负方向)
                            addtional_singularities.push_back(CDT::V2d(u_ri, v_le));
                            if(fabs(v_le - tempPS[i].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_ri, (v_le + tempPS[i].y) / 2.0));
                            }
                            addtional_singularities.push_back(CDT::V2d(u_ri, v_ri));
                            if(fabs(v_ri - tempPS[j].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_ri, (v_ri + tempPS[j].y) / 2.0));
                            }
                            tempPS.insert(tempPS.begin() + i + 1, addtional_singularities.begin(), addtional_singularities.end());
                        } else if(fabs(tempPS[i].x - u_le) <= SPAresabs && tempPS[i].y > tempPS[j].y) {
                            addtional_singularities.push_back(CDT::V2d(u_le, v_ri));
                            if(fabs(v_ri - tempPS[i].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_le, (v_ri + tempPS[i].y) / 2.0));
                            }
                            addtional_singularities.push_back(CDT::V2d(u_le, v_le));
                            if(fabs(v_le - tempPS[j].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_le, (v_le + tempPS[j].y) / 2.0));
                            }
                            tempPS.insert(tempPS.begin() + i + 1, addtional_singularities.begin(), addtional_singularities.end());
                        }
                    }

                    i = j + 1 + addtional_singularities.size();  // 下一个需要查看的索引
                } else if(i == tempPS.size() - 1) {              // the last point is singularity
                    // 找到tempPS中的奇点索引范围[j,i]
                    int j = i - 1;
                    // 向前找奇点
                    while(j >= 0) {
                        SPAposition following_pt = f->geometry()->equation().eval_position(SPApar_pos(tempPS[j].x, tempPS[j].y));
                        if(following_pt == matched_singularity) {
                            j--;
                        } else {
                            break;
                        }
                    }
                    // 如果tempPS[j]和tempPS[j]这两个点的参数v的插值大于容差，说明跨越了边界
                    if(fabs(tempPS[j].y - tempPS[i].y) > v_tolerance) {
                        cross_pairs.push_back({j, tempPS[j], tempPS[i]});
                    }
                    // delete every duplicate points in bewteen that are singularities when mapped back to 3D space
                    tempPS.erase(tempPS.begin() + j + 1, tempPS.end() - 1);
                    i = j + 1;
                    // only adjust the param position of the last point
                    tempPS[i].y = tempPS[j].y;

                    point_vector addtional_singularities;
                    if(lptype == loop_periphery) {  // 因为是periphery，右手面内，所以点的方向是在u_le是向上(v轴正方向)，在u_ri是向下(u轴负方向)
                        if(fabs(tempPS[i].x - u_ri) <= SPAresabs && tempPS[j].y < tempPS[i].y) {
                            addtional_singularities.push_back(CDT::V2d(u_ri, v_le));
                            if(fabs(v_le - tempPS[j].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_ri, (v_le + tempPS[j].y) / 2.0));
                            }
                            addtional_singularities.push_back(CDT::V2d(u_ri, v_ri));
                            if(fabs(v_ri - tempPS[i].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_ri, (v_ri + tempPS[i].y) / 2.0));
                            }
                            tempPS.insert(tempPS.begin() + i + 1, addtional_singularities.begin(), addtional_singularities.end());
                        } else if(fabs(tempPS[i].x - u_le) <= SPAresabs && tempPS[i].y < tempPS[j].y) {
                            addtional_singularities.push_back(CDT::V2d(u_le, v_ri));
                            if(fabs(v_ri - tempPS[j].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_le, (v_ri + tempPS[j].y) / 2.0));
                            }
                            addtional_singularities.push_back(CDT::V2d(u_le, v_le));
                            if(fabs(v_le - tempPS[i].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_le, (v_le + tempPS[i].y) / 2.0));
                            }
                            tempPS.insert(tempPS.begin() + i + 1, addtional_singularities.begin(), addtional_singularities.end());
                        }
                    }

                    i = i + addtional_singularities.size() + 1;
                } else {  // a point in the middle is singularity
                    // the starting point will be adjusted based on the previous point
                    // 找到tempPS中的奇点索引范围[j,i]
                    int j = i - 1;

                    if(fabs(tempPS[j].y - tempPS[i].y) > v_tolerance) {
                        cross_pairs.push_back({j, tempPS[j], tempPS[i]});
                    }
                    tempPS[i].y = tempPS[j].y;
                    // if there is any following point that is a singularity, adjust it as well
                    j = i + 1;
                    while(j < tempPS.size()) {
                        SPAposition following_point = f->geometry()->equation().eval_position(SPApar_pos(tempPS[j].x, tempPS[j].y));
                        if(following_point == matched_singularity) {
                            j++;
                        } else {
                            break;
                        }
                    }
                    if(j == i + 1 && j < tempPS.size()) {
                        tempPS.insert(tempPS.begin() + j, pt);
                        j = i + 2;
                    }

                    while(j < tempPS.size()) {
                        SPAposition following_point = f->geometry()->equation().eval_position(SPApar_pos(tempPS[j].x, tempPS[j].y));
                        if(following_point == matched_singularity) {
                            j++;
                        } else {
                            break;
                        }
                    }

                    if(fabs(tempPS[j].y - tempPS[j - 1].y) > v_tolerance) {
                        cross_pairs.push_back({j - 1, tempPS[j - 1], tempPS[j]});
                    }
                    // then adjust the last singularity
                    tempPS[j - 1].y = tempPS[j].y;
                    //// next adjust the points in bewteen
                    // int cnt = j - i - 1;
                    // double delta_y = (tempPS[j].y - tempPS[i].y) / cnt;
                    // for(int k = i + 1; k < j - 1; k++) {
                    //     tempPS[k].y = tempPS[i].y + (k - i) * delta_y;
                    // }
                    tempPS.erase(tempPS.begin() + i + 1, tempPS.begin() + j - 1);
                    j = i + 1;

                    point_vector addtional_singularities;
                    if(lptype == loop_periphery) {
                        if(fabs(tempPS[i].x - u_ri) <= SPAresabs && tempPS[i].y < tempPS[j].y) {
                            addtional_singularities.push_back(CDT::V2d(u_ri, v_le));
                            if(fabs(v_le - tempPS[i].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_ri, (v_le + tempPS[i].y) / 2.0));
                            }
                            addtional_singularities.push_back(CDT::V2d(u_ri, v_ri));
                            if(fabs(v_ri - tempPS[j].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_ri, (v_ri + tempPS[j].y) / 2.0));
                            }
                            tempPS.insert(tempPS.begin() + i + 1, addtional_singularities.begin(), addtional_singularities.end());
                        } else if(fabs(tempPS[i].x - u_le) <= SPAresabs && tempPS[i].y > tempPS[j].y) {
                            addtional_singularities.push_back(CDT::V2d(u_le, v_ri));
                            if(fabs(v_ri - tempPS[i].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_le, (v_ri + tempPS[i].y) / 2.0));
                            }
                            addtional_singularities.push_back(CDT::V2d(u_le, v_le));
                            if(fabs(v_le - tempPS[j].y) >= v_tolerance) {
                                addtional_singularities.push_back(CDT::V2d(u_le, (v_le + tempPS[j].y) / 2.0));
                            }
                            tempPS.insert(tempPS.begin() + i + 1, addtional_singularities.begin(), addtional_singularities.end());
                        }
                    }

                    i = j + addtional_singularities.size() + 1;
                }
            } else {
                // auto next_y = tempPS[(i + 1) % tempPS.size()].y, curr_y = tempPS[i].y;
                if(fabs(tempPS[(i + 1) % tempPS.size()].y - tempPS[i].y) > v_tolerance) {
                    cross_pairs.push_back({i, tempPS[i], tempPS[(i + 1) % tempPS.size()]});
                }
                i++;
            }
        }

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
        loop_cross_pairs.push_back(cross_pairs);
    }

    return has_periphery_loop;
}

void sphere_faceter::peripheryProcess() {
    bool closed_u = f->geometry()->equation().closed_u();
    bool closed_v = f->geometry()->equation().closed_v();

    point curr, next;  // curr表示当前的点，next表示下一个点

    preprocessor(periphery, closed_u, closed_v, le_low, ri_high);  // 对periphery进行预处理
    // printLoop(periphery);

    std::vector<point_vector> small_loops;                                                   // 从跨越边界的loop中分出的不跨越边界的小loop
    boundary_loop_split(periphery, le_low, ri_high, u_tolerance, v_tolerance, small_loops);  // 对loop进行分解，分解成小的不跨越边界的loop，存在small_loops

    for(int j = 0; j < small_loops.size(); j++) {
        point_vector loop = small_loops[j];  // 当前处理的小loop中的点集
        size_t start = points.size();        // loop中点的数量
        bool boundary_point_end = false;     // 边界点结尾，
        for(int p = 0; p < loop.size(); p++) {
            points.push_back(loop[p]);
            curr = loop[p];
            next = loop[(p + 1) % loop.size()];
            if(next == curr) {  // 重复点需要丢弃
                points.pop_back();
                continue;
            }
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));  // 将两点连线的edge加入edges
            // 如果两点位于边界，需要添加在面内的点
            if((fabs(curr.x - u_ri) < EPSILON && fabs(next.x - u_ri) < EPSILON) || (fabs(curr.x - u_le) < EPSILON && fabs(next.x - u_le) < EPSILON)) {  // 如果两点在参数u的边界上
                if(fabs(curr.y - next.y) < dv) continue;                                                                                                // 如果curr和next的参数v差值小于每个分块的dv，则无需处理
                // printf("%lf %lf %lf %lf\n", curr.x, curr.y, next.x, next.y);
                if(curr.y < next.y) {
                    for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low < next.y - EPSILON; low += dv) {
                        if(low < curr.y) continue;                                     // 应该是get_cloest_boundary_point函数的精度问题，low有可能会比curr.y小。为防出现错误，判断一下
                        CDT::V2d<double> temp = {curr.x, low};                         // 需要添加的点
                        points.push_back(temp);                                        // 添加点
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));  // 添加边
                    }
                } else {
                    for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low > next.y + EPSILON; low -= dv) {
                        if(low > curr.y) continue;                                     // 同上，low有可能会比curr.y大(
                        CDT::V2d<double> temp = {curr.x, low};                         // 需要添加的点
                        points.push_back(temp);                                        // 添加点
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));  // 添加边
                    }
                }
            } else if((fabs(curr.y - v_ri) < EPSILON && fabs(next.y - v_ri) < EPSILON) || (fabs(curr.y - v_le) < EPSILON && fabs(next.y - v_le) < EPSILON)) {  // 如果两点在参数v的边界上
                if(fabs(curr.x - next.x) < du) continue;
                if(curr.x < next.x) {
                    for(double low = get_cloest_boundary_point(u_le, du, curr.x); low < next.x; low += du) {
                        if(low < curr.x) continue;
                        CDT::V2d<double> temp = {low, curr.y};
                        points.push_back(temp);
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                } else {
                    for(double low = get_cloest_boundary_point(u_le, du, curr.x); low > next.x; low -= du) {
                        if(low > curr.x) continue;
                        CDT::V2d<double> temp = {low, curr.y};
                        points.push_back(temp);
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                }
            }
        }
        edges.pop_back();
        edges.push_back(CDT::Edge(points.size() - 1, start));
    }
}

bool sphere_faceter::USeperationProcess() {
    // 存在uloop仅有一个的情况
    SPApar_pos param;                    // 参数域点
    point curr, next;                    // 当前点和后一个点
    point le_high(le_low.x, ri_high.y);  // 左上角点
    point ri_low(ri_high.x, le_low.y);   // 右下角点

    std::vector<std::vector<std::tuple<unsigned, double, bool>>> u_points;  // 记录跨域点的索引，及其y值，并且索引0记录y最小值，索引1记录y最大值
    std::vector<int> addeds;

    for(unsigned j = 0; j < u_loops.size(); j++) {
        std::vector<CDT::V2d<double>> u_loop = u_loops[j];  // loop的点集
        auto cross_pairs = loop_cross_pairs[j];             // 跨越点对

        std::vector<std::tuple<unsigned, double, bool>> u_point;  // 记录当前vloop的跨域点
        int added = 0;

        // 进行前处理，确保位于参数域边界的边不会反复横跳
        preprocessor(u_loop, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v(), le_low, ri_high);
        // printLoop(v_loop);
        // printf("\n");

        // 遍历寻找跨域点
        for(int idx = 0; idx < 1; idx++) {  // 由于preprocessor保证不会出现反复横跳，强制认为每个loop只存在一堆跨域点
            auto item = cross_pairs[idx];
            int i = get<0>(item);
            curr = get<1>(item);
            // next = get<2>(item);
            next = u_loop[(i + 1) % u_loop.size()];
            auto nxt_next = u_loop[(i + 2) % u_loop.size()];

            bool flag = i == u_loop.size() - 1;

            if(fabs(curr.y - next.y) < v_tolerance) {
                if((next.y == v_ri && nxt_next.y != v_le) || (next.y == v_le && nxt_next.y != v_ri)) {
                    u_loop.insert(u_loop.begin() + i + 2, CDT::V2d(curr.x, -next.y));
                }

                if(u_point.size() == 0) {
                    u_point.push_back({i, next.x, flag});
                    u_point.push_back({i, next.x, flag});
                } else {
                    if(next.y < get<1>(u_point[0])) {
                        get<0>(u_point[0]) = i + 1;
                        get<1>(u_point[0]) = next.x;
                    } else if(next.y > get<1>(u_point[1])) {
                        get<0>(u_point[1]) = i + 1;
                        get<1>(u_point[1]) = next.x;
                    }
                }
                added++;

                continue;
            }

            cross_type cross_v = cross_boundary(curr.y, next.y, v_tolerance);  // 跨越v边的类型

            if(cross_v == cross_upper) {
                // 指向上
                CDT::V2d<double> temp = {0, 0};                                     // 需要添加的跨边界点
                bool has_intersect = intersect(next, curr, le_low, le_high, temp);  // 查找next和curr是否和上边界有交点
                if(!has_intersect) {                                                // 如果不存在交点，就需要在边界上添加一个点
                    temp.x = curr.x;
                    temp.y = v_ri;
                }

                u_loop.insert(u_loop.begin() + i, temp);  // 将点插入u_loop

                // 仅记录最大和最小的即可
                if(u_point.size() == 0) {
                    u_point.push_back({i, temp.x, flag});
                    u_point.push_back({i, temp.x, flag});
                } else {
                    if(temp.y < get<1>(u_point[0])) {
                        get<0>(u_point[0]) = i + 1;
                        get<1>(u_point[0]) = temp.x;
                    } else if(temp.y > get<1>(u_point[1])) {
                        get<0>(u_point[1]) = i + 1;
                        get<1>(u_point[1]) = temp.x;
                    }
                }

                temp.y = v_le;                                // 在下边界也添加一个相对应的点
                u_loop.insert(u_loop.begin() + i + 1, temp);  // 将这个下边界点插入u_loop
                // i += 2;

                auto prev = u_loop[i - 1].y;
                auto nxt = u_loop[i].y;
                if(fabs(prev - nxt) >= v_tolerance && curr.x == u_ri) {
                    u_loop.insert(u_loop.begin() + i, CDT::V2d<double>(curr.x, (prev + nxt) / 2.0));
                    added++;
                }

                prev = u_loop[i + added + 1].y;
                nxt = u_loop[i + added + 2].y;
                if(fabs(prev - nxt) >= v_tolerance && curr.x == u_ri) {
                    u_loop.insert(u_loop.begin() + i + 1 + added + 1, CDT::V2d<double>(curr.x, (prev + nxt) / 2.0));
                }
            } else if(cross_v == cross_lower) {
                // 指向下
                CDT::V2d<double> temp = {0, 0};                                     // 需要添加的跨边界点
                bool has_intersect = intersect(next, curr, le_low, le_high, temp);  // 查找next和curr是否和上边界有交点
                if(!has_intersect) {                                                // 如果不存在交点，就需要在边界上添加一个点
                    temp.x = curr.x;
                    temp.y = v_le;
                }
                if(!flag)
                    u_loop.insert(u_loop.begin() + i, temp);
                else
                    u_loop.insert(u_loop.begin() + i + 1, temp);

                // 仅记录最大和最小的即可
                if(u_point.size() == 0) {
                    u_point.push_back({i, temp.x, flag});
                    u_point.push_back({i, temp.x, flag});
                } else {
                    if(temp.y < get<1>(u_point[0])) {
                        get<0>(u_point[0]) = i + 1;
                        get<1>(u_point[0]) = temp.x;
                    } else if(temp.y > get<1>(u_point[1])) {
                        get<0>(u_point[1]) = i + 1;
                        get<1>(u_point[1]) = temp.x;
                    }
                }

                temp.y = v_ri;
                if(!flag)
                    u_loop.insert(u_loop.begin() + i + 1, temp);
                else
                    u_loop.insert(u_loop.begin() + i + 1 + 1, temp);
                // i += 2;

                auto prev = u_loop[i - !flag].y;
                auto nxt = u_loop[i + flag].y;
                if(fabs(prev - nxt) >= v_tolerance && curr.x == u_ri) {
                    u_loop.insert(u_loop.begin() + i, CDT::V2d<double>(curr.x, (prev + nxt) / 2.0));
                    added++;
                }

                prev = u_loop[i + added + 1].y;
                nxt = u_loop[i + added + 2].y;
                if(fabs(prev - nxt) >= v_tolerance && curr.x == u_ri) {
                    u_loop.insert(u_loop.begin() + i + 1 + added + 1, CDT::V2d<double>(curr.x, (prev + nxt) / 2.0));
                }
            }
        }
        u_loops[j] = u_loop;
        u_points.push_back(u_point);
        addeds.push_back(added);
    }

    // 仅有一个uloop的情况
    if(u_loops.size() == 1) {
        periphery = u_loops[0];
        int added = addeds[0];
        if(u_loops[0][get<0>(u_points[0][0]) + added].y > u_loops[0][get<0>(u_points[0][0]) + added + 1].y) {
            point_vector temp_vector2;
            point_vector temp_vector;
            for(double x = get_cloest_boundary_point(u_le, du, u_loops[0][get<0>(u_points[0][0]) + added].x); x < u_ri; x += du) {
                if(x < u_loops[0][get<0>(u_points[0][0])].x) continue;
                temp_vector.push_back(point(x, v_ri));
                temp_vector2.push_back(point(x, v_le));
            }
            temp_vector.push_back(point(u_ri, v_ri));
            temp_vector.push_back(point(u_ri, 0.66 * (v_ri - v_le) + v_le));
            temp_vector.push_back(point(u_ri, 0.33 * (v_ri - v_le) + v_le));
            temp_vector.push_back(point(u_ri, v_le));
            temp_vector.insert(temp_vector.end(), temp_vector2.rbegin(), temp_vector2.rend());
            periphery.insert(periphery.begin() + get<0>(u_points[0][1]) + 1 + added, temp_vector.begin(), temp_vector.end());
        } else {
            point_vector temp_vector2;
            point_vector temp_vector;
            for(double x = get_cloest_boundary_point(u_le, du, u_loops[0][get<0>(u_points[0][0]) + added].x); x < u_ri; x += du) {
                if(x < u_loops[0][get<0>(u_points[0][0])].x) continue;
                temp_vector2.push_back(point(x, v_ri));
                temp_vector.push_back(point(x, v_le));
            }
            temp_vector.push_back(point(u_le, v_le));
            temp_vector.push_back(point(u_le, 0.33 * (v_ri - v_le) + v_le));
            temp_vector.push_back(point(u_le, 0.66 * (v_ri - v_le) + v_le));
            temp_vector.push_back(point(u_le, v_ri));
            temp_vector.insert(temp_vector.end(), temp_vector2.rbegin(), temp_vector2.rend());
            periphery.insert(periphery.begin() + get<0>(u_points[0][0]) + 1 + added, temp_vector.begin(), temp_vector.end());
        }

        return true;
    }

    // 如果第一个uloop不是指向上面的，则翻转之
    int idx1 = get<0>(u_points[0][0]) + addeds[0] + get<2>(u_points[0][0]), idx2 = get<0>(u_points[0][0]) + addeds[0] + 1 + get<2>(u_points[0][0]);
    if(u_loops[0][idx1 % u_loops[0].size()].y < u_loops[0][idx2 % u_loops[0].size()].y) {
        std::reverse(u_loops.begin(), u_loops.end());
        std::reverse(u_points.begin(), u_points.end());
        std::reverse(addeds.begin(), addeds.end());
    }

    // 有两个uloop的情况，由于其在u向上不跨越边界，故其必定顺序连接两处跨域点即可
    if(get<1>(u_points[0][1]) < get<1>(u_points[1][0])) {
        periphery = u_loops[0];
        point_vector temp_vector;
        // 不需要额外添加点以防止后续处理认为其跨越边界
        // cross_type left_right_cross = cross_boundary(u_points[0][1].second, u_points[1][0].second, v_tolerance);
        // if(left_right_cross != no_cross) temp_vector.push_back(point(u_le, (u_points[0][1].second + u_points[1][0].second) / 2));
        temp_vector.insert(temp_vector.end(), u_loops[1].begin() + get<0>(u_points[1][0]) + 1 + addeds[1] + get<2>(u_points[1][0]), u_loops[1].end());
        temp_vector.insert(temp_vector.end(), u_loops[1].begin(), u_loops[1].begin() + get<0>(u_points[1][0]) + addeds[1] + 1 + get<2>(u_points[1][0]));
        // if(left_right_cross != no_cross) temp_vector.push_back(point(u_ri, (u_points[0][1].second + u_points[1][0].second) / 2));
        periphery.insert(periphery.begin() + get<0>(u_points[0][1]) + 1 + addeds[0], temp_vector.begin(), temp_vector.end());
    } else {
        periphery = u_loops[1];
        point_vector temp_vector;
        temp_vector.insert(temp_vector.end(), u_loops[0].begin() + get<0>(u_points[0][0]) + 1 + addeds[0] + 1, u_loops[0].end());
        temp_vector.insert(temp_vector.end(), u_loops[0].begin(), u_loops[0].begin() + get<0>(u_points[0][0]) + 1 + addeds[0] + 1);
        periphery.insert(periphery.begin() + get<0>(u_points[1][0]) + 1 + addeds[1] + 1, temp_vector.begin(), temp_vector.end());
    }

    return true;
}

void sphere_faceter::holeProcess() {
    for(int i = 0; i < hole_loops.size(); i++) {
        std::vector<point_vector> small_loops;                                                       // 跨边界的loops分割成不跨边界的small_loops
        boundary_loop_split(hole_loops[i], le_low, ri_high, u_tolerance, v_tolerance, small_loops);  // 分割
        for(int j = 0; j < small_loops.size(); j++) {
            point_vector loop = small_loops[j];
            size_t start = points.size();
            for(int p = 0; p < loop.size(); p++) {
                points.push_back(loop[p]);
                edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            }
            // Edge(points.size() - 1, points.size())无效，所以需要pop出去
            edges.pop_back();
            // 因为是个loop，需要首尾相连
            edges.push_back(CDT::Edge(points.size() - 1, start));
        }
    }
}

void sphere_faceter::attachMesh() {
    // for(auto edge: edges) {
    //     point curr = points[edge.v1()];
    //     point next = points[edge.v2()];
    //     printf("%lf %lf %lf %lf\n", curr.x, curr.y, next.x, next.y);
    // }
    // printf("\n");
    REVBIT r = f->sense();
    SPApar_pos param;

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

    for(double v = v_le; v < v_ri + SPAresabs / 2; v += dv) {
        // 构造一条从y = v 的直线 求交
        // 求交结果记录在inters中，以两个点为一组，两个点中间的点为inside的点，添加到points中
        SPAposition a(u_le - 2 * SPAresabs, v, 0.0);
        SPAunit_vector dir(1.0, 0.0, 0.0);
        straight v_line(a, dir);
        SPAinterval line_range(0.0, u_ri - u_le + 4 * SPAresabs);
        v_line.limit(line_range);
        std::vector<scan_inter> inter_u;

        for(unsigned i = 0; i < edges.size(); i++) {
            point curr = points[edges[i].v1()];
            point next = points[edges[i].v2()];
            if(curr == next) continue;
            SPAposition b(curr.x, curr.y, 0.0);
            SPAposition c(next.x, next.y, 0.0);

            SPAvector edge_dir = c - b;
            straight edge_line(b, normalise(edge_dir));
            SPAinterval edge_range(0.0, (edge_dir).len());
            edge_line.limit(edge_range);
            curve_curve_int* results_head = int_cur_cur(v_line, edge_line);
            curve_curve_int* results = results_head;

            if(results != nullptr) {
                if(results->next != nullptr) {
                    delete_curve_curve_ints(results_head);
                    continue;
                }
                double p1 = results->param1;
                SPAposition inters = v_line.eval_position(p1);
                unsigned j = 0;
                for(; j < inter_u.size(); j++) {
                    if(inter_u[j].u == inters.x()) {
                        break;
                    }
                }
                if(j == inter_u.size()) {
                    inter_u.push_back(scan_inter(inters.x(), 0, 0));
                }
                if(equal(v, curr.y) && equal(v, next.y)) {
                    inter_u[j].low_inter += 1;
                    inter_u[j].up_inter += 1;
                } else if(equal(v, curr.y) && v > next.y || v > curr.y && equal(v, next.y)) {
                    inter_u[j].low_inter += 1;
                } else if(equal(v, curr.y) && v < next.y || v < curr.y && equal(v, next.y)) {
                    inter_u[j].up_inter += 1;
                } else {
                    inter_u[j].low_inter += 1;
                    inter_u[j].up_inter += 1;
                }
            }

            delete_curve_curve_ints(results_head);
        }
        std::sort(inter_u.begin(), inter_u.end(), [](const scan_inter& a, const scan_inter& b) { return a.u < b.u; });
        for(unsigned i = 0; i < inter_u.size(); i++) {
            if((inter_u[i].up_inter + inter_u[i].low_inter) % 2 == 1) {
                if(i == 0 || i == inter_u.size() - 1) {
                    inter_u.erase(inter_u.begin() + i);
                    i -= 1;
                    continue;
                }
            }
            if(inter_u[i].low_inter % 2 == 0 && inter_u[i].up_inter % 2 == 0) {
                inter_u.erase(inter_u.begin() + i);
                i -= 1;
            }
        }

        // if (inter_u.size() % 2 != 0)
        //     assert(inter_u.size() % 2 == 0);
        for(int i = 0; i < inter_u.size(); i += 2) {
            double inf = inter_u[i].u;
            if(i + 1 >= inter_u.size()) continue;
            double sup = inter_u[i + 1].u;
            for(double u = get_cloest_boundary_point(u_le, du, inf); u < sup; u += du) {
                if(u < inf) continue;
                points.push_back(point(u, v));
            }
        }
    }

    // 调用CDT
    CDT::Triangulation<double>* pcdt = nullptr;  // CDT进行三角剖分
    int _count = 0;
    // do {
    if(pcdt) delete pcdt;
    pcdt = new CDT::Triangulation<double>;
    CDT::RemoveDuplicatesAndRemapEdges(points, edges);

    pcdt->insertVertices(points);
    pcdt->insertEdges(edges);
    pcdt->eraseOuterTrianglesAndHoles();
    //} while(_count++ < 6 && !fixTriangle(f, nt, st, pcdt, &points));

    CDT::TriangleVec triangles = pcdt->triangles;
    std::vector<CDT::V2d<double>> vertice = pcdt->vertices;
    INDEXED_MESH* mesh = new INDEXED_MESH( vertice.size(), pcdt->triangles.size(), pcdt->triangles.size() * 3);
    for(int i = 0; i < vertice.size(); i++) {
        CDT::V2d<double> pp = vertice[i];
        SPApar_pos param(pp.x, pp.y);
        SPAposition pos = f->geometry()->equation().eval_position(param);
        SPAunit_vector normal = f->geometry()->equation().eval_normal(param);
        if(r) normal = -normal;
        mesh->add_vertex( pos, normal, param);
    }

    for(int i = 0; i < triangles.size(); i++) {
        int p1 = triangles[i].vertices[0];
        int p2 = triangles[i].vertices[1];
        int p3 = triangles[i].vertices[2];
        mesh->add_polygon( i, 3);
        indexed_polygon* poly0 = mesh->get_polygon( i);
        poly0->set_vertex( 0, &mesh->get_vertex( p1));
        poly0->set_vertex( 1, &mesh->get_vertex( p2));
        poly0->set_vertex( 2, &mesh->get_vertex( p3));
    }
    attach_indexed_mesh_to_face(f, mesh);
}

outcome sphere_faceter::facet() {
    init();
    decideUVlen();

    bool has_periphery_loop = getLoops();

    if(u_loops.size() > 0) {
        has_periphery_loop = USeperationProcess();
    }

    if(has_periphery_loop) {  // 如果存在periphery类型的loop，就需要对其进行处理
        // if(periphery.size() >= 2)
        peripheryProcess();
    } else {  // 如果不存在，就对整个球面增加一个periphery类型的loop。这个loop是参数u和参数v边界构成的。
        int index1, index2;
        for(int i = 0; i <= vlen; i++) {                     // 构建参数u边界上的edge
            CDT::V2d<double> tempL = {u_le, v_le + i * dv};  // u_le边界上的edge的点
            CDT::V2d<double> tempH = {u_ri, v_le + i * dv};  // u_ri边界上的edge的点
            points.push_back(tempL);                         // 将tempL加入点集points
            points.push_back(tempH);                         // 将tempH加入点集points
            index1 = points.size() - 1;                      // tempH的索引
            index2 = points.size() - 2;                      // tempL的索引
            edges.push_back(CDT::Edge(index2, index2 + 2));  // 连接当前的tempL和下一个tempL
            edges.push_back(CDT::Edge(index1, index1 + 2));  // 连接当前的tempH和下一个tempH
            // if(i < vlen) {
            //     printf("%lf %lf %lf %lf\n", u_le, v_le + i * dv, u_le, v_le + (i + 1) * dv);
            //     printf("%lf %lf %lf %lf\n", u_ri, v_le + i * dv, u_ri, v_le + (i + 1) * dv);
            // }
        }
        // 由于最后的CDT::Edge(index2, index2 + 2)和CDT::Edge(index1, index1 + 2)是无效的，所以需要删去
        edges.pop_back();
        edges.pop_back();
        for(int i = 0; i <= ulen; i++) {                     // 构建参数v边界上的edge
            CDT::V2d<double> tempL = {u_le + i * du, v_le};  // v_le边界上的edge的点
            CDT::V2d<double> tempH = {u_le + i * du, v_ri};  // u_ri边界上的edge的点
            points.push_back(tempL);                         // 将tempL加入点集points
            points.push_back(tempH);                         // 将tempH加入点集points
            index1 = points.size() - 1;                      // tempH的索引
            index2 = points.size() - 2;                      // tempL的索引
            edges.push_back(CDT::Edge(index2, index2 + 2));  // 连接当前的tempL和下一个tempL
            edges.push_back(CDT::Edge(index1, index1 + 2));  // 连接当前的tempH和下一个tempH
            // if(i < ulen) {
            //     printf("%lf %lf %lf %lf\n", u_le + i * du, v_le, u_le + (i + 1) * du, v_le);
            //     printf("%lf %lf %lf %lf\n", u_le + i * du, v_le, u_le + (i + 1) * du, v_le);
            // }
        }
        // 由于最后的CDT::Edge(index2, index2 + 2)和CDT::Edge(index1, index1 + 2)是无效的，所以需要删去
        edges.pop_back();
        edges.pop_back();
    }

    if(hole_loops.size() > 0) {  // 如果存在hole类型的loop也需要进行处理
        holeProcess();
    }

    attachMesh();
    return outcome();
}

outcome gme_facet_face_sphere(FACE* f, double nt, double st) {
    sphere_faceter faceter(f, nt, st);
    return faceter.facet();
}
