#include "acis/gme/faceter/gme_facet_face_torus.hxx"

#include <iostream>

#include "acis/gme/faceter/gme_facet_face_utils.hxx"
#include "acis/include/add_pcu.hxx"
#include "acis/include/cucuint.hxx"
#include "acis/include/getbox.hxx"
#include "acis/include/intrapi.hxx"
#include "acis/include/pcurve.hxx"
#include "acis/include/strdef.hxx"
#include "acis/include/torus.hxx"
#include "acis/include/transf.hxx"

/**
 * @brief 对可能跨越边界的闭合环进行预处理和分割，确保每个分割后的子环在指定边界内闭合
 *
 * @param loop 跨边界的loop
 * @param le_low 左下角边界点
 * @param ri_high 右上角边界点
 * @param u_tolerance 容差用于u跨边界的判断
 * @param v_tolerance 容差用于v跨边界的判断
 * @param small_loops 分割得到的不跨边界的loop
 */
void boundary_loop_split(point_vector loop, CDT::V2d<double> le_low, CDT::V2d<double> ri_high, double u_tolerance, double v_tolerance, std::vector<point_vector>& small_loops) {
    // 计算左上角边界点和右下角边界点的坐标
    CDT::V2d<double> le_high = {le_low.x, ri_high.y};
    CDT::V2d<double> ri_low = {ri_high.x, le_low.y};
    // 用于判断当前线段是否跨越u或v方向的边界的标记
    bool segSplitU = false;
    bool segSplitV = false;
    // 一个临时二维向量用于存储中间计算结果
    CDT::V2d<double> temp = {0, 0};
    // 用于标记是否存在未知的跨边界情况
    bool has_unknown = false;
    // 两个二维向量用于遍历loop中的点
    CDT::V2d<double> next, curr;

    // icount用于计数，i用于遍历loop中的点
    int icount = 0, i = 0;

    /**
     * 找到第一个跨越边界的线段
     */
    while(i < loop.size()) {
        // 获取下一个点和当前点
        next = loop[(i + 1) % loop.size()];
        curr = loop[i];

        // 计算当前点和下一个点在u(v)方向上是否跨越了边界
        bool cross_u = fabs(next.x - curr.x) > u_tolerance;
        bool cross_v = fabs(next.y - curr.y) > v_tolerance;
        // 当前线段在在u(v)方向上跨越了边界就退出循环
        if(cross_u || cross_v) {
            break;
        }
        i++;
    }
    // 如果没有跨边界的线段
    if(i == loop.size()) {
        // 将整个loop添加到small_loops中
        small_loops.push_back(loop);
    } else {  // 找到了跨边界的线段
        // 用于标记是否是第一次处理
        bool ifFirst = true;

        while(icount < loop.size()) {
            // 用于存储当前分割的线段
            std::vector<CDT::V2d<double>> segment = {};

            if(segSplitU || segSplitV) {
                // 如果跨越u方向上的边界，即需要在u方向上分割
                if(segSplitU) {
                    if(temp.x == le_low.x)
                        temp.x = ri_high.x;
                    else
                        temp.x = le_low.x;
                }
                // 如果跨越v方向上的边界，即需要在v方向上分割
                if(segSplitV) {
                    if(temp.y == le_low.y)
                        temp.y = ri_high.y;
                    else
                        temp.y = le_low.y;
                }
                // 将调整后的点加入segment中
                segment.push_back(temp);
                // 重置两个判断标记
                segSplitU = segSplitV = false;
            }

            while(true) {
                // next和curr就是正在检查的跨越边界的线段的两个端点
                next = loop[(i + 1) % loop.size()];
                curr = loop[i % loop.size()];
                // 递增i，准备检查下一个点
                i++;
                // 如果不是第一次处理
                if(!ifFirst)
                    icount++;
                else
                    ifFirst = false;
                // 将当前点加入segment中，表示当前分割的线段的一部分
                segment.push_back(curr);

                // 计算当前点和下一个点在u(v)方向上是否跨越了边界
                bool cross_u = fabs(next.x - curr.x) > u_tolerance;
                bool cross_v = fabs(next.y - curr.y) > v_tolerance;

                if(cross_u || cross_v) {  // 如果当前两个点由跨界
                    if(fabs(next.x - curr.x) <= u_tolerance && fabs(next.y - curr.y) <= v_tolerance) {
                    } else if(fabs(next.x - curr.x) > u_tolerance && fabs(next.y - curr.y) > v_tolerance) {  // 同时跨越u，v两处边界，分先u后v 和 先v后u两种情况
                        // 定义两个二维向量用于存储分割线段的交点
                        CDT::V2d<double> tempx = {0, 0};
                        CDT::V2d<double> tempy = {0, 0};
                        if(next.x < curr.x && next.y < curr.y) {  // 检查下一个点是否位于当前点的上右方向
                            // 将next的x和y坐标都加上x和y方向的边界宽度
                            next.x += ri_high.x - le_low.x;
                            next.y += ri_high.y - le_low.y;
                            // 在左上和右上边界之间的交点存储在tempx
                            intersect(curr, next, le_high, ri_high, tempx);
                            // 在右下和右上边界之间的交点存储在tempy
                            intersect(curr, next, ri_low, ri_high, tempy);
                        } else if(next.x < curr.x && next.y > curr.y) {  // 检查下一个点是否位于当前点的下右方向
                            // 将next的x坐标加上x方向的边界宽度
                            next.x += ri_high.x - le_low.x;
                            // 将next的y坐标减去y方向的边界宽度
                            next.y -= ri_high.y - le_low.y;
                            // 在左下和右下边界之间的交点存储在tempx
                            intersect(curr, next, le_low, ri_low, tempx);
                            // 在右下和右上边界之间的交点存储在tempy
                            intersect(curr, next, ri_low, ri_high, tempy);
                        } else if(next.x > curr.x && next.y < curr.y) {  // 检查下一个点是否位于当前点的左上方向
                            // 将next的x坐标减去x方向的边界宽度
                            next.x -= ri_high.x - le_low.x;
                            // 将next的y坐标减去y方向的边界宽度
                            next.y += ri_high.y - le_low.y;
                            // 在左上和右上边界之间的交点存储在tempx
                            intersect(curr, next, le_high, ri_high, tempx);
                            // 在左下和左上边界之间的交点存储在tempy
                            intersect(curr, next, le_low, le_high, tempy);
                        } else if(next.x > curr.x && next.y > curr.y) {  // 检查下一个点是否位于当前点的左下方向
                            // 将next的x和y坐标都减去x和y方向的边界宽度
                            next.x -= ri_high.x - le_low.x;
                            next.y -= ri_high.y - le_low.y;
                            // 在左下和右下边界之间的交点存储在tempx
                            intersect(curr, next, le_low, ri_low, tempx);
                            // 在左下和左上边界之间的交点存储在tempy
                            intersect(curr, next, le_low, le_high, tempy);
                        }

                        // 判断哪个交点在x方向上更接近当前点
                        // tempx更接近curr
                        if(fabs(tempx.x - curr.x) < fabs(tempy.x - curr.x)) {
                            temp = tempx;
                            segment.push_back(temp);
                            temp = tempy;
                        } else {  // tempy更接近curr
                            temp = tempy;
                            segment.push_back(temp);
                            temp = tempx;
                        }
                        // 将已经被分割成不跨越边界的线段添加到small_loops中
                        small_loops.push_back(segment);
                        // 设置标记，表示在u和v方向上进行了分割
                        segSplitU = true;
                        segSplitV = true;
                        break;
                    } else if(fabs(next.x - curr.x) > u_tolerance) {  // 仅跨越u边界
                        // 下一个点在当前点的右边
                        if(next.x < curr.x) {
                            // 给next的x坐标增加一个边界宽度，模拟跨越右边界情况
                            next.x += ri_high.x - le_low.x;
                            // 计算当前点和调整后的next与右边界的交点(存储在temp中)
                            bool has_intersect = intersect(curr, next, ri_low, ri_high, temp);
                            // 如果没有交点，设置一个交点
                            if(!has_intersect) {
                                temp.y = curr.y;
                                temp.x = ri_high.x;
                            }
                        } else {  // 下一个点在当前点的左边
                            // 给next的x坐标减少一个边界宽度，模拟跨越左边界情况
                            next.x -= ri_high.x - le_low.x;
                            // 计算当前点和调整后的next与左边界的交点(存储在temp中)
                            bool has_intersect = intersect(curr, next, le_low, le_high, temp);
                            // 如果没有交点，设置一个交点
                            if(!has_intersect) {
                                temp.y = curr.y;
                                temp.x = le_low.x;
                            }
                        }

                        // 设置标记，表示在u方向上进行了分割
                        segSplitU = true;

                    } else if(fabs(next.y - curr.y) > v_tolerance) {  // 仅跨越v边界
                        // 下一个点在当前点的上边
                        if(next.y < curr.y) {
                            // 给next的y坐标增加一个边界高度，模拟跨越上边界情况
                            next.y += ri_high.y - le_low.y;
                            // 计算当前点和调整后的next与上边界的交点(存储在temp中)
                            bool has_intersect = intersect(curr, next, le_high, ri_high, temp);
                            // 如果没有交点，设置一个交点
                            if(!has_intersect) {
                                temp.x = curr.x;
                                temp.y = ri_high.y;
                            }
                        } else {  // 下一个点在当前点的下边
                            // 给next的y坐标减少一个边界高度，模拟跨越下边界情况
                            next.y -= ri_high.y - le_low.y;
                            // 计算当前点和调整后的next与下边界的交点(存储在temp中)
                            bool has_intersect = intersect(curr, next, le_low, ri_low, temp);
                            // 如果没有交点，设置一个交点
                            if(!has_intersect) {
                                temp.x = curr.x;
                                temp.y = le_low.y;
                            }
                        }

                        // 设置标记，表示在v方向上进行了分割
                        segSplitV = true;
                    }

                    // 将交点加入segment中
                    segment.push_back(temp);
                    // 判断当前正在处理的环的类型
                    if(i == 110) {
                        int kk = 1;
                    }
                    small_loop_type type = judge_type(segment, le_low, ri_high);
                    // 用于控制是否将segment添加到最终结果中
                    bool append = true;
                    switch(type) {
                        // 如果当前segment是退化的，即不构成有效的环
                        case small_loop_type::degenerate: {
                            // 不将该segment添加到最终结果中
                            append = false;
                            // if(segment.size() >= 2 && (u_tolerance > 0.0 || v_tolerance > 0.0)) {
                            //     // 边界边的点需要添加到特定与边界相交的loop中以保持水密性，通常只有封闭曲面需要添加
                            //     double low, high;
                            //     boundary_edge_loop_type kind = judge_edge_type(segment, le_low, ri_high, low, high);
                            //     if(kind != not_boundary_edge) {
                            //         boundary_edge_loops.push_back(boundary_edge_loop(segment, kind, low, high));
                            //     }
                            // }
                            break;
                        }
                        // 如果当前segment跨越了u和v两个方向上的边界
                        case small_loop_type::crossUV: {
                            // 构造一个新点(第一个点的x坐标，最后一个点的y坐标)
                            CDT::V2d<double> tempCorner = {segment[0].x, segment[segment.size() - 1].y};
                            segment.push_back(tempCorner);
                            break;
                        }
                        // 如果当前segment跨越了v和u两个方向上的边界
                        case small_loop_type::crossVU: {
                            // 构造一个新点(最后一个点的x坐标，第一个点的y坐标)
                            CDT::V2d<double> tempCorner = {segment[segment.size() - 1].x, segment[0].y};
                            segment.push_back(tempCorner);
                            break;
                        }
                        // 如果当前segment类型未知
                        case small_loop_type::unknown: {
                            has_unknown = true;
                        }
                        case small_loop_type::crossU:
                        case small_loop_type::crossV:
                        default:
                            break;
                    }
                    // 符合条件将环添加进small_loops中
                    if(append) {
                        small_loops.push_back(segment);
                        segment.clear();
                    }
                    break;
                }
            }
        }
        // small_loops[len - 1].push_back(small_loops[0][1]);
        // small_loops.erase(small_loops.begin());
    }
    // 获取不跨边界的环的数量
    size_t len = small_loops.size();
    // printf("total segments: %u\n", small_loops.size());
}

/**
 * @brief 预处理可能跨边界的loop，调整点使其在u或v方向上闭合
 *
 * @param loop 跨边界的loop
 * @param le_low 左下角边界点
 * @param ri_high 右上角边界点
 * @param closed_u 布尔值，表示环是否在u方向上闭合
 * @param closed_v 布尔值，表示环是否在v方向上闭合
 */
void loop_preprocess(point_vector& loop, point le_low, point ri_high, bool closed_u, bool closed_v) {
    // 当前点和下一个点
    point curr, next;
    // 当前点所在区域和下一个点所在区域
    int curr_area, next_area;
    // next_area = ((loop[0].x > (le_low.x - SPAresabs)) + (loop[0].x > (ri_high.x + SPAresabs))) * 3 + ((loop[0].y > (le_low.y - SPAresabs)) + (loop[0].y > (ri_high.y + SPAresabs)));
    // 遍历loop中的每一个点
    for(int i = 0; i < loop.size(); i++) {
        curr = loop[i];
        next = loop[(i + 1) % loop.size()];
        /*
       2 5 8
       1 4 7
       0 3 6
       */
        // if(curr_area == 4 && next_area != 4)  // 跨边界了
        // {
        // }
        // curr_area = next_area;
        // next_area = ((next.x > (le_low.x - SPAresabs)) + (next.x > (ri_high.x + SPAresabs))) * 3 + ((next.y > (le_low.y - SPAresabs)) + (next.y > (ri_high.y + SPAresabs)));

        // switch(curr_area) {
        //     case 4:
        //         break;
        //     case 1:
        //         curr.x += ri_high.x - le_low.x;
        //         break;
        //     case 3:
        //         curr.y += ri_high.y - le_low.y;
        //         break;
        //     case 5:
        //         curr.y -= ri_high.y - le_low.y;
        //         break;
        //     case 7:
        //         curr.x -= ri_high.x - le_low.x;
        //         break;
        //     default:
        //         printf("error");
        //         break;
        // }
        // 使用floor函数将curr.x调整到指定边界范围内
        if(closed_u) curr.x = curr.x - (ri_high.x - le_low.x) * floor((curr.x - le_low.x) / (ri_high.x - le_low.x));
        // 使用floor函数将curr.y调整到指定边界范围内
        if(closed_v) curr.y = curr.y - (ri_high.y - le_low.y) * floor((curr.y - le_low.y) / (ri_high.y - le_low.y));
        // 更新当前点curr
        loop[i] = curr;
    }
}

/**
 * @brief 初始化环面对象及其相关参数
 *
 */
void torus_faceter::init() {
    // 将f对象中的几何形状更新方程的引用转换为环面类型的指针
#ifdef FACETER_USE_ACIS
    tori = (torus*)&f->geometry()->equation_for_update();
    SPAinterval u_interval = f->geometry()->equation().param_range_u();
    SPAinterval v_interval = f->geometry()->equation().param_range_v();
#else
    tori = (torus*)&f->gme_geometry()->gme_equation_for_update();
    // 获取参数u和v的取值范围
    SPAinterval u_interval = f->gme_geometry()->gme_equation().gme_param_range_u();
    SPAinterval v_interval = f->gme_geometry()->gme_equation().gme_param_range_v();
    // 参数u的起始点和结束点
#endif
    u_le = u_interval.start_pt();
    u_ri = u_interval.end_pt();
    // 参数v的起始点和结束点
    v_le = v_interval.start_pt();
    v_ri = v_interval.end_pt();
    // 设置容差值，用于判断是否跨越边界
    u_tolerance = 1e15;
    v_tolerance = 1e15;
    // 如果在u(v)方向上是封闭的，将容差值设为u(v)参数范围的一半
#ifdef FACETER_USE_ACIS
    if(f->geometry()->equation().closed_u()) u_tolerance = tori->param_range_u().length() / 2;  // 判断跨边界所用值
    if(f->geometry()->equation().closed_v()) v_tolerance = tori->param_range_v().length() / 2;
#else
    if(f->gme_geometry()->gme_equation().gme_closed_u()) u_tolerance = tori->gme_param_range_u().gme_length() / 2;  // 判断跨边界所用值
    if(f->gme_geometry()->gme_equation().gme_closed_v()) v_tolerance = tori->gme_param_range_v().gme_length() / 2;
#endif

    // 设置左下角点的坐标
    le_low.x = u_le;
    le_low.y = v_le;
    // 设置右上角点的坐标
    ri_high.x = u_ri;
    ri_high.y = v_ri;

#ifdef FACETER_USE_ACIS
    SPAunit_vector zeroNormal = f->geometry()->equation().eval_normal(SPApar_pos(0, 0));
#else
    // 计算参数(0,0)处的表面法向量
    SPAunit_vector zeroNormal = f->gme_geometry()->gme_equation().gme_eval_normal(SPApar_pos(0, 0));
    // 如果面方向为真，则取反表面法向量
#endif
    if(f->sense()) zeroNormal = -zeroNormal;
    if(zeroNormal != tori->uv_oridir) normalR = REVERSED;
}

/**
 * @brief 通过调整ulen和vlen参数来确定一个给定环面的最优网格密度配置
 *
 * @param f 待处理的环面
 * @param ulen 环面沿着u方向上的网格数
 * @param vlen 环面沿着v方向上的网格数
 * @param nt 法向容差
 * @param st 距离容差
 */
void get_uvlen_torus(FACE* f, int& ulen, int& vlen, double nt, double st) {
#ifdef FACETER_USE_ACIS
    torus* tori = (torus*)&f->geometry()->equation_for_update();
#else
    // 将f对象中的几何形状更新方程的引用转换为环面类型的指针
    torus* tori = (torus*)&f->gme_geometry()->gme_equation_for_update();
    // 获取环面的次半径
#endif
    double mir = tori->minor_radius;
    // 获取环面的主半径
    double mar = tori->major_radius;
    // 设定初始值
    ulen = 2;
    vlen = 4;
    while(true) {
        // printf("====================%d %d==================\n", ulen, vlen);
        ulen *= 2, vlen *= 2;

        double theta0;
        double theta1 = M_PI / ulen;
        // 根据环面的主半径和次半径的比较关系，计算theta0的值
        if(mar > mir) {
            theta0 = 2 * M_PI / ulen;
        } else if(mar > 0) {
            theta0 = 2 * (M_PI - acos(mar / mir)) / ulen;
        } else {
            theta0 = 2 * acos(-mar / mir) / ulen;
        }
        // 环面上的两个三角形的顶点位置
        SPApar_pos* tmp_p0 = new SPApar_pos[3];
        SPApar_pos* tmp_p1 = new SPApar_pos[3];
        tmp_p0[0] = SPApar_pos(0, theta1 / 2);
        tmp_p0[1] = SPApar_pos(0, -theta1 / 2);
        tmp_p0[2] = SPApar_pos(theta0, theta1 / 2);
        tmp_p1[0] = SPApar_pos(theta0, -theta1 / 2);
        tmp_p1[1] = SPApar_pos(0, -theta1 / 2);
        tmp_p1[2] = SPApar_pos(theta0, theta1 / 2);
        // 检查两个三角形是否满足指定的法向容差nt和距离容差st
        if(checkTriangle(f, nt, st, tmp_p0) && checkTriangle(f, nt, st, tmp_p1)) {
            break;
        }
    }

    /**
     * 不断细化ulen和vlen的值，来寻找最优的网格密度配置
     */
    ulen /= 2, vlen /= 2;
    // 每次循环中调整的步长
    int bias_u = ulen / 2, bias_v = vlen / 2;
    while(bias_u >= 1 && bias_v >= 1) {
        // printf("====================%d %d==================\n", ulen, vlen);
        //  printf("====================%d %d==================\n", bias_u, bias_v);
        // 加网格密度
        ulen += bias_u, vlen += bias_v;
        // 调整每次循环中的步长
        bias_u /= 2, bias_v /= 2;

        double theta0;
        double theta1 = M_PI / ulen;
        if(mar > mir) {
            theta0 = 2 * M_PI / ulen;
        } else if(mar > 0) {
            theta0 = 2 * (M_PI - acos(mar / mir)) / ulen;
        } else {
            theta0 = 2 * acos(-mar / mir) / ulen;
        }
        SPApar_pos* tmp_p0 = new SPApar_pos[3];
        SPApar_pos* tmp_p1 = new SPApar_pos[3];
        if(ulen % 2 == 0) {
            tmp_p0[0] = SPApar_pos(0, theta1 / 2);
            tmp_p0[1] = SPApar_pos(0, -theta1 / 2);
            tmp_p0[2] = SPApar_pos(theta0, theta1 / 2);
            tmp_p1[0] = SPApar_pos(theta0, -theta1 / 2);
            tmp_p1[1] = SPApar_pos(0, -theta1 / 2);
            tmp_p1[2] = SPApar_pos(theta0, theta1 / 2);
        } else {
            tmp_p0[0] = SPApar_pos(theta0 / 2, theta1 / 2);
            tmp_p0[1] = SPApar_pos(-theta0 / 2, theta1 / 2);
            tmp_p0[2] = SPApar_pos(-theta0 / 2, -theta1 / 2);
            tmp_p1[0] = SPApar_pos(theta0 / 2, theta1 / 2);
            tmp_p1[1] = SPApar_pos(theta0 / 2, -theta1 / 2);
            tmp_p1[2] = SPApar_pos(-theta0 / 2, -theta1 / 2);
        }
        // 检查两个三角形是否满足指定的法向容差nt和距离容差st
        if(checkTriangle(f, nt, st, tmp_p0) && checkTriangle(f, nt, st, tmp_p1)) {
            if(bias_u == 0) {
                ulen -= 1, vlen -= 2;
            } else {
                ulen -= 2 * bias_u, vlen -= 2 * bias_v;
            }
            continue;
        }
    }
    ulen += 1, vlen += 2;
    // printf("ulen vlen\n");
    // printf("%d %d\n", ulen, vlen);
}

/**
 * @brief 计算u和v方向上的步长
 *
 */
void torus_faceter::decideUVlen() {
    // 返回ulen和vlen表示u和v方向上的分割数量
    get_uvlen_torus(f, ulen, vlen, nt, st);
    // 计算u和v方向上的步长
    du = (u_ri - u_le) / ulen;
    dv = (v_ri - v_le) / vlen;
    // printf("%lf %lf\n", du, dv);
    // printf("====================%d %d==================\n", ulen, vlen);
}

// 暂时没有添加处理奇点
/**
 * @brief 该函数遍历给定环面的所有环，并根据环的类型将其存储在相应的容器中
 * 如果某个环是外环且经过奇点，函数会进行额外处理
 *
 * @return 是否存在外围环。如果存在外围环，返回true；否则返回false
 */
bool torus_faceter::getLoops() {
    // 用于知识是否存在外环
    bool has_periphery_loop = false;
    // 用于存储参数位置
    SPApar_pos param;
#ifdef FACETER_USE_ACIS
    for(LOOP* lp = f->loop(); lp; lp = lp->next()) {
#else
    // 遍历环面的所有环
    for(LOOP* lp = f->gme_loop(); lp; lp = lp->gme_next()) {
        // 用于存储临时的二维点
#endif
        std::vector<CDT::V2d<double>> tempPS;
        // 用于存储环的类型
        loop_type lptype;
        // 获取loop的类型
        api_loop_type(lp, lptype);
        // 用于计数
        unsigned ic = 0;
        // 用于只是外环是否有奇点
        bool periphery_has_singularity = false;
        // 获取组成该loop的所有的点
        // 根据先验知识需要对获取到的点集进行额外处理
        for(COEDGE* coedge = lp->start(); coedge != lp->start() || ic == 0; coedge = coedge->next(), ic++) {
            // 存储边的多线段和参数信息
            SPAposition* polyline;
            double* params;
            // 判断edge的指向与coedge指向是否相反，关系到点添加顺序问题
            REVBIT r = coedge->sense();
            // 用于存储临时的二维点(包含当前边的点的向量)
            std::vector<CDT::V2d<double>> tempps;
            // 存储点的数量
            int nP = 0;
            // 获取当前边的点和参数信息存储在polyline和params中，同时更新nP
            get_facet_edge_points_and_params(coedge->edge(), polyline, params, nP);
#ifdef FACETER_USE_ACIS
            // api_get_facet_edge_points();
#else
            // gme_api_get_facet_edge_points();
#endif
            // 可能造成内存泄漏；暂不解耦
            // 向当前coedge添加参数曲线
            sg_add_pcurve_to_coedge(coedge);
            // 获取coedge的几何信息存储在PC指针中
            PCURVE* PC = coedge->geometry();
            // 该边为null edge
            if(PC == nullptr) {
                // 获取第一个点的参数位置
#ifdef FACETER_USE_ACIS
                param = f->geometry()->equation().param(polyline[0]);
#else
                param = f->gme_geometry()->gme_equation().gme_param(polyline[0]);
#endif
                tempps.push_back(point(param.u, le_low.y));
                tempps.push_back(point(param.u, ri_high.y));
                // 检查param.u是否等于ri_high.x，相等则反转tempps向量中的点顺序
                if(equal(param.u, ri_high.x)) {
                    std::reverse(tempps.begin(), tempps.end());
                }
            } else {  // 该边为非null edge
                // 获取当前coedge的参数曲线方程存储在pc变量
                pcurve pc = PC->equation();
                // 遍历当前边的所有点(polyline中的点)
                for(int i = 0; i < nP; i++) {
                    // 评估当前点polyline[i]在参数曲线pc中的位置
#ifdef FACETER_USE_ACIS
                    param = pc.eval_position(polyline[i], params[i], 1);
                    SPApar_pos param1 = f->geometry()->equation().param(polyline[i]);
#else
                    param = pc.gme_eval_position(polyline[i], params[i], 1);
                    // 获取polyline[i]在几何方程中的参数位置
                    SPApar_pos param1 = f->gme_geometry()->gme_equation().gme_param(polyline[i]);
#endif

                    // printf("param: %lf \t pc:%lf %lf\tf: %lf %lf\n", params[i], param.u, param.v, param1.u, param1.v);
                    // if(tori->degenerate() && lptype == loop_type::loop_u_separation) param = coedge->geometry()->equation().eval_position(polyline[i], 1);  // @todo: 为该接口第二个参数寻找一个合适的参数par
                    // printf("%lf %lf %lf\t%lf %lf\t%lf\n", polyline[i].x(), polyline[i].y(), polyline[i].z(), param.u, param.v, params[i]);
#ifdef FACETER_USE_ACIS
                    if(f->geometry()->equation().closed_u() && equal(param1.u, u_le) && equal(param.u, u_ri)) {
                        param1.u = u_ri;
                    }
                    if(f->geometry()->equation().closed_v() && equal(param1.v, v_le) && equal(param.v, v_ri)) {
#else
                    if(f->gme_geometry()->gme_equation().gme_closed_u() && equal(param1.u, u_le) && equal(param.u, u_ri)) {
                        param1.u = u_ri;
                    }
                    if(f->gme_geometry()->gme_equation().gme_closed_v() && equal(param1.v, v_le) && equal(param.v, v_ri)) {
#endif
                        param1.v = v_ri;
                    }
                    // 创建一个点(param1.u，param1.v)
                    CDT::V2d<double> temp = {param1.u, param1.v};
                    tempps.push_back(temp);
                }
                // 判断边类型
                boundary_edge_loop_type edge_type = judge_edge_type(tempps, le_low, ri_high);
                // 如果是边界边
                if(edge_type != not_boundary_edge) {
                    // 将此边作为boundary_edge_loop对象添加到boundary_edge_loops向量中
                    boundary_edge_loops.push_back(boundary_edge_loop(tempps, edge_type));
                }
                // 打印输出
                // if(r) {
                //     for(int i = nP - 1; i >= 0; i--) {
                //         param = pc.eval_position(polyline[i], params[i], 1);
                //         SPApar_pos param1 = f->geometry()->equation().param(polyline[i]);
                //         printf("param: %lf \t pc:%lf %lf\tf: %lf %lf\n", params[i], param.u, param.v, param1.u, param1.v);
                //     }
                // } else {
                //     for(int i = 0; i < nP; i++) {
                //         param = pc.eval_position(polyline[i], params[i], 1);
                //         SPApar_pos param1 = f->geometry()->equation().param(polyline[i]);
                //         printf("param: %lf \t pc:%lf %lf\tf: %lf %lf\n", params[i], param.u, param.v, param1.u, param1.v);
                //     }
                // }
                // printf("\n");

                // 如果当前边是外环且起始于奇点终止于奇点
#ifdef FACETER_USE_ACIS
                if(lptype == loop_periphery && coedge->starts_at_singularity() && coedge->ends_at_singularity()) {  // 经过奇点的外环
#else
                if(lptype == loop_periphery && coedge->gme_starts_at_singularity() && coedge->gme_ends_at_singularity()) {  // 经过奇点的外环
#endif
                    if(tempps[0].x < tempps[1].x) {
                        // 在tempps的开头插入一个新点(tempps[0].x,le_low.y)
                        tempps.insert(tempps.begin(), point(tempps[0].x, le_low.y));
                        // 在tempps的末尾添加一个新点(最后一个点.x,le_low.y)
                        tempps.push_back(point(tempps[tempps.size() - 1].x, le_low.y));
                    } else {
                        // 在tempps的开头插入一个新点(tempps[0].x,ri_high.y)
                        tempps.insert(tempps.begin(), point(tempps[0].x, ri_high.y));
                        // 在tempps的末尾添加一个新点(最后一个点.x,ri_high.y)
                        tempps.push_back(point(tempps[tempps.size() - 1].x, ri_high.y));
                    }
                    // 标记当前外环经过奇点
                    periphery_has_singularity = true;
                }
            }
            // 如果当前边的方向与coedge的方向相反，则反转tempps向量中点的顺序
            if(r) std::reverse(tempps.begin(), tempps.end());
            // 将tempps向量中的所有点添加到tempPS向量的末尾
            tempPS.insert(tempPS.end(), tempps.begin(), tempps.end());
        }

        // 如果表面法向量为真，则反转tempPS向量中点的顺序
        if(normalR) {
            std::reverse(tempPS.begin(), tempPS.end());
        }
        // loop_preprocess(tempPS, le_low, ri_high, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v());

        switch(lptype) {
            // 如果当前环是孔环
            case loop_type::loop_hole: {
                hole_loops.push_back(tempPS);
                break;
            }
            // 如果当前环是外围环
            case loop_type::loop_periphery: {
                // printf("periphery\n");
                // 检查是否经过奇点
                if(periphery_has_singularity) {
                    // 断言tempPS的大小至少为 4，以确保后续操作不会越界
                    assert(tempPS.size() >= 4);
                    // 遍历当前环的每个点
                    for(int i = 0; i < tempPS.size(); i++) {
                        // 获取当前点
                        point curr = tempPS[i];
#ifdef FACETER_USE_ACIS
                        if(tori->singular_u(curr.x)) {
#else
                        // 检查当前点是否为奇点
                        if(tori->gme_singular_u(curr.x)) {
                            // 获取第i+3个点
#endif
                            point next3 = tempPS[(i + 3) % tempPS.size()];
                            // 如果next3.x也是奇点
                            if(tori->singular_u(next3.x)) {
                                // 获取第i+1和第i+2个点
                                point next1 = tempPS[(i + 1) % tempPS.size()];
                                point next2 = tempPS[(i + 2) % tempPS.size()];
                                if(next1.y < next2.y && curr.y > next3.y || next1.y > next2.y && curr.y < next3.y) {
                                    // 如果i是倒数第二个点，则删除tempPS的最后一个点和第一个点
                                    if(i == tempPS.size() - 2) {
                                        tempPS.pop_back();
                                        tempPS.erase(tempPS.begin());
                                        i--;
                                    } else {  // 否则断言i+3不超过tempPS的大小删除第i+1到第i+3个点
                                        assert((i + 3) <= tempPS.size());
                                        tempPS.erase(tempPS.begin() + i + 1, tempPS.begin() + i + 3);
                                    }
                                }
                            }
                        }
                    }
                }
                periphery = tempPS;
                // 标记存在外环
                has_periphery_loop = true;

                break;
            }
            // 如果当前环是u方向上分离的环
            case loop_type::loop_u_separation: {
                // 检查当前环是否退化且tempPS的大小不等于2
                if(tori->degenerate() && tempPS.size() != 2) {
                    // 用于标记是否存在穿过奇点的情况
                    bool singularity_cross = false;
                    // 遍历环中的每一个点
                    for(int i = 0; i < tempPS.size() - 1; i++) {
                        // 获取当前点和下一个点
                        point curr = tempPS[i];
                        point next = tempPS[i + 1];
                        // 检查从curr.y到next.y是否穿过了边界
                        cross_type cross_u = cross_boundary(curr.y, next.y, v_tolerance);
                        // 如果穿过边界
                        if(cross_u != no_cross) {
                            // 标记存在穿过奇点的情况
                            singularity_cross = true;
                            break;
                        }
                    }
                    // 当环不穿过奇点时
                    if(!singularity_cross) {
                        // 获取点集的起点和终点
                        point start = tempPS[0];
                        point end = tempPS[tempPS.size() - 1];
                        // 检查从end.y到start.y是否穿过了边界
                        cross_type cross_u = cross_boundary(end.y, start.y, v_tolerance);
                        // 如果没穿过边界
                        if(cross_u == no_cross) {
                            if(start.y > end.y) {
                                if(fabs(start.y - ri_high.y) < v_tolerance / 2)
                                    tempPS.push_back(point(start.x, le_low.y));
                                else
                                    tempPS.push_back(point(start.x, ri_high.y));
                            } else if(end.y > start.y) {
                                if(fabs(end.y - ri_high.y) < v_tolerance / 2)
                                    tempPS.push_back(point(start.x, le_low.y));
                                else
                                    tempPS.push_back(point(start.x, ri_high.y));
                            }
                        }
                    }
                }
                // printf("useperation \n");
                // printLoop(tempPS);
                // printf("\n");
                // 将tempPS向量添加到u_loops向量中
                u_loops.push_back(tempPS);
                break;
            }
            // 如果该环是v方向上分离的环
            case loop_type::loop_v_separation: {
                // 将当前点集加入v_loops向量中
                v_loops.push_back(tempPS);
                // printLoop(tempPS);
                // printf("\n");
                break;
            }
            // 如果该环是u和v方向上同时分离的环
            case loop_type::loop_uv_separation: {
                // 将当前点集添加到uv_loops向量中
                uv_loops.push_back(tempPS);
                error("loop uv seperation");
                break;
            }
            default: {
                error("unknown loop");
            }
        }
    }
    // 返回是否存在外环
    return has_periphery_loop;
}

/**
 * @brief 处理环面外环的分割和边界点的处理
 *
 * 该函数依据边界信息和容差值分割外环，并处理在边界上的点和边。
 * 同时添加在边界上需要的附加点，并连接这些点以形成边。
 */
void torus_faceter::peripheryProcess() {
    // 获取几何方程的u方向是否闭合
#ifdef FACETER_USE_ACIS
    bool closed_u = f->geometry()->equation().closed_u();
    bool closed_v = f->geometry()->equation().closed_v();
#else
    bool closed_u = f->gme_geometry()->gme_equation().gme_closed_u();
    // 获取几何方程的v方向是否闭合
    bool closed_v = f->gme_geometry()->gme_equation().gme_closed_v();
#endif

    point curr, next;
    // 对外围环periphery进行预处理
    preprocessor(periphery, closed_u, closed_v, le_low, ri_high);
    // printLoop(periphery);
    // 用于存储分割后的小环
    std::vector<point_vector> small_loops;
    // 调用函数对外围环按照边界信息和容差值进行分割
    boundary_loop_split(periphery, le_low, ri_high, u_tolerance, v_tolerance, small_loops);

    // 需要使用coedge存储boundary_edge_loops
    if(boundary_edge_loops.size() > 0)
        // 边界边处理

        // 遍历分割出来的小loop
        for(int j = 0; j < small_loops.size(); j++) {
            // 初始化计数器k为0
            int k = 0;
            // 将当前小环赋值给loop变量
            point_vector loop = small_loops[j];
            // 遍历小环的每个点
            while(k < loop.size()) {
                curr = loop[k];
                next = loop[(k + 1) % loop.size()];
                // 判断前后两点是否为在u方向上的边界边
                if(u_tolerance > 0.0 && ((equal(curr.y, le_low.y) && equal(next.y, le_low.y)) || (equal(curr.y, ri_high.y) && equal(next.y, ri_high.y)))) {
                    // 遍历每个边界边
                    for(int i = 0; i < boundary_edge_loops.size(); i++) {
                        if(boundary_edge_loops[i].kind == be_vle || boundary_edge_loops[i].kind == be_vri) {
                            // 获取当前边的起点和终点
                            point start = boundary_edge_loops[i].edge_loop[0];
                            point end = boundary_edge_loops[i].edge_loop[boundary_edge_loops[i].edge_loop.size() - 1];
                            // 根据不同的情况进行插入边界边操作
                            if(curr.x < start.x && start.x < end.x && end.x < next.x) {
                                loop.insert(loop.begin() + k, boundary_edge_loops[i].edge_loop.begin(), boundary_edge_loops[i].edge_loop.end());
                                k += boundary_edge_loops[i].edge_loop.size();
                            } else if(curr.x < end.x && end.x < start.x && start.x < next.x) {
                                loop.insert(loop.begin() + k, boundary_edge_loops[i].edge_loop.rbegin(), boundary_edge_loops[i].edge_loop.rend());
                                k += boundary_edge_loops[i].edge_loop.size();
                            } else if(next.x < end.x && end.x < start.x && start.x < curr.x) {
                                loop.insert(loop.begin() + k, boundary_edge_loops[i].edge_loop.begin(), boundary_edge_loops[i].edge_loop.end());
                                k += boundary_edge_loops[i].edge_loop.size();
                            } else if(next.x < start.x && start.x < end.x && end.x < curr.x) {
                                loop.insert(loop.begin() + k, boundary_edge_loops[i].edge_loop.rbegin(), boundary_edge_loops[i].edge_loop.rend());
                                k += boundary_edge_loops[i].edge_loop.size();
                            }
                        }
                    }
                    // 判断前后两点是否为在v方向上的边界边
                } else if(v_tolerance > 0.0 && ((equal(curr.x, le_low.x) && equal(next.x, le_low.x)) || (equal(curr.x, ri_high.x) && equal(next.x, ri_high.x)))) {
                }
                k++;
            }
            // 更新small_loops[j]为处理后的loop
            small_loops[j] = loop;
        }

    // 遍历分割出来的小环
    for(int j = 0; j < small_loops.size(); j++) {
        // 将当前小环赋值给loop变量
        point_vector loop = small_loops[j];
        // 标记当前环在points向量中的起始位置
        size_t start = points.size();
        // 边界点结尾，用于标记环是否以边界点结束
        bool boundary_point_end = false;
        // 遍历loop中的每个点
        for(int p = 0; p < loop.size(); p++) {
            // 将当前点添加到points向量中
            points.push_back(loop[p]);
            // 获取当前点和下一个点
            curr = loop[p];
            next = loop[(p + 1) % loop.size()];
            if(next == curr) {  // 重复点需要丢弃
                points.pop_back();
                continue;
            }
            // 将points向量的最后一个点和新点连接成边加入edges向量中
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            // 如果两点位于边界，需要添加在面内的点
            if((fabs(curr.x - u_ri) < EPSILON && fabs(next.x - u_ri) < EPSILON) || (fabs(curr.x - u_le) < EPSILON && fabs(next.x - u_le) < EPSILON)) {
                // 如果两点在y方向上的距离小于v方向的步长，不做处理
                if(fabs(curr.y - next.y) < dv) continue;
                // printf("%lf %lf %lf %lf\n", curr.x, curr.y, next.x, next.y);
                if(curr.y < next.y) {
                    // 调用函数获取接近v_le的点，以步长dv遍历，直到到达next.y的上方一点(减去EPSILON)
                    for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low < next.y - EPSILON; low += dv) {
                        if(low < curr.y) continue;
                        // 将点(curr.x, low)加入points向量中
                        CDT::V2d<double> temp = {curr.x, low};
                        points.push_back(temp);
                        // 将points向量的最后一个点和新点连接成边加入edges向量中
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                } else {
                    // 调用函数获取接近v_le的点，以步长dv遍历，直到low的值小于next.y的上方一点(加上EPSILON)
                    for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low > next.y + EPSILON; low -= dv) {
                        if(low > curr.y) continue;
                        // 将点(curr.x, low)加入points向量中
                        CDT::V2d<double> temp = {curr.x, low};
                        points.push_back(temp);
                        // 将points向量的最后一个点和新点连接成边加入edges向量中
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                }
            } else if((fabs(curr.y - v_ri) < EPSILON && fabs(next.y - v_ri) < EPSILON) || (fabs(curr.y - v_le) < EPSILON && fabs(next.y - v_le) < EPSILON)) {  // 如果两点位于v_ri或v_le边界的情况
                // 如果两点在x方向上的距离小于u方向的步长，不做处理
                if(fabs(curr.x - next.x) < du) continue;
                if(curr.x < next.x) {
                    // 调用函数获取接近u_le的点，以步长du遍历，直到low的值小于next.x
                    for(double low = get_cloest_boundary_point(u_le, du, curr.x); low < next.x; low += du) {
                        if(low < curr.x) continue;
                        // 将点(low, curr.y)加入points向量中
                        CDT::V2d<double> temp = {low, curr.y};
                        points.push_back(temp);
                        // 将points向量的最后一个点和新点连接成边加入edges向量中
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                } else {
                    // 调用函数获取接近u_le的点，以步长du遍历，直到low的值大于next.x
                    for(double low = get_cloest_boundary_point(u_le, du, curr.x); low > next.x; low -= du) {
                        if(low > curr.x) continue;
                        // 将点(low, curr.y)加入points向量中
                        CDT::V2d<double> temp = {low, curr.y};
                        points.push_back(temp);
                        // 将points向量的最后一个点和新点连接成边加入edges向量中
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                }
            }
        }
        // 从edges向量中移除最后一个边
        edges.pop_back();
        // 连接points向量的最后一个点和起始点成一个新的边添加到edges向量中
        edges.push_back(CDT::Edge(points.size() - 1, start));
    }
}

/**
 * @brief 处理u分离环，生成网格结构并连接边界
 *
 * @return 返回false
 */
bool torus_faceter::USeperationProcess() {
    // 变量param用于存储参数位置
    SPApar_pos param;
    // 存储当前点和下一个点
    point curr, next;
    // 获取左上角和右下角的坐标
    point le_high(le_low.x, ri_high.y);
    point ri_low(ri_high.x, le_low.y);

    // 断言u_loops的大小必须为 2
    assert(u_loops.size() == 2);

    // 记录跨边界的点
    std::vector<unsigned> u_points;
    // 遍历u_loops中的每一个环
    for(unsigned i = 0; i < u_loops.size(); i++) {
        // 将当前的环赋值给 u_loop
        std::vector<CDT::V2d<double>> u_loop = u_loops[i];
        // 对当前环进行预处理
#ifdef FACETER_USE_ACIS
        preprocessor(u_loop, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v(), le_low, ri_high);
#else
        preprocessor(u_loop, f->gme_geometry()->gme_equation().gme_closed_u(), f->gme_geometry()->gme_equation().gme_closed_v(), le_low, ri_high);
#endif
        // printLoop(u_loop);
        // preprocessor(u_loop, closed_u, closed_v, le_low, ri_high);
        // printf("\n");
        //
        int istart = 0;
        // 将当前环的第一个点赋值给curr，然后istart自增
        curr = u_loop[istart++];

        unsigned index = 0;
        // 试图取到第一个u值不为-PI的点，来判断PI该取正值还是负值
        // bool PI_in_right = false;
        // while(fabs(curr.x - u_le) <= SPAresabs && istart < u_loop.size()) {
        //     curr = u_loop[istart++];
        // }
        // if(istart == u_loop.size()) {
        //     curr.x += 2 * SPAresabs;
        //     param.u = curr.x;
        //     param.v = curr.y;
        // }
        // if(curr.x > 0) {
        //     PI_in_right = true;
        // }

        // 遍历当前环的点
        for(int i = 0; i < u_loop.size(); i++) {
            // 如果是环中的第一个点。获取points向量的当前大小(即已存储点的数量)，用于记录这个环的起始点
            if(i == 0) istart = points.size();
            // 获取当前点和下一个点
            curr = u_loop[i];
            next = u_loop[(i + 1) % u_loop.size()];
            // if(PI_in_right && fabs(curr.x - u_le) <= SPAresabs) {  // 根据pi所处的位置更改方向
            //     curr.x = u_ri;
            // }

            // 将当前点添加到points向量中并获取points向量的新大小
            points.push_back(curr);
            index = points.size();
            // 连接最后一个点和当前点成边添加到edges向量中
            edges.push_back(CDT::Edge(index - 1, index));
            // 仅有两个点的情况，通常是奇点，为了避免重复计算交点，选择跳过
            if(u_loop.size() == 2 && i == 1) {
                continue;
            }

            // 判断当前点和下一个点的y坐标之差是否大于v容差
            if(fabs(next.y - curr.y) > v_tolerance) {
                // printf("%lf %lf\n", curr.x, curr.y);
                CDT::V2d<double> temp = {0, 0};
                // 用于记录跨边界的点的索引
                u_points.push_back(index);
                if(next.y > curr.y) {  // 下穿上，指向下
                    // 使next进入下方区域
                    next.y -= (v_ri - v_le);
                    if(curr == next) next.y += 10;
                    // 检查curr和next之间的线段是否与边界线相交，交点存储在temp中
                    bool has_intersect = intersect(curr, next, le_low, ri_low, temp);
                    // printf("%d\n", has_intersect);
                    // 如果没有相交
                    if(!has_intersect) {
                        temp.x = curr.x;
                        temp.y = v_le;
                    }
                    // printf("%d\n", has_intersect);
                    points.push_back(temp);
                    temp.y = v_ri;
                    points.push_back(temp);
                } else {  // 上穿下，指向上
                    // 使next进入上方区域
                    next.y += (v_ri - v_le);
                    if(curr == next) next.y += 10;
                    // 检查curr和next之间的线段是否与边界线相交，交点存储在temp中
                    bool has_intersect = intersect(curr, next, le_high, ri_high, temp);
                    // 如果没有相交
                    if(!has_intersect) {
                        temp.x = curr.x;
                        temp.y = v_ri;
                    }
                    // printf("%d\n", has_intersect);
                    points.push_back(temp);
                    temp.y = v_le;
                    points.push_back(temp);
                }
                // 将index+1和index+2连接成一条边并添加到edges向量中
                edges.push_back(CDT::Edge(index + 1, index + 2));
            }
        }
        // 移除最后一个边
        edges.pop_back();
        index = points.size();
        // 将最后一个点和初始点连成边，完成环的闭合
        edges.push_back(CDT::Edge(index - 1, istart));
    }
    // 必须先经过排序，确保点的顺序正确
    // assert(u_points.size() == 2);
    // printf("%d\n", u_points.size());
    // 对u_points向量中的点按照x坐标进行排序
    for(int i = 0; i < u_points.size(); i++) {
        unsigned ind = u_points[i];
        for(int j = i + 1; j < u_points.size(); j++) {
            unsigned ind2 = u_points[j];
            if(points[ind].x > points[ind2].x) {
                unsigned tt = u_points[i];
                u_points[i] = u_points[j];
                u_points[j] = tt;
                break;
            }
        }
    }

    for(int i = 0; i < u_points.size(); i++) {
        curr = points[u_points[i]];
        next = points[u_points[i] + 1];
        if(next.y > curr.y) {  // 指向下
            // 找到下一个指向上的点或者边界
            // 如果找不到
            if(i == u_points.size() - 1) {
                // 创建两个新的点加入points向量中
                CDT::V2d<double> temp = {u_ri, v_le};
                points.push_back(temp);
                temp.y = v_ri;
                points.push_back(temp);

                // 连接当前处理点和新添加的点(u_ri, v_le)
                edges.push_back(CDT::Edge(u_points[i], points.size() - 2));
                // 连接新添加的两个点(u_ri, v_le)和(u_ri, v_ri)
                edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));
                // 连接新添加的点(u_ri, v_ri)与当前处理点的下一个点
                edges.push_back(CDT::Edge(points.size() - 1, u_points[i] + 1));
            } else {
                // u_points[i]是points向量中的索引
                edges.push_back(CDT::Edge(u_points[i], u_points[i + 1] + 1));
                edges.push_back(CDT::Edge(u_points[i] + 1, u_points[i + 1]));
            }
        }
        // 若当前点是在边界的底部，并且下一个点是在边界的顶部。
        if(i == 0 && next.y < curr.y) {
            // 创建两个新的点加入points向量中
            CDT::V2d<double> temp = {u_le, v_ri};
            points.push_back(temp);
            temp.y = v_le;
            points.push_back(temp);

            // 连接当前点与第一个新点之间的边
            edges.push_back(CDT::Edge(u_points[i], points.size() - 2));
            // 连接第一个新点与第二个新点之间的边
            edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));
            // 连接第二个新点与当前点的下一个点之间的边
            edges.push_back(CDT::Edge(points.size() - 1, u_points[i] + 1));
        }
    }
    return false;
}

/**
 * @brief 将vseperation转化为periphery以进行periphery的处理
 *
 * @return 返回true，表示成功处理
 */
bool torus_faceter::VSeperationProcess() {
    // 参数位置
    SPApar_pos param;
    // 当前点和下一个点
    point curr, next;
    // 获取左上角点和右下角点的位置
    point le_high(le_low.x, ri_high.y);
    point ri_low(ri_high.x, le_low.y);

    // 断言v_loops的大小必须为 2
    assert(v_loops.size() == 2);

    // 记录跨域点的索引及其y值的向量
    std::vector<std::vector<std::pair<unsigned, double>>> v_points;  // 记录跨域点的索引，及其y值，并且索引0记录y最小值，索引1记录y最大值

    // 遍历v_loops中的每一个环
    for(unsigned j = 0; j < v_loops.size(); j++) {
        // 将当前环赋值给v_loop
        std::vector<CDT::V2d<double>> v_loop = v_loops[j];
        // 记录当前vloop的跨域点
        std::vector<std::pair<unsigned, double>> v_point;

        // 进行前处理，确保位于参数域边界的边不会反复横跳
#ifdef FACETER_USE_ACIS
        preprocessor(v_loop, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v(), le_low, ri_high);
#else
        preprocessor(v_loop, f->gme_geometry()->gme_equation().gme_closed_u(), f->gme_geometry()->gme_equation().gme_closed_v(), le_low, ri_high);
#endif
        // printLoop(v_loop);
        // printf("\n");

        // 遍历寻找跨域点
        for(int i = 0; i < v_loop.size(); i++) {
            // 获取当前点和下一个点
            curr = v_loop[i];
            next = v_loop[(i + 1) % v_loop.size()];

            // 如果下一个点和当前点的x坐标差值大于u容差
            if(fabs(next.x - curr.x) > u_tolerance) {
                CDT::V2d<double> temp = {0, 0};  // 临时点

                if(next.x > curr.x) {       // 左穿右，指向左
                    next.x -= u_ri - u_le;  // 调整next点的x坐标使其进入左侧区域
                    // 获取交点
                    bool has_intersect = intersect(next, curr, le_low, le_high, temp);
                    if(!has_intersect) {
                        temp.x = u_le;
                        temp.y = curr.y;
                    }
                    // printf("%d\n", has_intersect);
                    // 将交点temp插入到v_loop中的当前点后面
                    v_loop.insert(v_loop.begin() + i + 1, temp);

                    // 仅记录最大和最小的即可
                    if(v_point.size() == 0) {
                        v_point.push_back(std::make_pair(i + 1, temp.y));  // 记录y最小值
                        v_point.push_back(std::make_pair(i + 1, temp.y));  // 记录y最大值
                    } else {
                        // 如果temp的y值比当前记录的最小值还小，则更新最小值的索引和值
                        if(temp.y < v_point[0].second) {
                            v_point[0].first = i + 1;
                            v_point[0].second = temp.y;
                        } else if(temp.y > v_point[1].second) {  // 如果temp的y值比当前记录的最大值还大，则更新最大值的索引和值
                            v_point[1].first = i + 1;
                            v_point[1].second = temp.y;
                        }
                    }

                    temp.x = u_ri;
                    // 将右侧边界点temp插入到v_loop中的当前点后面
                    v_loop.insert(v_loop.begin() + i + 2, temp);
                    // 跳过插入的两个点
                    i += 2;
                } else {  // 右穿左，指向右
                    // 调整next点的x坐标使其进入右侧区域
                    next.x += u_ri - u_le;
                    // 检查curr和next之间的线段是否与右侧边界线相交，交点存储在temp中
                    bool has_intersect = intersect(next, curr, ri_low, ri_high, temp);
                    if(!has_intersect) {
                        temp.x = u_ri;
                        temp.y = curr.y;
                    }
                    // 将交点temp插入到v_loop中的当前点后面
                    v_loop.insert(v_loop.begin() + i + 1, temp);

                    // 仅记录最大和最小的即可
                    if(v_point.size() == 0) {
                        v_point.push_back(std::make_pair(i + 1, temp.y));  // 记录y最小值
                        v_point.push_back(std::make_pair(i + 1, temp.y));  // 记录y最大值
                    } else {
                        if(temp.y < v_point[0].second) {
                            v_point[0].first = i + 1;
                            v_point[0].second = temp.y;
                        } else if(temp.y > v_point[1].second) {
                            v_point[1].first = i + 1;
                            v_point[1].second = temp.y;
                        }
                    }
                    temp.x = u_le;
                    // 将左侧边界点temp插入到v_loop中的当前点后面
                    v_loop.insert(v_loop.begin() + i + 2, temp);
                    // 跳过插入的两个点
                    i += 2;
                }
            }
        }
        // 更新v_loops中的第j个环
        v_loops[j] = v_loop;
        // 将当前环的跨域点信息v_point添加到v_points中
        v_points.push_back(v_point);
    }

    // 如果第一个vloop不是指向左边的，则翻转之
    if(v_loops[0][v_points[0][0].first].x > v_loops[0][v_points[0][0].first + 1].x) {
        std::reverse(v_loops.begin(), v_loops.end());
        std::reverse(v_points.begin(), v_points.end());
    }

    // if(equal(v_points[0][1].second, v_points[1][0].second)) {
    //     if(v_points[0][1].second < v_points[1][0].second) return true;
    //     return false;
    // } else
    if(v_points[0][1].second < v_points[1][0].second) {
        // 将第一个vloop的所有点复制给periphery
        periphery = v_loops[0];
        point_vector temp_vector;
        // 添加一些点防止后续处理认为其跨越边界
        cross_type left_right_cross = cross_boundary(v_points[0][1].second, v_points[1][0].second, v_tolerance);  // 检查左右跨界类型
        // 如果存在左右跨越边界，在temp_vector中添加一个点，用于模拟跨越边界
        if(left_right_cross != no_cross) temp_vector.push_back(point(u_le, (v_points[0][1].second + v_points[1][0].second) / 2));
        temp_vector.insert(temp_vector.end(), v_loops[1].begin() + v_points[1][0].first + 1, v_loops[1].end());
        temp_vector.insert(temp_vector.end(), v_loops[1].begin(), v_loops[1].begin() + v_points[1][0].first + 1);
        if(left_right_cross != no_cross) temp_vector.push_back(point(u_ri, (v_points[0][1].second + v_points[1][0].second) / 2));
        // 将temp_vector中的所有点插入到periphery中的第v_points[0][1].first + 1的位置之后
        periphery.insert(periphery.begin() + v_points[0][1].first + 1, temp_vector.begin(), temp_vector.end());
    } else if(v_points[0][1].second > v_points[1][0].second) {
        periphery = v_loops[0];
        point_vector temp_vector;
        // 添加一些点，防止后续认为其与上界或下界相连需要跨越边界
        // 指向左边的loop是否会判断为跨边界
        cross_type left_cross = cross_boundary(v_points[0][1].second, v_ri, v_tolerance);

        // 检查指向右边的loop是否会跨越边界
        cross_type right_cross = cross_boundary(v_points[1][0].second, v_le, v_tolerance);
        // 如果存在左跨越边界，在temp_vector中添加一个点，用于模拟跨越边界
        if(left_cross != no_cross) temp_vector.push_back(point(u_le, (v_ri + v_points[0][1].second) / 2));
        // if(v_points[0][1].second < )
        temp_vector.push_back(le_high);
        temp_vector.push_back(le_low);
        // 如果存在右跨越边界，在temp_vector中添加一个点，用于模拟跨越边界
        if(right_cross != no_cross) temp_vector.push_back(point(u_le, (v_le + v_points[1][0].second) / 2));

        temp_vector.insert(temp_vector.end(), v_loops[1].begin() + v_points[1][0].first + 1, v_loops[1].end());
        temp_vector.insert(temp_vector.end(), v_loops[1].begin(), v_loops[1].begin() + v_points[1][0].first + 1);

        if(right_cross != no_cross) temp_vector.push_back(point(u_ri, (v_le + v_points[1][0].second) / 2));
        temp_vector.push_back(ri_low);
        temp_vector.push_back(ri_high);

        if(left_cross != no_cross) temp_vector.push_back(point(u_ri, (v_ri + v_points[0][1].second) / 2));

        // 将temp_vector中的所有点插入到periphery中的第v_points[0][1].first + 1的位置之后
        periphery.insert(periphery.begin() + v_points[0][1].first + 1, temp_vector.begin(), temp_vector.end());
    }
    return true;
    // if(isSame) {
    //     printf("same\n");
    //     v_points.pop_back();
    //     v_points.pop_back();
    //     // 设置一整个域为空洞
    //     points.push_back(le_low);
    //     points.push_back(CDT::V2d<double>(le_low.x, ri_high.y));
    //     points.push_back(ri_high);
    //     points.push_back(CDT::V2d<double>(ri_high.x, le_low.y));
    //     edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));
    //     edges.push_back(CDT::Edge(points.size() - 2, points.size() - 3));
    //     edges.push_back(CDT::Edge(points.size() - 3, points.size() - 4));
    //     edges.push_back(CDT::Edge(points.size() - 1, points.size() - 4));
    // }
    // for(int i = 0; i < v_points.size(); i++) {
    //     curr = points[v_points[i]];
    //     next = points[v_points[i] + 1];
    //     if(next.x < curr.x) {  // 指向右
    //                            // 找到下一个指向左的点或者边界
    //                            // 如果找不到
    //         if(i == v_points.size() - 1) {
    //             points.push_back(ri_high);
    //             points.push_back(le_high);

    // edges.push_back(CDT::Edge(v_points[i], points.size() - 2));
    // edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));
    // edges.push_back(CDT::Edge(points.size() - 1, v_points[i] + 1));
    // } else {
    // edges.push_back(CDT::Edge(v_points[i], v_points[i + 1] + 1));
    // edges.push_back(CDT::Edge(v_points[i] + 1, v_points[i + 1]));
    // }
    // }
    // if(i == 0 && next.x > curr.x) {
    // points.push_back(le_low);
    // points.push_back(ri_low);

    // edges.push_back(CDT::Edge(v_points[i], points.size() - 2));
    // edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));
    // edges.push_back(CDT::Edge(points.size() - 1, v_points[i] + 1));
    // }
    // }
}

/**
 * @brief 根据给定边界和容差将孔环分割成小环，并将小环添加到顶点和边列表中
 *
 */
void torus_faceter::holeProcess() {
    // 遍历每一个孔环
    for(int i = 0; i < hole_loops.size(); i++) {
        std::vector<point_vector> small_loops;
        // 将当前孔环按照给定的边界和容差分割为小环存储在small_loops中
        boundary_loop_split(hole_loops[i], le_low, ri_high, u_tolerance, v_tolerance, small_loops);
        // 遍历每一个小环
        for(int j = 0; j < small_loops.size(); j++) {
            point_vector loop = small_loops[j];
            // 记录当前points向量的大小，作为新环的起始索引
            size_t start = points.size();
            // 遍历小环的每一个点
            for(int p = 0; p < loop.size(); p++) {
                points.push_back(loop[p]);
                // 连接新加入的点的前一个点成边加入edges向量
                edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            }
            // 之前最后一个点与第一个点形成的边是多余的，要移除
            edges.pop_back();
            // 重新连接最后一个点和起始点的边形成闭合环形
            edges.push_back(CDT::Edge(points.size() - 1, start));
        }
    }
}

/**
 * 存储扫描线与边界交点的信息
 */
struct scan_inter {
    double u;       // 扫描线与边界交点的 u 坐标值
    int low_inter;  // 扫描线与边界的交点中位于扫描线下方的交点数量
    int up_inter;   // 扫描线与边界的交点中位于扫描线上方的交点数量
};

/**
 * @brief 根据给定的参数和边界条件在环面上进行扫描和交点计算，然后使用CDT进行三角剖分，生成索引化的网格，并附加到指定的面对象上
 *
 */
void torus_faceter::attachMesh() {
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
        SPAunit_vector dir(1.0, 0.0, 0.0);                         // x轴正方向
        straight v_line(a, dir);                                   // 从点a出发、方向为dir的直线
        SPAinterval line_range(0.0, u_ri - u_le + 4 * SPAresabs);  // 创建直线v_line的范围
#ifdef FACETER_USE_ACIS
        v_line.limit(line_range);
#else
        v_line.gme_limit(line_range);  // 将直线的范围限制为line_range
#endif
        // 用于存储扫描线与曲面的交点信息
        std::vector<scan_inter> inter_u;

        for(unsigned i = 0; i < edges.size(); i++) {  // 遍历edges向量中的每条边
            // 获取当前边的起始点和终点
            point curr = points[edges[i].v1()];
            point next = points[edges[i].v2()];
            // 跳过长度为零的边
            if(curr == next) continue;

            SPAposition b(curr.x, curr.y, 0.0);
            SPAposition c(next.x, next.y, 0.0);

            SPAvector edge_dir = c - b;  // 从b到c的方向向量
#ifdef FACETER_USE_ACIS
            straight edge_line(b, normalise(edge_dir));
            SPAinterval edge_range(0.0, (edge_dir).len());
            edge_line.limit(edge_range);
            curve_curve_int* results_head = int_cur_cur(v_line, edge_line);
#else
            straight edge_line(b, gme_normalise(edge_dir));  // 从b出发、方向为edge_dir的直线
            SPAinterval edge_range(0.0, (edge_dir).len());
            edge_line.gme_limit(edge_range);  // 将直线的范围限制为edge_range
            // 计算扫描线v_line与边edge_line的曲线的交点
            curve_curve_int* results_head = gme_int_cur_cur(v_line, edge_line);
            // 指针用于遍历保存的交点链表
#endif
            curve_curve_int* results = results_head;

            // 检查是否有曲线与曲线的交点
            if(results != nullptr) {
                // 如果有多个交点，就删除曲线与曲线的交点链表，继续下一个边的处理
                if(results->next != nullptr) {
#ifdef FACETER_USE_ACIS
                    delete_curve_curve_ints(results_head);
                    continue;
                }
                double p1 = results->param1;
                SPAposition inters = v_line.eval_position(p1);
#else
                    gme_delete_curve_curve_ints(results_head);
                    continue;
                }
                double p1 = results->param1;                        // 获取第一个交点的参数值
                SPAposition inters = v_line.gme_eval_position(p1);  // 获取交点位置
#endif
                unsigned j = 0;
                for(; j < inter_u.size(); j++) {
                    if(inter_u[j].u == inters.x()) {
                        break;
                    }
                }
                if(j == inter_u.size()) {
                    // 将一个新的scan_inter结构体对象添加到inter_u向量中
                    inter_u.push_back(scan_inter(inters.x(), 0, 0));
                }
                // 如果交点在边界上方或下方
                if(equal(v, curr.y) && equal(v, next.y)) {
                    inter_u[j].low_inter += 1;                                                 // 位于下方的交点数量增加
                    inter_u[j].up_inter += 1;                                                  // 位于上方的交点数量增加
                } else if(equal(v, curr.y) && v > next.y || v > curr.y && equal(v, next.y)) {  // 如果交点在边界上方
                    inter_u[j].low_inter += 1;
                } else if(equal(v, curr.y) && v < next.y || v < curr.y && equal(v, next.y)) {  // 如果交点在边界下方
                    inter_u[j].up_inter += 1;
                } else {  // 交点同时在边界的上方和下方
                    inter_u[j].low_inter += 1;
                    inter_u[j].up_inter += 1;
                }
            }
            // 释放results_head指向的曲线与曲线的交点链表内存

#ifdef FACETER_USE_ACIS
            delete_curve_curve_ints(results_head);
#else
            gme_delete_curve_curve_ints(results_head);
#endif
        }

        // 对结构体对象按照u值进行升序排序
        std::sort(inter_u.begin(), inter_u.end(), [](const scan_inter& a, const scan_inter& b) { return a.u < b.u; });
        // 遍历每一个scan_inter
        for(unsigned i = 0; i < inter_u.size(); i++) {
            if((inter_u[i].up_inter + inter_u[i].low_inter) % 2 == 1) {  // 奇数次交点情况
                // 第一个或最后一个scan_inter对象
                if(i == 0 || i == inter_u.size() - 1) {
                    // 移除该奇数次交点的scan_inter对象
                    inter_u.erase(inter_u.begin() + i);
                    i -= 1;
                    continue;
                }
            }
            if(inter_u[i].low_inter % 2 == 0 && inter_u[i].up_inter % 2 == 0) {  // 偶数次交点情况
                // 移除该偶数次交点的 scan_inter 对象
                inter_u.erase(inter_u.begin() + i);
                i -= 1;
            }
        }

        // if (inter_u.size() % 2 != 0)
        //     assert(inter_u.size() % 2 == 0);
        // 在指定的扫描线上生成新的点坐标，并将这些点坐标添加到points向量中
        for(int i = 0; i < inter_u.size(); i += 2) {
            double inf = inter_u[i].u;  // 获取当前段的起始u坐标
            if(i + 1 >= inter_u.size()) continue;
            double sup = inter_u[i + 1].u;  // 获取当前段的结束u坐标
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
    pcdt = new CDT::Triangulation<double>;              // 新建一个三角剖分对象。
    CDT::RemoveDuplicatesAndRemapEdges(points, edges);  // 处理点和边的重复，重新映射边的索引

    // 将处理过的点集和边集插入pcdt对象中
    pcdt->insertVertices(points);
    pcdt->insertEdges(edges);
    pcdt->eraseOuterTrianglesAndHoles();  // 清除外部三角形和孔洞
    //} while(_count++ < 6 && !fixTriangle(f, nt, st, pcdt, &points));

    CDT::TriangleVec triangles = pcdt->triangles;                                                                      // 获取三角形的集合
    std::vector<CDT::V2d<double>> vertice = pcdt->vertices;                                                            // 获取顶点的集合
    INDEXED_MESH* mesh = new INDEXED_MESH("gme", vertice.size(), pcdt->triangles.size(), pcdt->triangles.size() * 3);  // 创建一个索引化的网格对象mesh
    // 遍历每个顶点
    for(int i = 0; i < vertice.size(); i++) {
        CDT::V2d<double> pp = vertice[i];
        SPApar_pos param(pp.x, pp.y);  // 获取参数位置
#ifdef FACETER_USE_ACIS
        SPAposition pos = f->geometry()->equation().eval_position(param);
        SPAunit_vector normal = f->geometry()->equation().eval_normal(param);
        if(tori->degenerate() && fabs(param.u - u_ri) < SPAresabs) {
            param.u -= 2 * SPAresabs;
            normal = f->geometry()->equation().eval_normal(param);
            param.u += 2 * SPAresabs;
        } else if(tori->degenerate() && fabs(param.u - u_le) < SPAresabs) {
            param.u += 2 * SPAresabs;
            normal = f->geometry()->equation().eval_normal(param);
#else
        SPAposition pos = f->gme_geometry()->gme_equation().eval_position(param);          // 计算几何位置
        SPAunit_vector normal = f->gme_geometry()->gme_equation().gme_eval_normal(param);  // 计算法向量
        if(tori->degenerate() && fabs(param.u - u_ri) < SPAresabs) {                       // 如果是退化状态且u接近u_ri或u_le
            param.u -= 2 * SPAresabs;
            normal = f->gme_geometry()->gme_equation().gme_eval_normal(param);
            param.u += 2 * SPAresabs;
        } else if(tori->degenerate() && fabs(param.u - u_le) < SPAresabs) {
            param.u += 2 * SPAresabs;
            normal = f->gme_geometry()->gme_equation().gme_eval_normal(param);
#endif
            param.u -= 2 * SPAresabs;
        }
        if(r) normal = -normal;
        // 将顶点信息添加到mesh对象中
        mesh->add_vertex("gme", pos, normal, param);
    }

    // 遍历每个三角形
    for(int i = 0; i < triangles.size(); i++) {
        // 获取当前三角形的三个顶点索引
        int p1 = triangles[i].vertices[0];
        int p2 = triangles[i].vertices[1];
        int p3 = triangles[i].vertices[2];
        mesh->add_polygon("gme", i, 3);
        // 获取刚刚添加的多边形poly0的指针，并设置其顶点信息
        indexed_polygon* poly0 = mesh->get_polygon("gme", i);
        poly0->set_vertex("gme", 0, &mesh->get_vertex("gme", p1));
        poly0->set_vertex("gme", 1, &mesh->get_vertex("gme", p2));
        poly0->set_vertex("gme", 2, &mesh->get_vertex("gme", p3));
    }
    // 将索引化网格mesh附加到面f上
    attach_indexed_mesh_to_face(f, mesh);
}

/**
 * @brief 处理不同的边界和环，对环面进行网格化处理，并将生成的网格附加到面对象上
 *
 * @return 处理结果
 */
outcome torus_faceter::facet() {
    init();
    decideUVlen();

    // 获取环的信息
    bool has_periphery_loop = getLoops();
    // 有v分离环，调用VSeperationProcess()函数处理
    if(v_loops.size() > 0) {
        has_periphery_loop = VSeperationProcess();
    }

    // 有外环，调用peripheryProcess()函数处理
    if(has_periphery_loop) {
        // if(periphery.size() >= 2)
        peripheryProcess();
    } else {  // 如果没有周边环，就创建u、v方向上的边界点和边，这些边界点和边将作为网格的边界
        int index1, index2;
        // 创建v方向上的边界点和边
        for(int i = 0; i <= vlen; i++) {
            // 创建边界点
            CDT::V2d<double> tempL = {u_le, v_le + i * dv};
            CDT::V2d<double> tempH = {u_ri, v_le + i * dv};
            points.push_back(tempL);
            points.push_back(tempH);
            index1 = points.size() - 1;
            index2 = points.size() - 2;
            // 创建边
            edges.push_back(CDT::Edge(index2, index2 + 2));
            edges.push_back(CDT::Edge(index1, index1 + 2));
            // if(i < vlen) {
            //     printf("%lf %lf %lf %lf\n", u_le, v_le + i * dv, u_le, v_le + (i + 1) * dv);
            //     printf("%lf %lf %lf %lf\n", u_ri, v_le + i * dv, u_ri, v_le + (i + 1) * dv);
            // }
        }
        edges.pop_back();
        edges.pop_back();
        // 创建水平方向上的边界点和边
        for(int i = 0; i <= ulen; i++) {
            // 创建边界点
            CDT::V2d<double> tempL = {u_le + i * du, v_le};
            CDT::V2d<double> tempH = {u_le + i * du, v_ri};
            points.push_back(tempL);
            points.push_back(tempH);
            index1 = points.size() - 1;
            index2 = points.size() - 2;
            // 创建边
            edges.push_back(CDT::Edge(index2, index2 + 2));
            edges.push_back(CDT::Edge(index1, index1 + 2));
            // if(i < ulen) {
            //     printf("%lf %lf %lf %lf\n", u_le + i * du, v_le, u_le + (i + 1) * du, v_le);
            //     printf("%lf %lf %lf %lf\n", u_le + i * du, v_le, u_le + (i + 1) * du, v_le);
            // }
        }
        edges.pop_back();
        edges.pop_back();
    }

    // 有u分离环，调用USeperationProcess()函数处理
    if(u_loops.size() > 0) {
        USeperationProcess();
    }

    // 有孔环，调用holeProcess()函数处理
    if(hole_loops.size() > 0) {
        holeProcess();
    }

    // 调用attachMesh()函数生成网格并附加到面对象上
    attachMesh();
    return outcome();
}

/**
 * @brief 对环面进行离散
 *
 * @param f 待处理的环面
 * @param nt 法向容差
 * @param st 距离容差
 * @return 返回网格化处理结果
 */
outcome gme_facet_face_torus(FACE* f, double nt, double st) {
    torus_faceter faceter(f, nt, st);
    return faceter.facet();
}