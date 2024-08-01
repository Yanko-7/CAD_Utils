#include "acis/gme/faceter/gme_facet_face_spline.hxx"

#include <CDT/include/CDTUtils.h>

#include <acis/include/geometry.hxx>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <kdtree/include/KDTree.hpp>
#include <set>

#include "acis/gme/faceter/gme_facet_face_torus.hxx"
#include "acis/include/acistol.hxx"
#include "acis/include/add_pcu.hxx"
#include "acis/include/cucuint.hxx"
#include "acis/include/curdef.hxx"
#include "acis/include/getbox.hxx"
#include "acis/include/getowner.hxx"
#include "acis/include/intrapi.hxx"
#include "acis/include/pcurve.hxx"
#include "acis/include/spldef.hxx"
#include "acis/include/spline.hxx"
#include "acis/include/strdef.hxx"
#include "acis/include/transfrm.hxx"

// #define DEBUG
#ifdef DEBUG
#    define TIMERSTART(tag) auto tag##_start = std::chrono::high_resolution_clock::now();
#else
#    define TIMERSTART(tag)
#endif  // DEBUG

#ifdef DEBUG
#    define TIMEREND(tag)                                           \
        auto tag##_end = std::chrono::high_resolution_clock::now(); \
        printf("|| %-40s|| time costs %f ms\n", #tag, std::chrono::duration<double, std::milli>(tag##_end - tag##_start).count());
#else
#    define TIMEREND(tag)
#endif  // DEBUG

// #define TIMERSTART(tag) auto tag##_start = std::chrono::high_resolution_clock::now();
// #def ine TIMEREND(tag)                                           \
//    auto tag##_end = std::chrono::high_resolution_clock::now(); \
//    printf("|| %-40s|| time costs %f ms\n", #tag, std::chrono::duration<double, std::milli>(tag##_end - tag##_start).count());

std::map<FACE*, my_mesh> spline_faceter::spline_map_mesh = {};
singular_type spline_faceter::is_singular(point a) {
    if(singular_ule && equal(a.x, u_le))
        return ule;
    else if(singular_uri && equal(a.x, u_ri))
        return uri;
    else if(singular_vle && equal(a.y, v_le))
        return vle;
    else if(singular_vri && equal(a.y, v_ri))
        return vri;
    return no_singular;
}

point_boundary_type spline_faceter::is_on_boundary(point a) {
    bool bule = equal(a.x, u_le);
    bool buri = equal(a.x, u_ri);
    bool bvle = equal(a.y, v_le);
    bool bvri = equal(a.y, v_ri);
    if(bule) {
        if(bvle)
            return b_lele;
        else if(bvri)
            return b_leri;
        return b_ule;
    }
    if(buri) {
        if(bvle)
            return b_rile;
        else if(bvri)
            return b_riri;
        return b_uri;
    }
    if(bvle) return b_vle;
    if(bvri) return b_vri;
    return not_boundary;
}
SPApar_dir spline_faceter::interior_vector(point a, point b) {
    assert(a != b);
    SPAposition a3 = spl->eval_position(SPApar_pos(a.x, a.y));
    SPAposition b3 = spl->eval_position(SPApar_pos(b.x, b.y));
    SPApar_pos a2(a.x, a.y);
    SPAvector ab3 = b3 - a3;
    SPAunit_vector an = spl->eval_normal(a2);
    if(f->sense()) an = -an;
    SPAvector dir3 = an * ab3;
    SPAunit_vector udir3(dir3.x(), dir3.y(), dir3.z());
    return spl->param_dir(udir3, a2);
}
void spline_faceter::singular_preprocess(point_vector& loop) {
    point curr, next;
    singular_type at, bt;
    for(int i = 1; i < loop.size(); i++) {
        curr = loop[i];
        next = loop[(i + 1) % loop.size()];
        at = is_singular(curr);
        bt = is_singular(next);
        if(at != no_singular && bt != no_singular) {
            // if(at != bt) printf("error singularity\n");
            switch(at) {
                case ule:

                    if(next.y < curr.y && !equal(next.y, curr.y)) {
                        loop.insert(loop.begin() + i + 1, point(u_le, v_ri));
                        loop.insert(loop.begin() + i + 2, point(u_le, v_le));
                        i += 3;
                    }
                    break;
                case uri:
                    if(next.y > curr.y && !equal(next.y, curr.y)) {
                        loop.insert(loop.begin() + i + 1, point(u_ri, v_le));
                        loop.insert(loop.begin() + i + 2, point(u_ri, v_ri));
                        i += 3;
                    }
                    break;
                case vle:
                    if(next.x > curr.x && !equal(next.x, curr.x)) {
                        loop.insert(loop.begin() + i + 1, point(u_le, v_le));
                        loop.insert(loop.begin() + i + 2, point(u_ri, v_le));
                        i += 3;
                    }
                    break;
                case vri:
                    if(next.x < curr.x && !equal(next.x, curr.x)) {
                        loop.insert(loop.begin() + i + 1, point(u_ri, v_ri));
                        loop.insert(loop.begin() + i + 2, point(u_le, v_ri));
                        i += 3;
                    }
                    break;
            }
        }
    }
}

QuadTree::QuadTree(FACE* f, double nt, double st, std::vector<SPApar_vec> vertex_set = std::vector<SPApar_vec>(0)): face_(f), distance_tolerance_(st), angle_tolerance_(nt), mRoot_(std::make_unique<Node>()) {
    SPAinterval u_interval = f->geometry()->equation().param_range_u();
    SPAinterval v_interval = f->geometry()->equation().param_range_v();
    mBox_ = Box(u_interval.start_pt(), v_interval.start_pt(), u_interval.length(), v_interval.length());
    mRoot_->vertex_set = std::move(vertex_set);
    BuildTreeByBin(mBox_, 1);
}

void QuadTree::BuildTree(const Box& box, Node* root, int dep) {
    // 如果深度大于1，且按照其中一条对角线切成两个三角形拟合可以满足条件，则不再继续分割
    if(dep > 1 && CheckFitAndJudgeDiagonal(face_, box, angle_tolerance_, distance_tolerance_) != -1) {
        return;
    }
    // 矩形中点坐标
    double&& mid_u = (box.u_ + box.u_ + box.width_) / 2;
    double&& mid_v = (box.v_ + box.v_ + box.height_) / 2;
    split(root, box);  // 将当前节点存的点集划分到子节点中
    //// 0 1
    //// 2 3 放入
    for(int i = 0; i < 4; i++) {
        BuildTree(ComputeBox(box, i), root->children[i].get(), dep + 1);
    }
    // 生成五个采样点
    result_points_.emplace_back(box.u_, mid_v);
    result_points_.emplace_back(mid_u, mid_v);
    result_points_.emplace_back(mid_u, box.v_ + box.height_);
    result_points_.emplace_back(box.u_ + box.width_, box.v_ + box.height_ / 2);
    result_points_.emplace_back(box.u_ + box.width_ / 2, box.v_);
}

void QuadTree::split(Node* node, const Box& box) {
    assert(node != nullptr);
    assert(isLeaf(node) && "Only leaves can be split");

    for(auto& child: node->children) child = std::make_unique<Node>();

    auto new_vertex_set = std::vector<SPApar_vec>();
    for(const auto& vertex: node->vertex_set) {
        auto i = GetQuadrant(box, vertex);
        if(i != -1) {
            node->children[static_cast<std::size_t>(i)]->vertex_set.push_back(vertex);
        } else {
            new_vertex_set.push_back(vertex);
        }
    }
    node->vertex_set = std::move(new_vertex_set);
}

inline int CheckFitAndJudgeDiagonal(FACE* face, const Box& box, const double angle_tolerance, const double distance_tolerance) {
    SPApar_pos tmp_p0[3];
    SPApar_pos tmp_p1[3];
    SPApar_pos tmp_p2[3];
    SPApar_pos tmp_p3[3];
    // tmp_p0为左下三角形
    tmp_p0[0] = {box.u_, box.v_ + box.height_};
    tmp_p0[1] = {box.u_ + box.width_, box.v_};
    tmp_p0[2] = {box.u_, box.v_};
    // tmp_p1为右上三角形
    tmp_p1[0] = {box.u_, box.v_ + box.height_};
    tmp_p1[1] = {box.u_ + box.width_, box.v_ + box.height_};
    tmp_p1[2] = {box.u_ + box.width_, box.v_};
    // tmp_p2为左上三角形
    tmp_p2[0] = {box.u_, box.v_ + box.height_};
    tmp_p2[1] = {box.u_ + box.width_, box.v_ + box.height_};
    tmp_p2[2] = {box.u_, box.v_};
    // tmp_p3为右下三角形
    tmp_p3[0] = {box.u_, box.v_};
    tmp_p3[1] = {box.u_ + box.width_, box.v_};
    tmp_p3[2] = {box.u_ + box.width_, box.v_ + box.height_};

    const SPAposition& pos_botLeft = face->geometry()->equation().eval_position({box.u_, box.v_});
    const SPAposition& pos_topRight = face->geometry()->equation().eval_position({box.u_ + box.width_, box.v_ + box.height_});
    const SPAposition& pos_botRight = face->geometry()->equation().eval_position({box.u_ + box.width_, box.v_});
    const SPAposition& pos_topLeft = face->geometry()->equation().eval_position({box.u_, box.v_ + box.height_});

    const SPAposition& pos_center = face->geometry()->equation().eval_position({(box.u_ + box.width_ + box.u_) / 2, (box.v_ + box.height_ + box.v_) / 2});

    // 主对角线上的中点
    const SPAposition& main = pos_botRight + (pos_topLeft - pos_botRight) / 2;

    // 副对角线上的中点
    const SPAposition& second = pos_botLeft + (pos_topRight - pos_botLeft) / 2;
    // 0为主对角线，1为副对角线
    if((pos_center - main).len() <= (pos_center - second).len()) {
        if(checkTriangle_(face, angle_tolerance, distance_tolerance, tmp_p0) && checkTriangle_(face, angle_tolerance, distance_tolerance, tmp_p1)) {
            return 0;
        }
    } else {
        if(checkTriangle_(face, angle_tolerance, distance_tolerance, tmp_p2) && checkTriangle_(face, angle_tolerance, distance_tolerance, tmp_p3)) {
            return 1;
        }
    }
    return -1;
}
const Box QuadTree::ComputeBox(const Box& box, int i) const {
    auto origin = box.getBotLeft();
    auto childSize = box.getSize() / 2;
    switch(i) {
        // North West
        case 0:
            return Box(SPApar_vec(origin.du, origin.dv + childSize.dv), childSize);
        // Norst East
        case 1:
            return Box(origin + childSize, childSize);
        // South West
        case 2:
            return Box(origin, childSize);
        // South East
        case 3:
            return Box(SPApar_vec(origin.du + childSize.du, origin.dv), childSize);
        default:
            assert(false && "Invalid child index");
            return Box();
    }
}
int QuadTree::GetQuadrant(const Box& nodeBox, const SPApar_vec& vertex) const {
    if(vertex.du < nodeBox.u_ || vertex.du > nodeBox.u_ + nodeBox.width_ || vertex.dv < nodeBox.v_ || vertex.dv > nodeBox.v_ + nodeBox.height_) {
        return -1;
    }
    auto center = nodeBox.getCenter();
    // West
    if(vertex.du < center.du) {
        // North West
        if(vertex.dv >= center.dv) return 0;
        // South West
        else
            return 2;
    }
    // East
    else {
        // North East
        if(vertex.dv >= center.dv) return 1;
        // South East
        else
            return 3;
    }
}
bool QuadTree::isLeaf(const Node* node) const {
    return !static_cast<bool>(node->children[0]);
}

void QuadTree::BuildTreeByBin(const Box& box, int dep) {
    if(dep > 1 && CheckFitAndJudgeDiagonal(face_, box, angle_tolerance_, distance_tolerance_) != -1) {
        // leaf_rect_.push_back(box);
        result_points_.emplace_back(box.u_, box.v_);
        result_points_.emplace_back(box.u_ + box.width_, box.v_);
        result_edges_.emplace_back(result_points_.size() - 1, result_points_.size() - 2);
        result_points_.emplace_back(box.u_ + box.width_, box.v_ + box.height_);
        result_edges_.emplace_back(result_points_.size() - 1, result_points_.size() - 2);
        result_points_.emplace_back(box.u_, box.v_ + box.height_);
        result_edges_.emplace_back(result_points_.size() - 1, result_points_.size() - 2);
        result_edges_.emplace_back(result_points_.size() - 1, result_points_.size() - 4);
        return;
    }
    // 矩形中点坐标
    auto& spl = face_->geometry()->equation();
    auto lineu = spl.u_param_line(box.v_ + box.height_ / 2);
    double ul = lineu->length(box.u_, box.u_ + box.width_);
    auto linev = spl.v_param_line(box.u_ + box.width_ / 2);
    double vl = linev->length(box.v_, box.v_ + box.height_);
    double gu = (spl.eval_position({box.u_, box.v_ + box.height_ / 2}) - spl.eval_position({box.u_ + box.width_, box.v_ + box.height_ / 2})).len();
    double gv = (spl.eval_position({box.u_ + box.width_ / 2, box.v_}) - spl.eval_position({box.u_ + box.width_ / 2, box.v_ + box.height_})).len();
    double ur = 0;
    double vr = 0;
    if(gu != 0) {
        ur = (ul - gu) / gu;
    }
    if(gv != 0) {
        vr = (vl - gv) / gv;
    }
    if(gv < SPAresabs && gu < SPAresabs) {
        BuildTreeByBin(Box(box.u_, box.v_, box.width_ / 2, box.height_ / 2), dep + 1);
        BuildTreeByBin(Box(box.u_ + box.width_ / 2, box.v_, box.width_ / 2, box.height_ / 2), dep + 1);
        BuildTreeByBin(Box(box.u_ + box.width_ / 2, box.v_ + box.height_ / 2, box.width_ / 2, box.height_ / 2), dep + 1);
        BuildTreeByBin(Box(box.u_, box.v_ + box.height_ / 2, box.width_ / 2, box.height_ / 2), dep + 1);

        double&& mid_u = (box.u_ + box.u_ + box.width_) / 2;
        double&& mid_v = (box.v_ + box.v_ + box.height_) / 2;
        // result_points_.emplace_back(box.u_, mid_v);
        // result_points_.emplace_back(mid_u, mid_v);
        // result_points_.emplace_back(mid_u, box.v_ + box.height_);
        // result_points_.emplace_back(box.u_ + box.width_, box.v_ + box.height_ / 2);
        // result_points_.emplace_back(box.u_ + box.width_ / 2, box.v_);
        return;
    }
    // double ur = get_box_angle_mean(Box(box.u_, box.v_, box.width_ / 2, box.height_)) + get_box_angle_mean(Box(box.u_ + box.width_ / 2, box.v_, box.width_ / 2, box.height_));
    // double vr = get_box_angle_mean(Box(box.u_, box.v_, box.width_, box.height_ / 2)) + get_box_angle_mean(Box(box.u_, box.v_ + box.height_ / 2, box.width_, box.height_ / 2));
    if(gu < SPAresabs || (ur > vr && gv > SPAresabs)) {
        BuildTreeByBin(Box(box.u_, box.v_, box.width_ / 2, box.height_), dep + 1);
        BuildTreeByBin(Box(box.u_ + box.width_ / 2, box.v_, box.width_ / 2, box.height_), dep + 1);
        // result_points_.emplace_back(box.u_ + box.width_ / 2, box.v_);
        // result_points_.emplace_back(box.u_ + box.width_ / 2, box.v_ + box.height_);
        // result_edges_.emplace_back(result_points_.size() - 1, result_points_.size() - 2);
    } else {
        BuildTreeByBin(Box(box.u_, box.v_, box.width_, box.height_ / 2), dep + 1);
        BuildTreeByBin(Box(box.u_, box.v_ + box.height_ / 2, box.width_, box.height_ / 2), dep + 1);
        // result_points_.emplace_back(box.u_, box.v_ + box.height_ / 2);
        // result_points_.emplace_back(box.u_ + box.width_, box.v_ + box.height_ / 2);
        // result_edges_.emplace_back(result_points_.size() - 1, result_points_.size() - 2);
    }
    return;
}
double QuadTree::get_box_angle_mean(Box box) {
    const int num = 5;
    const double du = box.width_ / (num + 1);
    const double dv = box.height_ / (num + 1);
    auto& spl = face_->geometry()->equation();
    static std::vector<std::vector<SPAunit_vector>> vec(num + 2, std::vector<SPAunit_vector>(num + 2));
    for(int i = 0; i <= num + 1; i++) {
        for(int j = 0; j <= num + 1; j++) {
            vec[i][j] = spl.eval_normal({box.u_ + du * (i), box.v_ + dv * (j)});
        }
    }
    double sum = 0;
    for(int i = 1; i <= num + 1; i++) {
        for(int j = 0; j <= num + 1; j++) {
            auto& one = vec[i - 1][j];
            auto& two = vec[i][j];
            auto theta = std::min(one % two, 1.0);
            sum += acos(theta);
        }
    }
    for(int i = 0; i <= num + 1; i++) {
        for(int j = 1; j <= num + 1; j++) {
            auto& one = vec[i][j - 1];
            auto& two = vec[i][j];
            auto theta = std::min(one % two, 1.0);
            sum += acos(theta);
        }
    }
    return sum;
}
void spline_faceter::init() {
    spl = (spline*)&f->geometry()->equation();

    SPApar_box face_par_box = f->geometry()->equation().param_range();
    SPAinterval u_interval = f->geometry()->equation().param_range_u();
    SPAinterval v_interval = f->geometry()->equation().param_range_v();
    u_le = u_interval.start_pt();
    u_ri = u_interval.end_pt();
    v_le = v_interval.start_pt();
    v_ri = v_interval.end_pt();

    u_tolerance = 1e15;
    v_tolerance = 1e15;

    if(f->geometry()->equation().closed_u()) {
        u_tolerance = spl->param_range_u().length() / 2;  // 判断跨边界所用值
        if(spl->singular_u(v_le)) {
            singular_vle = true;
        } else {  // B_7  B14的样例需要通过下列判断才能确定奇点情况
            curve* vle = f->geometry()->equation().v_param_line(v_le);
            SPAinterval interval = vle->param_range();
            if(is_equal(vle->length(interval.start_pt(), interval.end_pt()), 0.0)) {
                singular_vle = true;
            }
        }
        if(spl->singular_u(v_ri)) {
            singular_vri = true;
        } else {
            curve* vri = f->geometry()->equation().v_param_line(v_ri);
            SPAinterval interval = vri->param_range();
            if(is_equal(vri->length(interval.start_pt(), interval.end_pt()), 0.0)) {
                singular_vri = true;
            }
        }
    }
    if(f->geometry()->equation().closed_v()) {
        v_tolerance = spl->param_range_v().length() / 2;

        if(spl->singular_u(u_le)) {
            singular_ule = true;
        } else {  // B_7  B14的样例需要通过下列判断才能确定奇点情况
            curve* ule = f->geometry()->equation().v_param_line(u_le);
            SPAinterval interval = ule->param_range();
            if(is_equal(ule->length(interval.start_pt(), interval.end_pt()), 0.0)) {
                singular_ule = true;
            }
        }
        if(spl->singular_u(u_ri)) {
            singular_uri = true;
        } else {
            curve* uri = f->geometry()->equation().v_param_line(u_ri);
            SPAinterval interval = uri->param_range();
            if(is_equal(uri->length(interval.start_pt(), interval.end_pt()), 0.0)) {
                singular_uri = true;
            }
        }
    }
    le_low.x = u_le;  // 左下角点
    le_low.y = v_le;
    ri_high.x = u_ri;
    ri_high.y = v_ri;
}

void spline_faceter::decideUVlen() {
    ulen = 5;
    vlen = 5;
    du = (u_ri - u_le) / ulen;
    dv = (v_ri - v_le) / vlen;
}

void spline_faceter::peripheryProcess() {
    point curr, next;
    bool closed_u = f->geometry()->equation().closed_u();
    bool closed_v = f->geometry()->equation().closed_v();
    // std::vector<point_vector> temp = {periphery};
    // writeLoopsToFile(temp, "spline2.txt");
    // 奇点跳变需要处理
    FaceterFaceInfo info(u_ri, u_le, v_ri, v_le, u_tolerance, v_tolerance, closed_u, closed_v);
    info.du = du;
    info.dv = dv;
    // FPLoopPreprocessor(periphery, info);
    if((singular_ule || singular_uri) && f->geometry()->equation().closed_v() || (singular_vle || singular_vri) && f->geometry()->equation().closed_u()) {
        preprocessor(periphery, closed_u, closed_v, le_low, ri_high);
        // 取相邻的两个点
        for(int i = 0; i < periphery.size(); i++) {
            curr = periphery[i];
            next = periphery[(i + 1) % periphery.size()];
            // 奇点处理
            // 如果两个点都位于u_le奇点边界上,则补全该两点之间的点，以下同理
            if(singular_ule && is_equal(curr.x, u_le) && is_equal(next.x, u_le)) {
                point_vector res = SingularityCompletion(curr.y, next.y, SingularityType::kUle, info);
                periphery.insert(periphery.begin() + i + 1, res.begin(), res.end());
            } else if(singular_uri && is_equal(curr.x, u_ri) && is_equal(next.x, u_ri)) {
                point_vector res = SingularityCompletion(curr.y, next.y, SingularityType::kUri, info);
                periphery.insert(periphery.begin() + i + 1, res.begin(), res.end());
            } else if(singular_vle && is_equal(curr.y, v_le) && is_equal(curr.y, v_le)) {
                point_vector res = SingularityCompletion(curr.x, next.x, SingularityType::kVle, info);
                periphery.insert(periphery.begin() + i + 1, res.begin(), res.end());
            } else if(singular_vri && is_equal(curr.y, v_ri) && is_equal(curr.y, v_ri)) {
                point_vector res = SingularityCompletion(curr.x, next.x, SingularityType::kVri, info);
                periphery.insert(periphery.begin() + i + 1, res.begin(), res.end());
            }
        }
    }
    // 打个补丁，保证可靠性
    if(closed_u || closed_v) {
        preprocessor(periphery, closed_u, closed_v, le_low, ri_high);
    }
    // 分割跨界loop成不跨界的小loop
    std::vector<point_vector> small_loops;
    boundary_loop_split(periphery, le_low, ri_high, u_tolerance, v_tolerance, small_loops);
    for(int j = 0; j < small_loops.size(); j++) {  // 遍历每个loop集
        point_vector loop = std::move(small_loops[j]);
        size_t start = points.size();           // 起始位置从最后开始
        bool boundary_point_end = false;        // 边界点结尾标志，此处未使用
        for(int p = 0; p < loop.size(); p++) {  // 遍历loop中的每个点
            points.push_back(loop[p]);          // 将当前点添加到点集中
            // 设置相邻点
            curr = loop[p];
            next = loop[(p + 1) % loop.size()];
            if(next == curr) {  // 如果下一个点与当前点相同，则丢弃当前点（重复点），这样可以处理连续相同点
                points.pop_back();
                continue;
            }
            // 如果两个点不同，添加边
            edges.push_back(CDT::Edge(points.size() - 1, points.size()));
            // 如果两点位于边界，需要添加在面内的点，这一段与奇点处理比较像
            if((fabs(curr.x - u_ri) < EPSILON && fabs(next.x - u_ri) < EPSILON) || (fabs(curr.x - u_le) < EPSILON && fabs(next.x - u_le) < EPSILON)) {
                if(fabs(curr.y - next.y) < dv) continue;  // 如果y方向的距离小于阈值，则忽略
                if(curr.y < next.y) {
                    for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low < next.y - EPSILON; low += dv) {
                        if(low < curr.y) continue;
                        CDT::V2d<double> temp = {curr.x, low};
                        points.push_back(temp);
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                } else {
                    for(double low = get_cloest_boundary_point(v_le, dv, curr.y); low > next.y + EPSILON; low -= dv) {
                        if(low > curr.y) continue;
                        CDT::V2d<double> temp = {curr.x, low};
                        points.push_back(temp);
                        edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                    }
                }
            } else if((fabs(curr.y - v_ri) < EPSILON && fabs(next.y - v_ri) < EPSILON) || (fabs(curr.y - v_le) < EPSILON && fabs(next.y - v_le) < EPSILON)) {
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
        edges.pop_back();                                      // 移除最后一条边
        edges.push_back(CDT::Edge(points.size() - 1, start));  // 添加一条从最后一个点到起始点的边，闭合loop
    }
}

bool spline_faceter::USeperationProcess() {
    point curr, next;                               // 定义当前点和下一个点
    SPApar_pos param;                               // 定义参数位置
    point le_high(u_le, v_ri), ri_low(u_ri, v_le);  // 定义左上和右下的点

    std::vector<unsigned> u_points;  // 记录跨边界的点的索引

    for(unsigned i = 0; i < u_loops.size(); i++) {
        std::vector<CDT::V2d<double>> u_loop = u_loops[i];
        preprocessor(u_loop, f->geometry()->equation().closed_u(), f->geometry()->equation().closed_v(), le_low, ri_high);
        // printLoop(u_loop);
        // preprocessor(u_loop, closed_u, closed_v, le_low, ri_high);
        // printf("\n");
        //
        int istart = 0;
        curr = u_loop[istart++];

        // istart = 4
        unsigned index = 0;                                // 点的索引
        for(int i = 0; i < u_loop.size(); i++) {           // 遍历loop中的每个点
            if(i == 0) istart = points.size();             // 如果是loop的第一个点，记录起始索引
            curr = u_loop[i];                              // 获取当前点
            points.push_back(curr);                        // 将当前点添加到点集中
            index = points.size();                         // 更新点的索引
            edges.push_back(CDT::Edge(index - 1, index));  // 添加一条边
            if(u_loop.size() == 2 && i == 1) {             // 如果loop只有两个点，跳过以避免重复计算
                continue;
            }
            next = u_loop[(i + 1) % u_loop.size()];                                    // 获取下一个点
            if(fabs(next.y - curr.y) > v_tolerance) {                                  // 如果y方向的差异超过容差
                CDT::V2d<double> temp = {0, 0};                                        // 临时点
                u_points.push_back(index);                                             // 记录点的索引
                if(next.y > curr.y) {                                                  // 如果下一个点在当前点的上方
                    next.y -= (v_ri - v_le);                                           // 调整y坐标
                    bool has_intersect = intersect(curr, next, le_low, ri_low, temp);  // 计算交点
                    if(!has_intersect) {                                               // 如果没有交点
                        temp.x = curr.x;                                               // 使用当前点的x坐标
                        temp.y = v_le;                                                 // 使用下边界的y坐标
                    }
                    // printf("%d\n", has_intersect);
                    points.push_back(temp);                                            // 添加交点
                    temp.y = v_ri;                                                     // 调整y坐标到上边界
                    points.push_back(temp);                                            // 添加交点
                } else {                                                               // 如果下一个点在当前点的下方
                    next.y += (v_ri - v_le);                                           // 调整y坐标
                    bool has_intersect = intersect(curr, next, le_low, ri_low, temp);  // 计算交点
                    if(!has_intersect) {                                               // 如果没有交点
                        temp.x = curr.x;                                               // 使用当前点的x坐标
                        temp.y = v_ri;                                                 // 使用上边界的y坐标
                    }
                    points.push_back(temp);  // 添加交点
                    temp.y = v_le;           // 调整y坐标到下边界
                    points.push_back(temp);  // 添加交点
                }
                edges.push_back(CDT::Edge(index + 1, index + 2));  // 添加边
            }
        }
        edges.pop_back();                               // 移除最后一条边
        index = points.size();                          // 更新点的索引
        edges.push_back(CDT::Edge(index - 1, istart));  // 添加一条边，闭合循环
    }
    // 对跨边界的点进行排序，确保顺序正确
    for(int i = 0; i < u_points.size(); i++) {
        unsigned ind = u_points[i];  // 获取点的索引
        for(int j = i + 1; j < u_points.size(); j++) {
            unsigned ind2 = u_points[j];          // 获取另一个点的索引
            if(points[ind].x > points[ind2].x) {  // 如果第一个点在第二个点的右侧
                unsigned tt = u_points[i];        // 交换两个点的索引
                u_points[i] = u_points[j];
                u_points[j] = tt;
                break;
            }
        }
    }

    for(int i = 0; i < u_points.size(); i++) {  // 遍历所有跨边界的点
        curr = points[u_points[i]];             // 获取当前点
        next = points[u_points[i] + 1];         // 获取下一个点
        if(next.y > curr.y) {                   // 如果下一个点在当前点的上方
            // 找到下一个指向上的点或者边界
            // 如果找不到
            if(i == u_points.size() - 1) {             // 如果是最后一个点
                CDT::V2d<double> temp = {u_ri, v_le};  // 创建一个临时点
                points.push_back(temp);                // 添加临时点
                temp.y = v_ri;                         // 调整y坐标
                points.push_back(temp);                // 添加临时点

                edges.push_back(CDT::Edge(u_points[i], points.size() - 2));  // 添加边
                edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));
                edges.push_back(CDT::Edge(points.size() - 1, u_points[i] + 1));
            } else {
                edges.push_back(CDT::Edge(u_points[i], u_points[i + 1] + 1));  // 添加边
                edges.push_back(CDT::Edge(u_points[i] + 1, u_points[i + 1]));
            }
        }
        if(i == 0 && next.y < curr.y) {            // 如果是第一个点且下一个点在当前点的下方
            CDT::V2d<double> temp = {u_le, v_ri};  // 创建一个临时点
            points.push_back(temp);                // 添加临时点
            temp.y = v_le;                         // 调整y坐标
            points.push_back(temp);                // 添加临时点

            edges.push_back(CDT::Edge(u_points[i], points.size() - 2));  // 添加边
            edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));
            edges.push_back(CDT::Edge(points.size() - 1, u_points[i] + 1));
        }
    }
    return false;  // 返回false
}

bool spline_faceter::VSeperationProcess() {
    point curr, next;                               // 当前点和下一个点
    SPApar_pos param;                               // 参数位置
    point le_high(u_le, v_ri), ri_low(u_ri, v_le);  // 左高和右低点
    std::vector<unsigned> v_points;                 // 记录跨边界的点的索引

    // 遍历所有的垂直循环
    for(unsigned i = 0; i < v_loops.size(); i++) {
        std::vector<CDT::V2d<double>> v_loop = v_loops[i];
        int istart = 0;
        curr = v_loop[istart++];
        bool PI_in_right = false;
        unsigned index = 0;
        // 试图取到第一个u值不为-PI的点，来判断PI该取正值还是负值
        if(spl->closed_v()) {
            while(fabs(curr.y - v_le) <= SPAresabs && istart < v_loop.size()) {
                curr = v_loop[istart++];
            }
            if(istart == v_loop.size()) {
                curr.y += 2 * SPAresabs;
                param.u = curr.x;
                param.v = curr.y;
            }
            if(curr.y > 0) {
                PI_in_right = true;
            }
        }

        // 遍历循环中的每个点
        for(int i = 0; i < v_loop.size(); i++) {
            if(i == 0) istart = points.size();  // 设置起始索引
            curr = v_loop[i];                   // 当前点

            // 如果PI在右侧并且当前点的y值接近v_le
            if(PI_in_right && fabs(curr.y - v_le) <= SPAresabs) {
                curr.y = v_ri;  // 更改y值到v_ri
            }

            points.push_back(curr);                        // 添加当前点到points
            index = points.size();                         // 更新索引
            edges.push_back(CDT::Edge(index - 1, index));  // 添加边
            next = v_loop[(i + 1) % v_loop.size()];        // 下一个点

            // 如果下一个点和当前点的x值差距大于容差
            if(fabs(next.x - curr.x) > u_tolerance) {
                CDT::V2d<double> temp = {0, 0};  // 临时点
                v_points.push_back(index);       // 记录索引
                // 左穿右
                if(next.x > curr.x) {
                    next.x -= u_ri - u_le;                                              // 调整x值
                    bool has_intersect = intersect(next, curr, le_low, le_high, temp);  // 检查交点
                    if(!has_intersect) {
                        temp.x = u_le;    // 设置x值为u_le
                        temp.y = curr.y;  // 设置y值
                    }
                    points.push_back(temp);                                             // 添加临时点
                    temp.x = u_ri;                                                      // 设置x值为u_ri
                    points.push_back(temp);                                             // 添加临时点
                } else {                                                                // 右穿左
                    next.x += u_ri - u_le;                                              // 调整x值
                    bool has_intersect = intersect(next, curr, ri_low, ri_high, temp);  // 检查交点
                    if(!has_intersect) {
                        temp.x = u_ri;    // 设置x值为u_ri
                        temp.y = curr.y;  // 设置y值
                    }
                    points.push_back(temp);  // 添加临时点
                    temp.x = u_le;           // 设置x值为u_le
                    points.push_back(temp);  // 添加临时点
                }
                edges.push_back(CDT::Edge(index + 1, index + 2));  // 添加边
            }
        }
        edges.pop_back();                               // 移除最后一条边
        index = points.size();                          // 更新索引
        edges.push_back(CDT::Edge(index - 1, istart));  // 添加边
    }

    bool isSame = false;  // 标记是否相同
    // 对v_points进行排序
    for(int i = 0; i < v_points.size(); i++) {
        unsigned ind = v_points[i];  // 当前索引
        for(int j = i + 1; j < v_points.size(); j++) {
            unsigned ind2 = v_points[j];  // 下一个索引
            // 如果y值接近相等
            if(fabs(points[ind].y - points[ind2].y) < SPAresabs) {
                isSame = true;  // 标记为相同
                break;
            } else if(points[ind].y > points[ind2].y) {  // 排序
                unsigned tt = v_points[i];
                v_points[i] = v_points[j];
                v_points[j] = tt;
                break;
            }
        }
    }

    // 如果有相同的点，移除最后两个点
    if(isSame) {
        v_points.pop_back();
        v_points.pop_back();
    }
    // 处理边界连接
    for(int i = 0; i < v_points.size(); i++) {
        curr = points[v_points[i]];      // 当前点
        next = points[v_points[i] + 1];  // 下一个点
        // 如果指向右
        if(next.x < curr.x) {
            if(i == v_points.size() - 1) {  // 如果是最后一个点
                points.push_back(ri_high);  // 添加ri_high
                points.push_back(le_high);  // 添加le_high

                edges.push_back(CDT::Edge(v_points[i], points.size() - 2));        // 添加边
                edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));  // 添加边
                edges.push_back(CDT::Edge(points.size() - 1, v_points[i] + 1));
            } else {
                edges.push_back(CDT::Edge(v_points[i], v_points[i + 1] + 1));  // 添加边
                edges.push_back(CDT::Edge(v_points[i] + 1, v_points[i + 1]));  // 添加边
            }
        }
        // 如果是第一个点且指向左
        if(i == 0 && next.x > curr.x) {
            points.push_back(le_low);  // 添加le_low
            points.push_back(ri_low);  // 添加ri_low

            edges.push_back(CDT::Edge(v_points[i], points.size() - 2));        // 添加边
            edges.push_back(CDT::Edge(points.size() - 2, points.size() - 1));  // 添加边
            edges.push_back(CDT::Edge(points.size() - 1, v_points[i] + 1));    // 添加边
        }
    }
    return false;  // 返回false
}
void spline_faceter::holeProcess() {
    // 检查是否有孔洞Loop需要处理
    if(hole_loops.size()) {
        // 遍历所有孔洞Loop进行实际处理
        for(int i = 0; i < hole_loops.size(); i++) {
            // 存储分割后的小Loop的容器
            std::vector<point_vector> small_loops;
            // 根据容差将当前孔洞Loop分割成小Loop
            boundary_loop_split(hole_loops[i], le_low, ri_high, u_tolerance, v_tolerance, small_loops);
            // 遍历所有小Loop
            for(int j = 0; j < small_loops.size(); j++) {
                // 当前的小Loop
                point_vector loop = small_loops[j];
                // 该Loop中点在points中开始的索引
                size_t start = points.size();
                // 将Loop中的点添加到points中，并创建它们之间的边
                for(int p = 0; p < loop.size(); p++) {
                    // 将当前点添加到points向量中
                    points.push_back(loop[p]);
                    // 创建从当前点到下一个点的边
                    edges.push_back(CDT::Edge(points.size() - 1, points.size()));
                }
                // 移除Loop中最后添加的边，因为它是无效的
                edges.pop_back();
                // 添加最后一条边以闭合Loop
                edges.push_back(CDT::Edge(points.size() - 1, start));
            }
        }
    }
}
// 比较两个浮点数是否相等
inline bool isEqual(double a, double b) {
    return std::fabs(a - b) < SPAresabs;
}

// 比较a是否小于等于b，考虑到epsilon
inline bool isLessThanOrEqual(double a, double b) {
    return a < b || isEqual(a, b);
}

// 比较a是否大于等于b，考虑到epsilon
inline bool isGreaterThanOrEqual(double a, double b) {
    return a > b || isEqual(a, b);
}
struct scan_inter {
    double u;
    int low_inter;
    int up_inter;
};
std::unordered_map<CDT::Edge, std::vector<double>> CalculatePolygonEdgeIntersections(const point_vector& points, edge_vector edges, const edge_vector& polygon_edges, const point_vector& polygon_points, const bool is_vertical_scan) {
    using Index = int;     // 定义索引类型为整数
    using Coord = double;  // 定义坐标类型为双精度浮点数
    using Type = char;     // 定义事件类型为字符
    struct EventNode {     // 定义事件节点结构
        Coord cord;        // 坐标
        Index edge_idx;    // 边的索引
        Type type;         // 事件类型
    };
    std::vector<EventNode> scan_events;                                 // 扫描线事件列表
    std::unordered_map<Coord, std::vector<Index>> sweep_line_to_edges;  // 扫描线到边的映射
    for(int i = 0; i < edges.size(); i++) {                             // 遍历所有边
        auto [idx1, idx2] = edges[i].verts();                           // 获取边的两个顶点索引
        if(!is_vertical_scan) {                                         // 如果不是垂直扫描
            if(points[idx1].x > points[idx2].x) {                       // 如果第一个点的x坐标大于第二个点的
                std::swap(idx1, idx2);                                  // 交换两个点，确保从左到右
            }
            sweep_line_to_edges[points[idx1].y].push_back(i);  // 将边索引添加到对应的扫描线y坐标下
        } else {                                               // 如果是垂直扫描
            if(points[idx1].y > points[idx2].y) {              // 如果第一个点的y坐标大于第二个点的
                std::swap(idx1, idx2);                         // 交换两个点，确保从下到上
            }
            sweep_line_to_edges[points[idx1].x].push_back(i);  // 将边索引添加到对应的扫描线x坐标下
        }
    }
    for(auto& x: sweep_line_to_edges) {                                        // 对每个扫描线上的边进行排序
        sort(x.second.begin(), x.second.end(), [&](Index one, Index two) {     // 根据边的另一个顶点的坐标排序
            if(!is_vertical_scan) {                                            // 如果不是垂直扫描
                return points[edges[one].v1()].x < points[edges[two].v1()].x;  // 按x坐标排序
            } else {                                                           // 如果是垂直扫描
                return points[edges[one].v1()].y < points[edges[two].v1()].y;  // 按y坐标排序
            }
        });
    }
    for(int i = 0; i < polygon_edges.size(); i++) {                // 遍历多边形的所有边
        auto [idx1, idx2] = polygon_edges[i].verts();              // 获取边的两个顶点索引
        if(is_vertical_scan) {                                     // 如果是垂直扫描
            if(polygon_points[idx1].x > polygon_points[idx2].x) {  // 如果第一个点的x坐标大于第二个点的
                std::swap(idx1, idx2);                             // 交换两个点，确保从左到右
            }
            scan_events.emplace_back(polygon_points[idx1].x, i, 0);  // 添加扫描开始事件
            scan_events.emplace_back(polygon_points[idx2].x, i, 1);  // 添加扫描结束事件
        } else {                                                     // 如果不是垂直扫描
            if(polygon_points[idx1].y > polygon_points[idx2].y) {    // 如果第一个点的y坐标大于第二个点的
                std::swap(idx1, idx2);                               // 交换两个点，确保从下到上
            }
            scan_events.emplace_back(polygon_points[idx1].y, i, 0);  // 添加多边形边的起点
            scan_events.emplace_back(polygon_points[idx2].y, i, 1);  // 添加多边形边的终点
        }
    }
    for(auto& sweep_line_cord: sweep_line_to_edges) {            // 对每个扫描线坐标添加扫描事件
        scan_events.emplace_back(sweep_line_cord.first, -1, 2);  // 添加扫描线事件
    }
    std::unordered_map<CDT::Edge, std::vector<double>> result;                      // 结果映射，边到交点坐标的列表
    std::set<Index> active_edges;                                                   // 活动边集合
    sort(scan_events.begin(), scan_events.end(), [&](EventNode& a, EventNode& b) {  // 对所有事件按坐标排序
        if(isEqual(a.cord, b.cord)) {                                               // 如果坐标相等
            return a.type < b.type;                                                 // 按事件类型排序
        }
        return a.cord < b.cord;  // 按坐标排序
    });
    for(int i = 0; i < scan_events.size(); i++) {              // 遍历所有事件
        if(scan_events[i].type == 2) {                         // 如果是扫描线事件
            Coord sweep_line_cord = scan_events[i].cord;       // 获取扫描线坐标
            for(auto x: active_edges) {                        // 遍历所有活动边
                auto [idx1, idx2] = polygon_edges[x].verts();  // 获取边的两个顶点索引
                auto [x1, y1] = polygon_points[idx1];
                auto [x2, y2] = polygon_points[idx2];
                double intersect_cord = 0;  // 待计算的交点坐标
                if(!is_vertical_scan) {     // 如果不是垂直扫描
                    if(isEqual(y1, y2)) {   // 如果两个顶点的y坐标相等
                        continue;           // 跳过
                    }
                    intersect_cord = x1 + (sweep_line_cord - y1) * (x2 - x1) / (y2 - y1);  // 计算交点的x坐标
                } else {                                                                   // 如果是垂直扫描
                    if(isEqual(x1, x2)) {                                                  // 如果两个顶点的x坐标相等
                        continue;                                                          // 跳过
                    }
                    intersect_cord = y1 + (sweep_line_cord - x1) * (y2 - y1) / (x2 - x1);  // 计算交点的y坐标
                }
                auto sweep_edges = sweep_line_to_edges[sweep_line_cord];                               // 获取当前扫描线上的边
                for(auto edge: sweep_edges) {                                                          // 遍历这些边
                    double low_cord, high_cord;                                                        // 边的两个坐标范围
                    if(!is_vertical_scan) {                                                            // 如果不是垂直扫描
                        low_cord = std::min(points[edges[edge].v1()].x, points[edges[edge].v2()].x);   // 获取x坐标的最小值
                        high_cord = std::max(points[edges[edge].v1()].x, points[edges[edge].v2()].x);  // 获取x坐标的最大值
                    } else {                                                                           // 如果是垂直扫描
                        low_cord = std::min(points[edges[edge].v1()].y, points[edges[edge].v2()].y);   // 获取y坐标的最小值
                        high_cord = std::max(points[edges[edge].v1()].y, points[edges[edge].v2()].y);  // 获取y坐标的最大值
                    }
                    // 此处可以采用单调剖分 + 链表 优化
                    if(isLessThanOrEqual(low_cord, intersect_cord) && isGreaterThanOrEqual(high_cord, intersect_cord)) {  // 如果新坐标在边的坐标范围内
                        result[edges[edge]].push_back(intersect_cord);                                                    // 将新坐标添加到结果中
                    }
                }
            }
        }
        if(scan_events[i].type == 0) {                     // 如果是扫描开始事件
            active_edges.insert(scan_events[i].edge_idx);  // 将边添加到活动边集合
        }
        if(scan_events[i].type == 1) {                    // 如果是扫描结束事件
            active_edges.erase(scan_events[i].edge_idx);  // 从活动边集合中移除边
        }
    }
    return result;  // 返回结果
}
std::unordered_set<int> JudgePointsInPolygon2D(const point_vector& points, const edge_vector& polygon_edges, const point_vector& polygon_points, const bool is_vertical_scan) {
    using Index = int;
    using Coord = double;
    using Type = char;
    struct EventNode {
        Coord cord;
        Index edge_idx;
        Type type;  // 0为polygon_edge的起点  1为polygon_edge的终点 2 为扫描线
    };
    std::unordered_set<Coord> sweep_lines;                         // 扫描线
    std::unordered_map<Coord, std::vector<Coord>> intersect_cord;  // 用于存扫描线上产生的交点
    std::vector<EventNode> scan_events;                            // 扫描器

    for(const auto& p: points) {  // 获取每个点所在的扫描线
        if(is_vertical_scan) {
            sweep_lines.insert(p.x);
        } else {
            sweep_lines.insert(p.y);
        }
    }
    for(int i = 0; i < polygon_edges.size(); i++) {
        auto [idx1, idx2] = polygon_edges[i].verts();
        if(is_vertical_scan) {
            if(polygon_points[idx1].x > polygon_points[idx2].x) {
                std::swap(idx1, idx2);
            }
            scan_events.emplace_back(polygon_points[idx1].x, i, 0);
            scan_events.emplace_back(polygon_points[idx2].x, i, 1);
        } else {
            if(polygon_points[idx1].y > polygon_points[idx2].y) {
                std::swap(idx1, idx2);
            }
            scan_events.emplace_back(polygon_points[idx1].y, i, 0);
            scan_events.emplace_back(polygon_points[idx2].y, i, 1);
        }
    }
    for(auto& sweep_line_cord: sweep_lines) {
        scan_events.emplace_back(sweep_line_cord, -1, 2);  // 将每条扫描线抽象成一个点
    }
    sort(scan_events.begin(), scan_events.end(), [&](EventNode& a, EventNode& b) {
        if(isEqual(a.cord, b.cord)) {
            return a.type < b.type;  // 保证扫描线正确的位置
        }
        return a.cord < b.cord;
    });  // 将所有的点按y进行排序

    std::unordered_set<Index> active_edges;  // 维护一个当前扫描线扫到的polygon
    for(int i = 0; i < scan_events.size(); i++) {
        if(scan_events[i].type == 0) {  // 如果当前的点为polygon的边的起点，则将该边加入到当前的扫描线中
            active_edges.insert(scan_events[i].edge_idx);
        }
        if(scan_events[i].type == 2) {  // 当前的点是扫描线，则计算当前扫描线与所有相交的边的交点
            Coord sweep_line_cord = scan_events[i].cord;
            for(auto x: active_edges) {
                auto [idx1, idx2] = polygon_edges[x].verts();
                auto [x1, y1] = polygon_points[idx1];
                auto [x2, y2] = polygon_points[idx2];
                double new_cord = 0;
                if(is_vertical_scan) {
                    if(x1 == x2) {
                        continue;
                    }
                    new_cord = y1 + (sweep_line_cord - x1) * (y2 - y1) / (x2 - x1);
                } else {
                    if(y1 == y2) {
                        continue;
                    }
                    new_cord = x1 + (sweep_line_cord - y1) * (x2 - x1) / (y2 - y1);
                }
                if(sweep_line_cord == 0) {
                    int gg = 0;
                }
                intersect_cord[sweep_line_cord].push_back(new_cord);  // 由于扫描线的x或y值固定，只需要存y或x即可
            }
        }
        if(scan_events[i].type == 1) {  // 如果当前的点为polygon的边的终点，则将该边从当前的扫描线中移除
            active_edges.erase(scan_events[i].edge_idx);
        }
    }
    std::unordered_set<Index> result_idx;  // 用于存储最终的结果
    for(auto& x: intersect_cord) {         // 对所有扫描线得到的交点进行排序 （这里可以用链表维护简单多边形的连续性，就可以不用对点进行排序了）
        std::sort(x.second.begin(), x.second.end());
    }
    for(int i = 0; i < points.size(); i++) {
        int cnt = 0;
        bool is_in_polygon = false;
        if(is_vertical_scan) {
            auto& cords = intersect_cord[points[i].x];  // 获取当前的点所在的扫描线的交点集合
            // cnt = std::lower_bound(cords.begin(), cords.end(), points[i].y) - cords.begin();  // 计算在当前点左边的交点的个数（实际上就是点向左边射线产生的交点个数）
            for(auto x: cords) {
                if(isEqual(x, points[i].y)) {
                    is_in_polygon = true;
                    break;
                }
                if(x > points[i].y) {
                    break;
                }
                cnt++;
            }
        } else {
            auto& cords = intersect_cord[points[i].y];
            // cnt = std::lower_bound(cords.begin(), cords.end(), points[i].x) - cords.begin();
            for(auto x: cords) {
                if(isEqual(x, points[i].x)) {
                    is_in_polygon = true;
                    break;
                }
                if(x > points[i].x) {
                    break;
                }
                cnt++;
            }
        }
        if(is_in_polygon || !(cnt % 2)) {  // 根据cnt判断是否在面内
            result_idx.insert(i);
        }
    }
    return result_idx;
}

void ProcessEdge(point_vector& face_points, edge_vector& face_edges) {
    auto get_len = [&](const CDT::Edge& edge) {
        return (face_points[edge.v1()].x - face_points[edge.v2()].x) * (face_points[edge.v1()].x - face_points[edge.v2()].x) + (face_points[edge.v1()].y - face_points[edge.v2()].y) * (face_points[edge.v1()].y - face_points[edge.v2()].y);
    };

    sort(face_edges.begin(), face_edges.end(), [&](CDT::Edge& a, CDT::Edge& b) { return get_len(a) < get_len(b); });
    edge_vector res;
    std::vector<int> parent(face_points.size(), 0);
    std::iota(parent.begin(), parent.end(), 0);  // 从0开始填充
    std::function<int(int)> find;
    find = [&](int x) -> int {
        if(x != parent[x]) {
            parent[x] = find(parent[x]);  // 路径压缩
        }
        return parent[x];
    };
    for(auto& x: face_edges) {
        int p1 = find(x.v1());
        int p2 = find(x.v2());
        if(p1 != p2) {
            res.push_back(x);
            parent[p1] = p2;
            res.push_back(x);
        }
    }
    face_edges = std::move(res);
    return;
}

DeBugOutputUVGraph debug_graph;
void spline_faceter::ProcessFaceAndLoopsIntersection(point_vector& face_points, edge_vector& face_edges, point_vector& loops_points, edge_vector& loops_edges, point_vector& res_points, edge_vector& res_edges, CDT::EdgeUSet& markedEdge) {
    RemoveDuilplicatePointsAndEdges(face_points, face_edges);            // 移除面点和边的重复项
    edge_vector h_edges, v_edges;                                        // 水平和垂直边集
    for(int i = 0; i < face_edges.size(); i++) {                         // 遍历面的边
        auto [idx1, idx2] = face_edges[i].verts();                       // 获取边的两个顶点索引
        if(fabs(face_points[idx1].x - face_points[idx2].x) < EPSILON) {  // 判断是否为垂直边
            v_edges.push_back(face_edges[i]);                            // 添加到垂直边集
        } else {
            h_edges.push_back(face_edges[i]);  // 添加到水平边集
        }
        // debug_graph.AddconstraintEdge(std::make_pair(face_points[idx1].x, face_points[idx1].y), std::make_pair(face_points[idx2].x, face_points[idx2].y));
    }
    ProcessEdge(face_points, h_edges);
    ProcessEdge(face_points, v_edges);                                                                          // 合并内外部判断结果
    auto v_result = CalculatePolygonEdgeIntersections(face_points, v_edges, loops_edges, loops_points, true);   // 计算垂直边与环边的交点
    auto h_result = CalculatePolygonEdgeIntersections(face_points, h_edges, loops_edges, loops_points, false);  // 计算水平边与环边的交点
    point_vector new_points;                                                                                    // 新点集
    std::vector<int> new_map(face_points.size(), -1);                                                           // 点索引映射
    auto res = JudgePointsInPolygon2D(face_points, loops_edges, loops_points, true);                            // 判断点是否在多边形内部
    auto res_ = JudgePointsInPolygon2D(face_points, loops_edges, loops_points, false);                          // 判断点是否在多边形外部
    res.insert(res_.begin(), res_.end());
    for(int i = 0; i < face_points.size(); i++) {  // 遍历面的点
        if(res.find(i) == res.end()) {             // 如果点不在结果集中
            new_points.push_back(face_points[i]);  // 添加到新点集
            new_map[i] = new_points.size() - 1;    // 更新索引映射
        }
    }
    pointVec tmp;                   // 临时点向量
    for(auto& x: loops_points) {    // 遍历环点集
        tmp.push_back({x.x, x.y});  // 添加到临时点向量
    }
    auto kdtree = KDTree(tmp);                                                  // 创建KD树
    edge_vector cross_edges;                                                    // 交叉边集
    res_points = std::move(loops_points);                                       // 初始化结果点集为环点集
    size_t loops_edges_siz = loops_edges.size();                                // 多边形边集大小
    res_edges = std::move(loops_edges);                                         // 初始化结果边集为环边集
    const int loops_siz = res_points.size();                                    // 环点集大小
    res_points.insert(res_points.end(), new_points.begin(), new_points.end());  // 将新点集添加到结果点集

    for(int i = 0; i < v_edges.size(); i++) {                                                      // 遍历面的边
        auto [idx1, idx2] = v_edges[i].verts();                                                    // 获取边的两个顶点索引
        if(new_map[idx1] != -1 && new_map[idx2] != -1) {                                           // 如果两个顶点都在新点集中
            res_edges.push_back(CDT::Edge(new_map[idx1] + loops_siz, new_map[idx2] + loops_siz));  // 添加到结果边集
            // debug_graph.AddconstraintEdge(std::make_pair(face_points[idx2].x, face_points[idx2].y), std::make_pair(face_points[idx1].x, face_points[idx1].y));

        } else if(new_map[idx1] != -1 || new_map[idx2] != -1) {  // 如果其中一个顶点在新点集中
            if(new_map[idx1] == -1) {                            // 如果第一个顶点不在新点集中
                std::swap(idx1, idx2);                           // 交换两个顶点
            }
            auto v_cords = v_result[v_edges[i]];                              // 获取交点坐标
            for(auto cord: v_cords) {                                         // 遍历交点坐标,每个cord代表一个y值
                point_t temp = {face_points[idx1].x, cord};                   // 创建临时点
                auto near_idx = kdtree.nearest_index(temp);                   // 查找最近点索引
                res_edges.emplace_back(new_map[idx1] + loops_siz, near_idx);  // 添加到结果边集
                // debug_graph.AddconstraintEdge(std::make_pair(face_points[idx1].x, face_points[idx1].y), std::make_pair(res_points[near_idx].x, res_points[near_idx].y));
            }
        } else {
            // 如果两个顶点都不在新点集中
            if(!v_result[v_edges[i]].empty()) {
                auto cords = v_result[v_edges[i]];
                point_t temp1 = {face_points[idx1].x, cords.front()};
                auto near_idx1 = kdtree.nearest_index(temp1);
                point_t temp2 = {face_points[idx2].x, cords.back()};
                auto near_idx2 = kdtree.nearest_index(temp2);
                res_edges.emplace_back(near_idx1, near_idx2);
                // debug_graph.AddconstraintEdge(std::make_pair(res_points[near_idx2].x, res_points[near_idx2].y), std::make_pair(res_points[near_idx1].x, res_points[near_idx1].y));
            }
        }
    }
    for(int i = 0; i < h_edges.size(); i++) {                                                      // 遍历面的边
        auto [idx1, idx2] = h_edges[i].verts();                                                    // 获取边的两个顶点索引
        if(new_map[idx1] != -1 && new_map[idx2] != -1) {                                           // 如果两个顶点都在新点集中
            res_edges.push_back(CDT::Edge(new_map[idx1] + loops_siz, new_map[idx2] + loops_siz));  // 添加到结果边集
        } else if(new_map[idx1] != -1 || new_map[idx2] != -1) {                                    // 如果其中一个顶点在新点集中
            if(new_map[idx1] == -1) {                                                              // 如果第一个顶点不在新点集中
                std::swap(idx1, idx2);                                                             // 交换两个顶点
            }
            auto h_cords = h_result[h_edges[i]];                              // 获取交点坐标
            for(auto cord: h_cords) {                                         // 遍历交点坐标,每个cord代表一个x值
                point_t temp = {cord, face_points[idx1].y};                   // 创建临时点
                auto near_idx = kdtree.nearest_index(temp);                   // 查找最近点索引
                res_edges.emplace_back(new_map[idx1] + loops_siz, near_idx);  // 添加到结果边集
                // debug_graph.AddconstraintEdge(std::make_pair(face_points[idx1].x, face_points[idx1].y), std::make_pair(res_points[near_idx].x, res_points[near_idx].y));
            }
        } else {
            // 如果两个顶点都不在新点集中
            if(!h_result[h_edges[i]].empty()) {  // 如果

                auto cords = h_result[h_edges[i]];
                point_t temp1 = {cords.front(), face_points[idx1].y};
                auto near_idx1 = kdtree.nearest_index(temp1);
                point_t temp2 = {cords.back(), face_points[idx2].y};
                auto near_idx2 = kdtree.nearest_index(temp2);
                res_edges.emplace_back(near_idx1, near_idx2);
                // debug_graph.AddconstraintEdge(std::make_pair(res_points[near_idx2].x, res_points[near_idx2].y), std::make_pair(res_points[near_idx1].x, res_points[near_idx1].y));
            }
        }
    }

    CDT::RemoveDuplicatesAndRemapEdges(res_points, res_edges);              // 移除重复的边并重新映射边的索引
    for(int i = 0; i < loops_edges_siz; i++) {                              // 遍历环边集
        markedEdge.insert({res_edges[i].v1() + 3, res_edges[i].v2() + 3});  // 标记边
    }
}

void spline_faceter::attachMeshByQuad(QuadTree& tree) {
#ifdef DEBUG

    debug_graph.init(f);
    for(auto& x: tree.result_points_) {
        debug_graph.AddSamplePoint(std::make_pair(x.x, x.y));
    }
    for(auto& x: points) {
        debug_graph.AddBoundaryPoint(std::make_pair(x.x, x.y));
    }
    for(auto& x: edges) {
        debug_graph.AddconstraintEdge(std::make_pair(points[x.v1()].x, points[x.v1()].y), std::make_pair(points[x.v2()].x, points[x.v2()].y));
    }

#endif  // DEBUG

    // 将四叉树的面采样点集合并到现有点集中
    CDT::EdgeUSet markedEdge;
    edge_vector final_edges;
    point_vector final_points;

    TIMERSTART(ProcessFaceAndLoopsIntersection)
    ProcessFaceAndLoopsIntersection(tree.result_points_, tree.result_edges_, points, edges, final_points, final_edges, markedEdge);

    TIMEREND(ProcessFaceAndLoopsIntersection)
    // 创建CDT三角剖分对象
    TIMERSTART(CDT_TRI)
    CDT::Triangulation<double>* pcdt = new CDT::Triangulation<double>;  // CDT进行三角剖分
    // 插入顶点
    pcdt->insertVertices(final_points);
    // 插入边
    pcdt->insertEdges(final_edges);
    auto&& mp = CDT::EdgeToPiecesMapping(pcdt->pieceToOriginals);
    for(auto& x: mp) {
        if(markedEdge.find(x.first) != markedEdge.end()) {
            for(auto& v: x.second) {
                markedEdge.insert(v);
            }
            // debug_graph.AddconstraintEdge(std::make_pair(final_points[x.first.v1() - 3].x, final_points[x.first.v1() - 3].y), std::make_pair(final_points[x.first.v2() - 3].x, final_points[x.first.v2() - 3].y));
        }
    }
    pcdt->eraseOuterTrianglesAndHoles(markedEdge);
    TIMEREND(CDT_TRI)
    // 对三角形进行后处理
    TIMERSTART(PostProcessTriangle)
    PostProcessTriangle(f, nt, st, pcdt);
    TIMEREND(PostProcessTriangle)
    // 获取三角形和顶点数据
    const CDT::TriangleVec triangles = std::move(pcdt->triangles);
    const std::vector<CDT::V2d<double>> vertice = std::move(pcdt->vertices);

    // 创建索引网格对象
    INDEXED_MESH* mesh = new INDEXED_MESH( vertice.size(), triangles.size(), triangles.size() * 3);
    // 遍历所有顶点，添加到网格中
    for(int i = 0; i < vertice.size(); i++) {
        CDT::V2d<double> pp = vertice[i];
        SPApar_pos param(pp.x, pp.y);
        // 计算顶点位置和法线
        const SPAposition pos = f->geometry()->equation().eval_position(param);
        SPAunit_vector normal = f->geometry()->equation().eval_normal(param);
        // 如果面的方向相反，则反转法线
        if(f->sense()) normal = -normal;
        // 添加顶点到网格
        mesh->add_vertex( pos, normal, param);
    }
#ifdef DEBUG

    auto Edges_set = CDT::extractEdgesFromTriangles(triangles);
    debug_graph.AddTriangles(vertice, Edges_set);
    debug_graph.output();

#endif  // DEBUG

    // 遍历所有三角形，添加到网格中
    for(int i = 0; i < triangles.size(); i++) {
        int p1 = triangles[i].vertices[0];
        int p2 = triangles[i].vertices[1];
        int p3 = triangles[i].vertices[2];
        // 添加多边形到网格
        mesh->add_polygon( i, 3);
        indexed_polygon* poly0 = mesh->get_polygon( i);
        // 设置多边形的顶点
        poly0->set_vertex( 0, &mesh->get_vertex( p1));
        poly0->set_vertex( 1, &mesh->get_vertex( p2));
        poly0->set_vertex( 2, &mesh->get_vertex( p3));
    }

    // 将索引网格附加到面上
    attach_indexed_mesh_to_face(f, mesh);
}

inline const double spline_faceter::GetFitParam(const double param) {
    if(equal(param, u_le)) {
        return u_le;
    }
    if(equal(param, u_ri)) {
        return u_ri;
    }
    if(equal(param, v_le)) {
        return v_le;
    }
    if(equal(param, v_ri)) {
        return v_ri;
    }
    return param;
}
void spline_faceter::attachMesh() {
    SPApar_pos param;

    /*SHELL* shell = (SHELL*)f->owner();
    LUMP* lump = (LUMP*)shell->owner();
    BODY* body = (BODY*)lump->owner();*/
    SPAtransf transf;
    auto body = get_owner(f);
    if(!is_BODY(body)) {
        transf = SPAtransf();
    } else {
        TRANSFORM* TRANS = ((BODY*)body)->transform();
        if(TRANS != NULL) transf = TRANS->transform();
    }

    surface* surf = &f->geometry()->equation_for_update();

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

            if(results != NULL) {
                if(results->next != NULL) {
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
    do {
        if(pcdt) delete pcdt;
        pcdt = new CDT::Triangulation<double>;
        CDT::RemoveDuplicatesAndRemapEdges(points, edges);

        // for(int i = 0; i < points.size(); i++) {
        //     curr = points[i];
        //     printf("%lf %lf\n", curr.x, curr.y);
        // }
        // printf("%d\n", points.size());
        // for(int i = 0; i < edges.size(); i++) {
        //     printf("%d %d\n", edges[i].v1(), edges[i].v2());
        // }

        pcdt->insertVertices(points);
        pcdt->insertEdges(edges);
        pcdt->eraseOuterTrianglesAndHoles();
        _count += 1;
        // points = pcdt->vertices;
    } while(_count < 4 && !fixTriangle(f, nt, st, pcdt, &points));

    CDT::TriangleVec triangles = pcdt->triangles;
    std::vector<CDT::V2d<double>> vertice = pcdt->vertices;

    INDEXED_MESH* mesh = new INDEXED_MESH( vertice.size(), pcdt->triangles.size(), pcdt->triangles.size() * 3);
    for(int i = 0; i < vertice.size(); i++) {
        CDT::V2d<double> pp = vertice[i];
        SPApar_pos param(pp.x, pp.y);
        SPAposition pos = f->geometry()->equation().eval_position(param);
        // SPAposition pos = SPAposition(pp.x * 10.0, pp.y * 10.0, 0.0f);
        SPAunit_vector normal = f->geometry()->equation().eval_normal(param);
        if(f->sense()) normal = -normal;
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
outcome spline_faceter::facet() {
    // 初始化
    init();
    // 计算du,dv长度
    decideUVlen();
    TIMERSTART(getLoops);
    bool has_periphery_loop = getLoops();
    TIMEREND(getLoops);
    if(has_periphery_loop) {
        peripheryProcess();
        /*auto curr = periphery[0];
        auto next = periphery[1];
        if((equal(curr.x, ri_high.x) || equal(curr.x, le_low.x)) && (equal(next.x, ri_high.x) || equal(next.x, le_low.x))) {
            if(curr.y > next.y) {
                periphery[0].x = le_low.x;
                periphery[1].x = le_low.x;
            } else {
                periphery[0].x = ri_high.x;
                periphery[1].x = ri_high.x;
            }
        }

        if((equal(curr.y, ri_high.y) || equal(curr.y, le_low.y)) && (equal(next.y, ri_high.y) || equal(next.y, le_low.y))) {
            if(curr.x < next.x) {
                periphery[0].y = le_low.y;
                periphery[1].y = le_low.y;
            } else {
                periphery[0].y = ri_high.y;
                periphery[1].y = ri_high.y;
            }
        }*/
    } else {
        // 若没有边缘Loop，自己生成边集
        // int index1, index2;
        // for(int i = 0; i <= vlen; i++) {
        //    CDT::V2d<double> tempL = {u_le, v_le + i * dv};
        //    CDT::V2d<double> tempH = {u_ri, v_le + i * dv};
        //    points.push_back(tempL);
        //    points.push_back(tempH);
        //    index1 = points.size() - 1;
        //    index2 = points.size() - 2;
        //    edges.push_back(CDT::Edge(index2, index2 + 2));
        //    edges.push_back(CDT::Edge(index1, index1 + 2));
        //    // if(i < vlen) {
        //    //     printf("%lf %lf %lf %lf\n", u_le, v_le + i * dv, u_le, v_le + (i + 1) * dv);
        //    //     printf("%lf %lf %lf %lf\n", u_ri, v_le + i * dv, u_ri, v_le + (i + 1) * dv);
        //    // }
        //}
        // edges.pop_back();
        // edges.pop_back();
        // for(int i = 0; i <= ulen; i++) {
        //    CDT::V2d<double> tempL = {u_le + i * du, v_le};
        //    CDT::V2d<double> tempH = {u_le + i * du, v_ri};
        //    points.push_back(tempL);
        //    points.push_back(tempH);
        //    index1 = points.size() - 1;
        //    index2 = points.size() - 2;
        //    edges.push_back(CDT::Edge(index2, index2 + 2));
        //    edges.push_back(CDT::Edge(index1, index1 + 2));
        //    // if(i < ulen) {
        //    //     printf("%lf %lf %lf %lf\n", u_le + i * du, v_le, u_le + (i + 1) * du, v_le);
        //    //     printf("%lf %lf %lf %lf\n", u_le + i * du, v_le, u_le + (i + 1) * du, v_le);
        //    // }
        //}
        // edges.pop_back();
        // edges.pop_back();
    }

    if(u_loops.size() > 0) {
        USeperationProcess();
    }

    if(v_loops.size() > 0) {
        VSeperationProcess();
    }

    if(hole_loops.size() > 0) {
        holeProcess();
    }
    // attachMesh();
    std::vector<SPApar_vec> vec;
    for(auto& x: hole_loops) {
        for(auto& y: x) {
            vec.push_back(SPApar_vec(y.x, y.y));
        }
    }
    for(auto& x: v_loops) {
        for(auto& y: x) {
            vec.push_back(SPApar_vec(y.x, y.y));
        }
    }
    for(auto& x: u_loops) {
        for(auto& y: x) {
            vec.push_back(SPApar_vec(y.x, y.y));
        }
    }
    for(auto& y: periphery) {
        vec.push_back(SPApar_vec(y.x, y.y));
    }
    // QuadTree tree(f, nt, st, vec);
    TIMERSTART(QUADTREE);
    QuadTree tree(f, nt, st);
    TIMEREND(QUADTREE);
    attachMeshByQuad(tree);
    return outcome();
}

bool singular_u(FACE* f, SPAparameter u, SPAinterval vrange) {
    SPApar_pos param1 = {u, vrange.mid_pt()};
    SPApar_pos param2 = {u, vrange.start_pt()};

    SPAposition pos1 = f->geometry()->equation().eval_position(param1);
    SPAposition pos2 = f->geometry()->equation().eval_position(param2);
    if(pos1 == pos2) {
        return true;
    }
    return false;
}

outcome gme_facet_face_spline(FACE* f, double nt, double st) {
    spline_faceter faceter(f, nt, st);
    return faceter.facet();
}
// 获取periphery点集, 其他loop点集
bool spline_faceter::getLoops() {
    SPApar_pos param;
    bool has_loop_periphery = false;
    for(LOOP* lp = f->loop(); lp; lp = lp->next()) {
        std::vector<CDT::V2d<double>> tempPS;
        loop_type lptype;
        // 获取loop类型
        api_loop_type(lp, lptype);
        unsigned ic = 0;
        bool periphery_has_singularity = false;

        for(COEDGE* coedge = lp->start(); coedge != lp->start() || ic == 0; coedge = coedge->next(), ic++) {
            point_vector tempps;
            SPAposition* polyline;
            double* params;
            EDGE* edge = coedge->edge();
            REVBIT r = coedge->sense();  // 判断edge的指向与coedge指向是否相反，关系到点添加顺序问题
            // sg_add_pcurve_to_coedge(coedge);

            int nP = 0;
            // 从边离散的结果中获取离散后的三维点集
            get_facet_edge_points_and_params(coedge->edge(), polyline, params, nP);
            // api_get_facet_edge_points();
            sg_add_pcurve_to_coedge(coedge);
            PCURVE* PC = coedge->geometry();
            if(PC == nullptr) {
                param = f->geometry()->equation().param(polyline[0]);
                param = {GetFitParam(param.u), GetFitParam(param.u)};
                tempps.push_back(point(param.u, le_low.y));
                tempps.push_back(point(param.u, ri_high.y));
                if(equal(param.u, ri_high.x)) {
                    std::reverse(tempps.begin(), tempps.end());
                }
            } else {
                pcurve pc = PC->equation();
                for(int i = 0; i < nP; i++) {
                    SPApar_pos param1 = f->geometry()->equation().param(polyline[i]);
                    CDT::V2d<double> temp = {GetFitParam(param1.u), GetFitParam(param1.v)};

                    // 此处处理是由于边离散给了重复点导致的
                    // if(tempps.size() != 0 && temp == *tempps.rbegin()) continue;
                    tempps.push_back(temp);
                }
            }
            if(r) std::reverse(tempps.begin(), tempps.end());

            // 此处取奇点处最近似的点的uv值近似替代，细分会使得该近似的误差减小
            singular_type start = is_singular(tempps[0]);
            singular_type startNext = is_singular(tempps[1]);
            point_boundary_type startN = is_on_boundary(tempps[1]);
            singular_type end = is_singular(tempps[tempps.size() - 1]);
            singular_type endPre = is_singular(tempps[tempps.size() - 2]);
            point_boundary_type endP = is_on_boundary(tempps[tempps.size() - 2]);
            // 如果第一个点为奇点，下个点非奇点并且不在边界上
            if(start != no_singular && startNext == no_singular && startN == not_boundary) {
                switch(start) {
                    case ule:
                    case uri:
                        tempps[0].y = tempps[1].y;  // 调整奇点位置与非奇点位置平行
                        break;
                    case vle:
                    case vri:
                        tempps[0].x = tempps[1].x;  // 同理
                        break;
                }
            }
            if(end != no_singular && endPre == no_singular && endP == not_boundary) {  // 同上
                switch(end) {
                    case ule:
                    case uri:
                        tempps[tempps.size() - 1].y = tempps[tempps.size() - 2].y;
                        break;
                    case vle:
                    case vri:
                        tempps[tempps.size() - 1].x = tempps[tempps.size() - 2].x;
                        break;
                }
            }

            // std::cout << "tempss_begin:" << tempps.front().x << ' ' << tempps.front().y << std::endl;
            tempPS.insert(tempPS.end(), tempps.begin(), tempps.end());  // 插入该loop点集中
        }
        // 通过法线与切线叉积的向量映射到参数域空间与loop序做叉积判断loop点序
        bool normalRev = false;
        if(tempPS.size() > 2) {
            point a, b;
            SPAposition a3, b3;
            for(int i = 0; i < tempPS.size(); i++) {
                a = tempPS[i];
                b = tempPS[(i + 1) % tempPS.size()];
                a3 = spl->eval_position(SPApar_pos(a.x, a.y));
                b3 = spl->eval_position(SPApar_pos(b.x, b.y));
                cross_type cu = cross_boundary(a.x, b.x, u_tolerance);
                cross_type cv = cross_boundary(a.y, b.y, v_tolerance);
                if(a3 == b3 || cu != no_cross || cv != no_cross) continue;
                break;
            }
            SPApar_pos a2(a.x, a.y);
            SPAvector ab3 = b3 - a3;
            SPAunit_vector an = spl->eval_normal(SPApar_pos(a.x, a.y));
            if(f->sense()) an = -an;
            SPAvector dir3 = an * ab3;
            SPAunit_vector udir3(dir3.x(), dir3.y(), dir3.z());
            SPApar_dir dir = spl->param_dir(udir3, a2);
            // std::cout << dir.du << " " << dir.dv << std::endl << std::endl;
            double temp = (b.x - a.x) * dir.dv - (b.y - a.y) * dir.du;
            if((b.x - a.x) * dir.dv - (b.y - a.y) * dir.du > 0.0) {
                normalRev = true;
            }
        }

        // if(normalRev) reverse(tempPS.begin(), tempPS.end());

        // 不同Loop类型插入不同类型点集合中
        switch(lptype) {
            case loop_type::loop_hole: {
                hole_loops.push_back(tempPS);
                // printf("hole\n");
                break;
            }
            case loop_type::loop_periphery: {
                periphery = tempPS;
                has_loop_periphery = true;

                // printf("periphery\n");
                break;
            }
            case loop_type::loop_u_separation: {
                // printf("loop_u_separation\n");
                //  if(!f->geometry()->equation().closed_v()) {
                //      periphery = tempPS;
                //      has_loop_periphery = true;
                //  } else {
                u_loops.push_back(tempPS);
                // }
                break;
            }
            case loop_type::loop_v_separation: {
                v_loops.push_back(tempPS);
                // printf("loop_v_separation\n");
                break;
            }
            case loop_type::loop_uv_separation: {
                uv_loops.push_back(tempPS);
                // printf("loop_uv_separation\n");
                break;
            }
            default: {
                // printf("unknown loop\n");
            }
        }
    }

    return has_loop_periphery;
}
