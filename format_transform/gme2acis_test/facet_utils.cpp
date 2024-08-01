#include "acis/gme/faceter/gme_facet_utils.hxx"

#include <acis/include/faceutil.hxx>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

#include "acis/gme/faceter/gme_idx_mesh_utils.hxx"
#include "acis/include/add_pcu.hxx"
#include "acis/include/af_api.hxx"
#include "acis/include/body.hxx"
#include "acis/include/condef.hxx"
#include "acis/include/curdef.hxx"
#include "acis/include/curve.hxx"
#include "acis/include/elldef.hxx"
#include "acis/include/getowner.hxx"
#include "acis/include/idx_mesh.hxx"
#include "acis/include/intrapi.hxx"
#include "acis/include/lists.hxx"
#include "acis/include/meshat.hxx"
#include "acis/include/pcudef.hxx"
#include "acis/include/pcurve.hxx"
#include "acis/include/pladef.hxx"
#include "acis/include/plane.hxx"
#include "acis/include/position.hxx"
#include "acis/include/ptlist.hxx"
#include "acis/include/sp3srtn.hxx"
#include "acis/include/sphdef.hxx"
#include "acis/include/spldef.hxx"
#include "acis/include/surdef.hxx"
#include "acis/include/surface.hxx"
#include "acis/include/torus.hxx"
#include "acis/include/transfrm.hxx"

logical get_facet_edge_points_and_params1(EDGE* edge, SPAposition*& pos_array, double*& param_array, int& num_pts) {
    AF_POINT* start;
    AF_POINT* end;
    AF_POINT* curr;
    int knt = 0;
    pos_array = nullptr;
    param_array = nullptr;

    // Find the list of AF_POINTS on the given edge, if one exists.
    if(AF_POINT::find(edge, 0, start, end)) {
        // Determine the number of points in the list.
        for(curr = start; curr != end; curr = curr->next(0)) {
            knt++;
        }
        knt++;

        // Allocate arrays of the proper size.
        pos_array = ACIS_NEW SPAposition[knt];
        param_array = ACIS_NEW double[knt];

        // Populate the arrays.
        int index = 0;
        for(curr = start; index < knt; curr = curr->next(0)) {
            pos_array[index] = curr->get_position();
            param_array[index] = curr->get_parameter();
            index++;
        }
    }
    num_pts = knt;
    return TRUE;
}

void get_edge_tolerance(double in_body_diagnoal, facet_options_visualization* fo, double& out_edge_distance_tolerance, double& out_edge_angle_tolerance) {
    double out_face_distance_tolerance = 0, out_face_angle_tolerance = 0;
    get_face_tolerance(in_body_diagnoal, fo, out_face_distance_tolerance, out_face_angle_tolerance);
    edge_quality_level this_edge_quality = fo->get_edge_quality();
    if(this_edge_quality == medium) {
        out_edge_distance_tolerance = out_face_distance_tolerance;
        out_edge_angle_tolerance = out_face_angle_tolerance;
    } else if(this_edge_quality == better) {  // default
        out_edge_distance_tolerance = out_face_distance_tolerance / 2.0;
        out_edge_angle_tolerance = out_face_angle_tolerance / 2.0;
    } else if(this_edge_quality == best) {
        out_edge_distance_tolerance = out_face_distance_tolerance / 4.0;
        out_edge_angle_tolerance = out_face_angle_tolerance / 4.0;
    }
}

void get_face_tolerance(double in_body_diagnoal, facet_options_visualization* fo, double& out_face_distance_tolerance, double& out_face_angle_tolerance) {
    face_quality_level this_face_quality = fo->get_face_quality();
    if(this_face_quality == coarse) {
        out_face_distance_tolerance = in_body_diagnoal * 0.004;
        out_face_angle_tolerance = 40.0;
    } else if(this_face_quality == medium_coarse) {  // default
        out_face_distance_tolerance = in_body_diagnoal * 0.002;
        out_face_angle_tolerance = 30.0;
    } else if(this_face_quality == medium_fine) {
        out_face_distance_tolerance = in_body_diagnoal * 0.001;
        out_face_angle_tolerance = 20.0;
    } else if(this_face_quality == fine) {
        out_face_distance_tolerance = in_body_diagnoal * 0.0005;
        out_face_angle_tolerance = 10.0;
    }
}

void FacetMesh::PeriodicFaceUVRepair(FACE* face, FacetTestMesh& face_data) {
    bool closed_u = face->geometry()->equation().closed_u();
    bool closed_v = face->geometry()->equation().closed_v();
    if(!closed_u && !closed_v) return;

    for(int i = 0; i < face_data.triangles.size(); i++) {
        FixTrianglePeriodic(face, face_data.triangles[i]);
    }
}

void _print_curve(curve* c) {
    SPAinterval c_range = c->param_range();
    printf("period %lf\n", c->param_period());
    printf("start %lf end %lf\n", c_range.start_pt(), c_range.end_pt());
    if(c->type() == 2) {
        ellipse* ec = (ellipse*)c;
        printf("type ellipse\n");
        printf("center %lf %lf %lf\n", ec->centre.x(), ec->centre.y(), ec->centre.z());
        printf("major_axis %lf %lf %lf\n", ec->major_axis.x(), ec->major_axis.y(), ec->major_axis.z());
        printf("major_axis length %lf \n", ec->major_axis_length);
        printf("normal %lf %lf %lf \n", ec->normal.x(), ec->normal.y(), ec->normal.z());
        printf("param_off %lf\n", ec->param_off);
        printf("radius_ratio %lf \n", ec->radius_ratio);
    }
}
void _print_edge(EDGE* e) {
    printf("closed %d periodic %d\n", e->closed(), e->periodic());
    SPAinterval p_range = e->param_range();
    printf("start %lf end %lf\n", p_range.start_pt(), p_range.end_pt());
    printf("startp %lf endp %lf\n", e->start_param(), e->end_param());
}
void _print_entity(ENTITY* ent) {
    BODY* body = (BODY*)ent;
    LUMP* lmp = body->lump();
    for(int ilmp = 0; lmp != nullptr; lmp = lmp->next(), ilmp++) {
        printf("LUMP %d\n", ilmp);
        SHELL* shl = lmp->shell();
        for(int is = 0; shl != nullptr; shl = shl->next(), is++) {
            printf(" SHELL %d\n", is);
            FACE* fc = shl->face();
            for(int ifc = 0; fc != nullptr; fc = fc->next(), ifc++) {
                printf("  FACE %d\n", ifc);
                std::cout << "Face Type? (1: Plane; 2: Cone; 3: Sphere; 4: Torus; 10: Spline): " << fc->geometry()->equation().type() << std::endl;

                SPAinterval u_interval = fc->geometry()->equation().param_range_u();
                SPAinterval v_interval = fc->geometry()->equation().param_range_v();
                std::cout << "principle u range:" << u_interval.start_pt() << " " << u_interval.end_pt() << std::endl;
                std::cout << "principle v range:" << v_interval.start_pt() << " " << v_interval.end_pt() << std::endl;
                // std::cout << fc->geometry()->equation().singular_u(u_interval.end_pt());
                if(fc->geometry()->equation().type() == 11 || fc->geometry()->equation().type() == 4) {
                    // SPAunit_vector nor = fc->geometry()->equation().eval_normal()

                    curve* u_le = fc->geometry()->equation().v_param_line(u_interval.end_pt());
                    SPAinterval interval = u_le->param_range();
                    // std::cout << interval.start_pt() << " " << interval.end_pt() << std::endl;
                    std::cout << u_interval.end_pt() << "\t" << u_le->length(interval.start_pt(), interval.end_pt()) << std::endl;

                    u_le = fc->geometry()->equation().v_param_line(0.99999999);
                    interval = u_le->param_range();
                    // std::cout << interval.start_pt() << " " << interval.end_pt() << std::endl;
                    std::cout << "0.99999999"
                              << "\t" << u_le->length(interval.start_pt(), interval.end_pt()) << std::endl;

                    u_le = fc->geometry()->equation().v_param_line(0.9999);
                    interval = u_le->param_range();
                    // std::cout << interval.start_pt() << " " << interval.end_pt() << std::endl;
                    std::cout << 0.9999 << "\t" << u_le->length(interval.start_pt(), interval.end_pt()) << std::endl;

                    SPAunit_vector normal = fc->geometry()->equation().eval_normal(SPApar_pos(u_interval.end_pt(), v_interval.start_pt()));
                    SPAunit_vector outdir = fc->geometry()->equation().eval_outdir(SPApar_pos(u_interval.end_pt(), v_interval.start_pt()));
                    if(fc->sense()) normal = -normal;
                    std::cout << "normal:" << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;
                    std::cout << "outdir:" << outdir.x() << " " << outdir.y() << " " << normal.z() << std::endl;
                    // spline* spl = (spline*)&fc->geometry()->equation();
                    // bs3_surface sur = spl->sur();
                    // if(sur != nullptr) {
                    //     std::cout << "bs3_surface exists" << std::endl;
                    //     std::cout << "bs3_surface_singular_u: " << bs3_surface_singular_u(u_interval.end_pt(), sur) << std::endl;
                    // }
                }
                printf("\tface sense: %d", fc->sense());
                printf("\tclosed_u: %d\t closed_v: %d\n", fc->geometry()->equation().closed_u(), fc->geometry()->equation().closed_v());
                LOOP* lp = fc->loop();
                for(int ilp = 0; lp != nullptr; lp = lp->next(), ilp++) {
                    loop_type lptype;
                    api_loop_type(lp, lptype);  // 获取loop的类型
                    printf("   LOOP %d %d(4: u_seperation)\n", ilp, lptype);
                    COEDGE *cedge = lp->start(), *fcedge = cedge;
                    for(int ic = 0; cedge != fcedge || ic == 0; cedge = cedge->next(), ++ic) {
                        printf("    COEDGE %d\n", ic);
                        REVBIT r = cedge->sense();
                        REVBIT curR = cedge->edge()->sense();
                        printf("\tcoedge_sense: %d edge_sense: %d\n", r, curR);
                        SPAposition point = cedge->start_pos();
                        printf("\tstart_point: (%.2f, %.2f, %.2f)\t%d\n", point.x(), point.y(), point.z(), cedge->starts_at_singularity());
                        point = cedge->end_pos();
                        printf("\tend_point: (%.2f, %.2f, %.2f)\t%d\n", point.x(), point.y(), point.z(), cedge->ends_at_singularity());
                        int _typeid = fc->geometry()->equation().type();
                        if(_typeid == 4 || _typeid == 10) {
                            int nP = 0;
                            SPAposition* polyline;
                            double* params;
                            get_facet_edge_points_and_params1(cedge->edge(), polyline, params, nP);
                            torus* tori = (torus*)&fc->geometry()->equation_for_update();
                            printf("reverse-v: %d\n", tori->reverse_v);
                            sg_add_pcurve_to_coedge(cedge);
                            PCURVE* PC = cedge->geometry();
                            if(PC == nullptr) continue;
                            pcurve pc = PC->equation();
                            //
                            if(!r) {  // 基于指向进行不同处理
                                for(int i = 0; i < nP; i++) {
                                    SPApar_pos param1 = pc.eval_position(polyline[i], params[i], 1);
                                    SPApar_pos param = fc->geometry()->equation().param(polyline[i]);
                                    SPAvector normal = fc->geometry()->equation().eval_normal(param);
                                    printf("3D %lf %lf %lf\t", polyline[i].x(), polyline[i].y(), polyline[i].z());
                                    // printf("param %lf\t", params[i]);
                                    // printf("w2D %lf %lf\t", param1.u, param1.v);
                                    // printf("2D %lf %lf\n", param.u, param.v);

                                    // printf("normal %lf %lf %lf\n", normal.x(), normal.y(), normal.z());
                                    // if(tori->degenerate() && lptype == loop_type::loop_u_separation) param = coedge->geometry()->equation().eval_position(polyline[i], 1);  // @todo: 为该接口第二个参数寻找一个合适的参数par

                                    printf("%lf %lf\n", param.u, param.v);
                                }
                            } else {
                                for(int i = nP - 1; i >= 0; i--) {
                                    SPApar_pos param1 = pc.eval_position(polyline[i], params[i], 1);
                                    SPApar_pos param = fc->geometry()->equation().param(polyline[i]);
                                    SPAvector normal = fc->geometry()->equation().eval_normal(param);
                                    printf("3D %lf %lf %lf\t", polyline[i].x(), polyline[i].y(), polyline[i].z());
                                    // printf("param %lf\t", params[i]);
                                    // printf("w2D %lf %lf\t", param1.u, param1.v);
                                    // printf("2D %lf %lf\n", param.u, param.v);

                                    printf("%lf %lf\n", param.u, param.v);
                                    // printf("%lf %lf\n", param.u, param.v);
                                    // printf("%lf %lf\n", param.u, param.v);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void _print_mesh_result(ENTITY* ent) {
    bool if_print_position = FALSE;  // 选择是否输出点的坐标值，默认为否

    printf("\t=====Begin print the MESH.=====\n");
    ENTITY_LIST eg_list;
    get_edges(ent, eg_list);
    printf("\t=====EDGE's num: %d=====\n", eg_list.count());
    EDGE* e = (EDGE*)eg_list.first();
    for(int i_e = 0; i_e < eg_list.count(); i_e++) {
        if(e->length() != 0) printf("===EDGE[%d]: \033[47;30m%s\033[0m\n", i_e, e->geometry()->equation().type_name());
        SPAposition* pos = nullptr;
        int n;
        api_get_facet_edge_points(e, pos, n);
        printf("EDGE[%d] has %d points.\n", i_e, n);
        if(if_print_position) {
            for(int j = 0; j < n; j++) {
                printf("point[%2d]: (%5.2f,%5.2f,%5.2f)\n", j, pos[j].x(), pos[j].y(), pos[j].z());
            }
        }
        e = (EDGE*)eg_list.next();
    }

    ENTITY_LIST f_list;
    get_faces(ent, f_list);
    printf("\t=====FACE's num: %d=====\n", f_list.count());
    FACE* f = (FACE*)f_list.first();
    INDEXED_MESH* mesh;
    for(int i_f = 0; i_f < f_list.count(); i_f++, f = (FACE*)f_list.next()) {
        printf("===FACE[%d]: \033[47;30m%s\033[0m\n", i_f, f->geometry()->equation().type_name());
        mesh = GetIndexedMesh(f);
        if(mesh != nullptr) {
            int num_polygon = mesh->get_num_polygon();
            printf("mesh[%d] has %d polygons.\n", i_f, num_polygon);
            if(if_print_position) {
                for(int i_plg = 0; i_plg < num_polygon; i_plg++) {
                    int num_vertex = mesh->get_polygon(i_plg)->num_vertex();
                    printf("mesh[%d] - polygon[%d] has %d vertex.\n", i_f, i_plg, num_vertex);
                    printf("mesh[%d] - polygon[%d]", i_f, i_plg);
                    for(int i_v = 0; i_v < num_vertex; i_v++) {
                        SPAposition pos = mesh->get_polygon(i_plg)->get_vertex(i_v)->get_position();
                        printf("-(%5.2f,%5.2f,%5.2f)", pos.x(), pos.y(), pos.z());
                    }
                    printf("\n");
                }
            }
        }
    }
}

double PointToEdgeDist(SPAposition P, SPAposition A, SPAposition B) {
    SPAvector AB = B - A;
    SPAvector AP = P - A;
    double AB_squared = AB % AB;
    double t = AP % AB / AB_squared;
    if(t < 0) return AP.len();
    if(t > 1) return (P - B).len();
    return (AP - AB * t).len();
}

double PointToTriangleDist(SPAposition P, SPAposition A, SPAposition B, SPAposition C) {
    if(A == B && B == C)
        return (P - A).len();
    else if(A == B)
        return PointToEdgeDist(P, B, C);
    else if(B == C)
        return PointToEdgeDist(P, C, A);
    else if(C == A)
        return PointToEdgeDist(P, A, B);

    SPAvector AB = B - A;
    SPAvector AC = C - A;

    SPAunit_vector n = normalise(AB * AC);
    SPAvector AP = P - A;
    double d = n % AP;

    SPAposition Q = P - d * n;

    SPAvector AQ = Q - A;
    SPAvector BC = C - B;
    SPAvector BQ = Q - B;
    SPAvector CQ = Q - C;

    SPAvector cross1 = AQ * AB;
    SPAvector cross2 = BQ * BC;
    SPAvector cross3 = CQ * (-AC);

    if(cross1 % cross2 > 0 && cross2 % cross3 > 0) return fabs(d);

    double d1 = PointToEdgeDist(P, A, B);
    double d2 = PointToEdgeDist(P, B, C);
    double d3 = PointToEdgeDist(P, C, A);
    return std::min({d1, d2, d3});
}

double PointToTriangleDist(SPAposition P, FacetTestTriangle t) {
    return PointToTriangleDist(P, t.pt[0].pos, t.pt[1].pos, t.pt[2].pos);
}

double PointToMeshDistance(SPAposition P, INDEXED_MESH* mesh) {
    int num_polygon = mesh->get_num_polygon();
    double min_dis = DBL_MAX;
    for(int i = 0; i < num_polygon; i++) {
        FacetTestTriangle t(mesh->get_polygon(i));
        min_dis = std::min(min_dis, PointToTriangleDist(P, t));
    }
    return min_dis;
}

double HausdorffDistance(INDEXED_MESH* mesh1, INDEXED_MESH* mesh2) {
    double max_dis = 0;
    int num_vertex1 = mesh1->get_num_vertex();
    int num_vertex2 = mesh2->get_num_vertex();
    for(int i = 0; i < num_vertex1; i++) {
        FacetTestPoint pt(mesh1->get_vertex(i));
        max_dis = std::max(max_dis, PointToMeshDistance(pt.pos, mesh2));
    }
    for(int i = 0; i < num_vertex2; i++) {
        FacetTestPoint pt(mesh2->get_vertex(i));
        max_dis = std::max(max_dis, PointToMeshDistance(pt.pos, mesh1));
    }
    return max_dis;
}

FacetInfo::FacetInfo(const FacetInfo& obj) {
    total_error = obj.total_error;
    max_error = obj.max_error;
    min_error = obj.min_error;
    triangle_number = obj.triangle_number;
    vertex_number = obj.vertex_number;
    time = obj.time;
}

// FacetInfo& FacetInfo::operator=(const FacetInfo& obj) {
//     if(this == &obj) return *this;
//     total_error = obj.total_error;
//     max_error = obj.max_error;
//     min_error = obj.min_error;
//     triangle_number = obj.triangle_number;
//     vertex_number = obj.vertex_number;
//     time = obj.time;
//     return *this;
// }

// FacetInfo& FacetInfo::operator+=(const FacetInfo& other) {
//     total_error += other.total_error;
//     max_error = std::max(max_error, other.max_error);
//     min_error = std::min(min_error, other.min_error);
//     triangle_number += other.triangle_number;
//     vertex_number += other.vertex_number;
//     time += other.time;
//     return *this;
// }

// std::ostream& operator<<(std::ostream& os, const FacetInfo& obj) {
//     os << obj.total_error << "," << obj.max_error << "," << obj.min_error << "," << obj.triangle_number << "," << obj.vertex_number << "," << obj.time << ",";
//     return os;
// }

FaceterEvaluator::FaceterEvaluator(ENTITY* entity, bool use_acis) {
    entity_ = entity;
    use_acis_ = use_acis;
    use_ = true;
};

FaceterEvaluator::FaceterEvaluator(ENTITY* entity) {
    entity_ = entity;
    use_acis_ = true;
    use_ = true;
};

FaceterEvaluator::FaceterEvaluator() {
    entity_ = (ENTITY*)NULL_REF;
    use_acis_ = false;
    use_ = false;
};

FaceterEvaluator::~FaceterEvaluator() {
}

void FaceterEvaluator::FacetEntity() {
    INDEXED_MESH* idx_mesh = NULL;
    ENTITY* owner;
    api_get_owner(entity_, owner);

    if(is_BODY(owner)) {
        BODY* body;
        body = dynamic_cast<BODY*>(owner);
        if(body->transform()) transf_ = body->transform()->transform();
    }

    std::chrono::steady_clock::time_point start, end;
    std::chrono::microseconds time, time_gme;
    if(use_acis_) {
        start = std::chrono::high_resolution_clock::now();
        api_facet_entity(entity_);
        end = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    }
    if(use_gme_) {
        start = std::chrono::high_resolution_clock::now();
        api_facet_entity(entity_);
        end = std::chrono::high_resolution_clock::now();
        time_gme = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    }
    ENTITY_LIST face_list;
    get_faces(entity_, face_list);
    FACE* face = dynamic_cast<FACE*>(face_list.first());
    for(int i = 0; i < face_list.count(); i++, face = dynamic_cast<FACE*>(face_list.next())) {
        if(use_acis_) {
            MESH* face_mesh = NULL;
            af_query(static_cast<ENTITY*>(face), IDX_MESH_APP, IDX_MESH_ID, face_mesh);
            idx_mesh = dynamic_cast<INDEXED_MESH*>(face_mesh);
            face_error_map_.emplace(face, FacetInfo());
            face_error_map_.at(face).time = time.count();
            face_mesh_map_.emplace(face, idx_mesh);
        }
        if(use_gme_) {
            // std::cout << time.count() << std::endl;
            idx_mesh = GetIndexedMesh(face);
            face_error_map_.emplace(face, FacetInfo());
            face_error_map_.at(face).time = time_gme.count();
            face_mesh_map_.emplace(face, idx_mesh);
        }
    }
}

void FaceterEvaluator::ComputeAllError() {
    for(auto it: face_mesh_map_) {
        // std::cout << it.first->geometry()->equation().type_name() << std::endl;
        // WriteMeshInObjFormat("test.obj", it.second);
        if(use_acis_) ComputeErrorBySample(it.first, 0);
        if(use_) ComputeErrorBySample(it.first, 1);
        face_hausdorff_map_.emplace(it.first, HausdorffDistance(face_mesh_map_.at(it.first), face_mesh_map_.at(it.first)));
    }
}

void FaceterEvaluator::ComputeError(FACE* face, int engine) {
    // 三角片面的误差为三角面片参数域重心在空间中的位置与三角面片重心位置的距离，以便于计算
    // 该误差应当不为最大误差，但是通过大数定律取平均似乎可以得到一个较好的结果
    // 通过上述误差乘以该三角面片在uv参数域上的面积，累计求和，得到该face的总的误差
    // 以下注释均为调试信息，后续可删除
    const int kPolygonNodeNb = 3;
    if(engine == 0) {
        face_error_map_.at(face).triangle_number = face_mesh_map_.at(face)->get_num_polygon();
        face_error_map_.at(face).vertex_number = face_mesh_map_.at(face)->get_num_vertex();
        face_mesh_map_.at(face)->set_par_pos_mapping_01(false);
        FacetTestMesh mesh(face_mesh_map_.at(face));
        FacetMesh::PeriodicFaceUVRepair(face, mesh);
        for(int i = 0; i < mesh.triangles.size(); i++) {
            FacetTestTriangle t = mesh.triangles[i];
            SPApar_pos center_2d = t.pt[0].uv + (t.pt[1].uv - SPApar_pos(.0, .0) + t.pt[2].uv - SPApar_pos(.0, .0));
            center_2d.u /= kPolygonNodeNb;
            center_2d.v /= kPolygonNodeNb;
            SPAposition pos = face->geometry()->equation().eval_position(center_2d) * transf_;
            double error = PointToTriangleDist(pos, t.pt[0].pos * transf_, t.pt[1].pos * transf_, t.pt[2].pos * transf_);
            face_error_map_.at(face).max_error = std::max(face_error_map_.at(face).max_error, error);
            face_error_map_.at(face).min_error = std::min(face_error_map_.at(face).min_error, error);
            SPApar_vec v01(t.pt[1].uv - t.pt[0].uv), v02(t.pt[2].uv - t.pt[0].uv);
            // 以三角面片参数域上的面积作为权重
            face_error_map_.at(face).total_error += fabs(v01 * v02) / 2 * error;
        }
    }
    if(engine == 1) {
        int polygon_nb = face_mesh_map_.at(face)->get_num_polygon();
        face_error_map_.at(face).triangle_number = polygon_nb;
        face_error_map_.at(face).vertex_number = face_mesh_map_.at(face)->get_num_vertex();

        face_mesh_map_.at(face)->set_par_pos_mapping_01(false);
        for(int i = 0; i < polygon_nb; i++) {
            int vertex_nb = face_mesh_map_.at(face)->get_polygon(i)->num_vertex();
            SPAposition center_3d(.0, .0, .0);
            SPApar_pos center_2d(.0, .0);
            SPAposition pos_3d[kPolygonNodeNb];
            SPApar_pos pos_2d[kPolygonNodeNb];
            // SPApar_pos pos_2d2[kPolygonNodeNb];
            assert(kPolygonNodeNb == vertex_nb);
            for(int j = 0; j < vertex_nb; j++) {
                pos_3d[j] = face_mesh_map_.at(face)->get_polygon(i)->get_vertex(j)->get_position() * transf_;
                center_3d += pos_3d[j] - SPAposition(.0, .0, .0);
                // int index = face_mesh_map_.at(face)->get_vertex_index(face_mesh_map_.at(face)->get_polygon(i)->get_vertex(j));
                pos_2d[j] = face_mesh_map_.at(face)->get_polygon(i)->get_vertex(j)->get_uv();
                // SPApar_pos param = face_mesh_map_.at(face)->get_uv_as_entered(index);
                // SPApar_pos param2 = face_mesh_map_.at(face)->get_uv_as_scaled(index);
                center_2d += pos_2d[j] - SPApar_pos(.0, .0);
            }
            for(int j = 0; j < 3; j++) {
                center_3d.set_coordinate(j, (center_3d.coordinate(j)) / vertex_nb);
            }
            center_2d.u /= vertex_nb;
            center_2d.v /= vertex_nb;
            SPAposition pos = face->geometry()->equation().eval_position(center_2d) * transf_;
            // SPAposition tmp = face->geometry()->equation().eval_position(pos_2d[0]);
            SPApar_pos tmpp = face->geometry()->equation().param(pos_3d[0]);
            // api_apply_transf((ENTITY*)&tmp, (*transf_).transform_data);
            double error = PointToTriangleDist(pos, pos_3d[0], pos_3d[1], pos_3d[2]);
            PointToTriangleDist(pos, pos_3d[0], pos_3d[1], pos_3d[2]);
            face_error_map_.at(face).max_error = std::max(face_error_map_.at(face).max_error, error);
            face_error_map_.at(face).min_error = std::min(face_error_map_.at(face).min_error, error);
            SPApar_vec v01(pos_2d[1] - pos_2d[0]), v02(pos_2d[2] - pos_2d[0]);
            // 以三角面片参数域上的面积作为权重
            double tmp = fabs(v01 * v02) / 2;
            face_error_map_.at(face).total_error += fabs(v01 * v02) / 2 * error;
        }
    }
}

void FaceterEvaluator::ComputeErrorBySample(FACE* face, int engine) {
    // std::cout << face->geometry()->equation().type_name() << std::endl;
    point_face_containment pfc;
    double du, dv;
    SPApar_box box;
    sg_get_face_par_box(face, box);
    int u_sample = 30;
    int v_sample = 30;
    du = box.u_range().length() / u_sample;
    dv = box.v_range().length() / v_sample;
    for(int i = 0; i < u_sample; i++) {
        for(int j = 0; j < v_sample; j++) {
            SPApar_pos pos(box.u_range().start_pt() + i * du, box.v_range().start_pt() + j * dv);
            SPAposition p = face->geometry()->equation().eval_position(pos);
            api_point_in_face(p, face, SPAtransf(), pfc, pos);
            if(pfc == point_face_containment::point_inside_face) {
                double error = 0;
                if(engine == 0) {
                    error = PointToMeshDistance(p, face_mesh_map_.at(face));
                    face_error_map_.at(face).max_error = std::max(face_error_map_.at(face).max_error, error);
                    face_error_map_.at(face).min_error = std::min(face_error_map_.at(face).min_error, error);
                }
                if(engine == 1) {
                    error = PointToMeshDistance(p, face_mesh_map_.at(face));
                    face_error_map_.at(face).max_error = std::max(face_error_map_.at(face).max_error, error);
                    face_error_map_.at(face).min_error = std::min(face_error_map_.at(face).min_error, error);
                }
            }
        }
    }
    if(engine == 0) {
        face_error_map_.at(face).triangle_number = face_mesh_map_.at(face)->get_num_polygon();
        face_error_map_.at(face).vertex_number = face_mesh_map_.at(face)->get_num_vertex();
    }
    if(engine == 1) {
        face_error_map_.at(face).triangle_number = face_mesh_map_.at(face)->get_num_polygon();
        face_error_map_.at(face).vertex_number = face_mesh_map_.at(face)->get_num_vertex();
    }
}

bool FacetMesh::IsPeriodicU(FACE* face, FacetTestPoint pt, PointBoundaryType1D& out) {
    int type = face->geometry()->equation().type();
    if(type == 1 || type == 2 || type == 3) {
        // 平面 圆锥 球面
        out = PointBoundaryType1D::kMiddle;
        return false;
    }
    bool closed_u = face->geometry()->equation().closed_u();
    double lower = face->geometry()->equation().param_range_u().start_pt();
    double upper = face->geometry()->equation().param_range_u().end_pt();
    if(closed_u) {
        out = GetPointBoundaryType1D(pt.uv.u, lower, upper);
        if(out == PointBoundaryType1D::kOnBoundary) return true;
    } else
        out = PointBoundaryType1D::kMiddle;
    return false;
}

bool FacetMesh::IsPeriodicV(FACE* face, FacetTestPoint pt, PointBoundaryType1D& out) {
    int type = face->geometry()->equation().type();
    if(type == 1) {
        // 平面
        out = PointBoundaryType1D::kMiddle;
        return false;  // 平面不需要判断
    }
    bool closed_v = face->geometry()->equation().closed_v();
    double lower = face->geometry()->equation().param_range_v().start_pt();
    double upper = face->geometry()->equation().param_range_v().end_pt();
    if(closed_v) {
        out = GetPointBoundaryType1D(pt.uv.v, lower, upper);
        if(out == PointBoundaryType1D::kOnBoundary) return true;
    } else
        out = PointBoundaryType1D::kMiddle;
    return false;
}

bool FacetMesh::FixTrianglePeriodic(FACE* face, FacetTestTriangle& t) {
    int type = face->geometry()->equation().type();
    if(type == 1) return false;
    bool u_flag = true, v_flag = true;

    const int kTriangleNodeNb = 3;
    PointBoundaryType1D out[kTriangleNodeNb];
    bool is_periodic[kTriangleNodeNb];
    int cnt = 0;
    int idx = 0;
    for(int i = 0; i < kTriangleNodeNb; i++) {
        is_periodic[i] = IsPeriodicU(face, t.pt[i], out[i]);
        if(is_periodic[i]) cnt += 1;
    }
    if(cnt == 1) {
        for(int i = 0; i < kTriangleNodeNb; i++) {
            if(is_periodic[i]) {
                idx = i;
                break;
            }
        }
        if(out[(idx + 1) % 3] != out[(idx + 2) % 3]) {
            assert(out[(idx + 1) % 3] == out[(idx + 2) % 3]);
        }
        if(out[(idx + 1) % 3] == PointBoundaryType1D::kUpper) {
            t.pt[idx].uv.u = face->geometry()->equation().param_range_u().end_pt();
        } else if(out[(idx + 1) % 3] == PointBoundaryType1D::kLower) {
            t.pt[idx].uv.u = face->geometry()->equation().param_range_u().start_pt();
        }
    } else if(cnt == 2) {
        for(int i = 0; i < kTriangleNodeNb; i++) {
            if(!is_periodic[i]) {
                idx = i;
                break;
            }
        }
        if(out[idx] == PointBoundaryType1D::kUpper) {
            t.pt[(idx + 1) % 3].uv.u = face->geometry()->equation().param_range_u().end_pt();
            t.pt[(idx + 2) % 3].uv.u = face->geometry()->equation().param_range_u().end_pt();
        } else if(out[idx] == PointBoundaryType1D::kLower) {
            t.pt[(idx + 1) % 3].uv.u = face->geometry()->equation().param_range_u().start_pt();
            t.pt[(idx + 2) % 3].uv.u = face->geometry()->equation().param_range_u().start_pt();
        }
    } else
        u_flag = false;

    cnt = 0;
    for(int i = 0; i < kTriangleNodeNb; i++) {
        is_periodic[i] = IsPeriodicV(face, t.pt[i], out[i]);
        if(is_periodic[i]) cnt += 1;
    }
    if(cnt == 1) {
        for(int i = 0; i < kTriangleNodeNb; i++) {
            if(is_periodic[i]) {
                idx = i;
                break;
            }
        }
        // assert(out[(idx + 1) % 3] == out[(idx + 2) % 3]);
        if(out[(idx + 1) % 3] == PointBoundaryType1D::kUpper) {
            t.pt[idx].uv.v = face->geometry()->equation().param_range_v().end_pt();
        } else if(out[(idx + 1) % 3] == PointBoundaryType1D::kLower) {
            t.pt[idx].uv.v = face->geometry()->equation().param_range_v().start_pt();
        }
    } else if(cnt == 2) {
        for(int i = 0; i < kTriangleNodeNb; i++) {
            if(!is_periodic[i]) {
                idx = i;
                break;
            }
        }
        if(out[idx] == PointBoundaryType1D::kUpper) {
            t.pt[(idx + 1) % 3].uv.v = face->geometry()->equation().param_range_v().end_pt();
            t.pt[(idx + 2) % 3].uv.v = face->geometry()->equation().param_range_v().end_pt();
        } else if(out[idx] == PointBoundaryType1D::kLower) {
            t.pt[(idx + 1) % 3].uv.v = face->geometry()->equation().param_range_v().start_pt();
            t.pt[(idx + 2) % 3].uv.v = face->geometry()->equation().param_range_v().start_pt();
        }
    } else
        v_flag = false;

    if(u_flag || v_flag) return true;
}

bool FacetMesh::IsApex(FACE* face, const FacetTestPoint& pt, SingularityType& out) {
    out = SingularityType::kNotSingularity;
    int type = face->geometry()->equation().type();
    if(type == 1 || type == 10) {
        return false;
    } else if(type == 2) {
        // 圆锥
        cone* c = dynamic_cast<cone*>(&face->geometry()->equation_for_update());
        if(!c->cylinder() && c->get_apex() == pt.pos) {
            if(pt.uv.u == c->param_range_u().end_pt())
                out = SingularityType::kUri;
            else if(pt.uv.u == c->param_range_u().start_pt())
                out = SingularityType::kUle;
            return true;
        }
    } else if(type == 3) {
        sphere* s = dynamic_cast<sphere*>(&face->geometry()->equation_for_update());
        SPAposition north_pole = s->centre + s->pole_dir * s->radius;
        SPAposition south_pole = s->centre - s->pole_dir * s->radius;
        if(north_pole == pt.pos) {
            out = SingularityType::kUri;
            return true;
        } else if(south_pole == pt.pos) {
            out = SingularityType::kUle;
            return true;
        }
    } else if(type == 4) {
        torus* t = dynamic_cast<torus*>(&face->geometry()->equation_for_update());
        if(t->degenerate()) {
            SPAposition north_apex = t->centre + t->normal * t->apex_dist();
            SPAposition south_apex = t->centre - t->normal * t->apex_dist();
            if(north_apex == pt.pos) {
                out = SingularityType::kUri;
                return true;
            } else if(south_apex == pt.pos) {
                out = SingularityType::kUle;
                return true;
            }
        }
    } else if(type == 10) {
    }
    return false;
}

const FacetInfo FaceterEvaluator::total_info() const {
    FacetInfo info;
    for(const auto& it: face_error_map_) {
        info.total_error += it.second.total_error;
        info.max_error = std::max(info.max_error, it.second.max_error);
        info.min_error = std::min(info.min_error, it.second.min_error);
        info.triangle_number += it.second.triangle_number;
        info.vertex_number += it.second.vertex_number;
        info.time = it.second.time;
    }
    return info;
};

const FacetInfo FaceterEvaluator::total_info_gme() const {
    FacetInfo info;
    for(const auto& it: face_error_map_gme_) {
        info.total_error += it.second.total_error;
        info.max_error = std::max(info.max_error, it.second.max_error);
        info.min_error = std::min(info.min_error, it.second.min_error);
        info.triangle_number += it.second.triangle_number;
        info.vertex_number += it.second.vertex_number;
        info.time = it.second.time;
    }
    return info;
};

const double FaceterEvaluator::hausdorff_distance() const {
    double res = 0.0;
    for(const auto& it: face_hausdorff_map_) {
        res += it.second;
    }
    return res;
}

void FaceterEvaluator::Evaluate() {
    FacetEntity();
    ComputeAllError();
}

static void get_polylines_from_faceted_edges(ENTITY_LIST& edges, std::vector<FacetMesh::EdgeData>& edges_data, bool use_acis) {
    int nE = 0;
    int nV = 0;
    int numEdges = edges.iteration_count();
    for(ENTITY* ent = edges.first(); ent; ent = edges.next()) {
        assert(nE < numEdges);
        assert(is_EDGE(ent));
        if(!is_EDGE(ent)) {
            continue;
        }

        SPAtransf tr = get_owner_transf(ent);
        EDGE* edge = (EDGE*)ent;
        SPAposition* pos = nullptr;

        FacetMesh::EdgeData edge_data;
        int nP = 0;
        outcome out;
        if(use_acis)
            out = api_get_facet_edge_points(edge, pos, nP);
        else
            out = api_get_facet_edge_points(edge, pos, nP);
        if(!out.ok()) {
            ACIS_DELETE[] pos;
            continue;
        }
        for(int ii = 0; ii < nP; ii++) {
            pos[ii] *= tr;
            edge_data.coords.push_back((float)pos[ii].x());
            edge_data.coords.push_back((float)pos[ii].y());
            edge_data.coords.push_back((float)pos[ii].z());
        }
        ACIS_DELETE[] pos;
        pos = nullptr;
        edges_data.push_back(edge_data);

        nV += 3 * nP;
        nE++;
    }
}
static void get_triangles_from_faceted_face_withUV(FACE* face, FacetMesh::FaceData& face_data) {
    af_serializable_mesh* sm = GetSerializableMesh(face);
    if(nullptr == sm) {
        // Application decision: do we throw for unfaceted faces?
        return;
    }
    SPAtransf tr = get_owner_transf(face);

    const int nv = sm->number_of_vertices();
    int ntri = sm->number_of_polygons();

    face_data.coords.resize(3 * nv);
    sm->serialize_positions(face_data.coords.data());  // if std::vector::data is not available, &(coords[0]) will also work.

    bool const has_normals = sm->has_normals() == TRUE;
    if(has_normals) {
        face_data.normal_coords.resize(3 * nv);
        sm->serialize_normals(face_data.normal_coords.data());
    }

    if(!tr.identity()) {
        for(int ii = 0; ii < nv; ii++) {
            int jj = 3 * ii;

            SPAposition pos(face_data.coords[jj], face_data.coords[jj + 1], face_data.coords[jj + 2]);
            SPAvector normal(face_data.normal_coords[jj], face_data.normal_coords[jj + 1], face_data.normal_coords[jj + 2]);
            pos *= tr;
            normal = normalise(normal * tr);
            face_data.coords[jj] = pos.x();
            face_data.coords[jj + 1] = pos.y();
            face_data.coords[jj + 2] = pos.z();
            face_data.normal_coords[jj] = normal.x();
            face_data.normal_coords[jj + 1] = normal.y();
            face_data.normal_coords[jj + 2] = normal.z();
        }
    }

    bool const has_uv = sm->has_uv() == TRUE;
    if(has_uv) {
        face_data.uv_coords.resize(2 * nv);
        sm->serialize_uv_data(face_data.uv_coords.data(), true);
    }

    face_data.triangles.resize(3 * ntri);
    int ntri_actual = sm->serialize_triangles(face_data.triangles.data());
    while(ntri_actual < ntri) {
        face_data.triangles.pop_back();
        ntri_actual = static_cast<int>(face_data.triangles.size());
    }

    face_data.vertex_Nb = nv;
    face_data.triangle_Nb = ntri;
}
static void get_triangles_from_faceted_faces_withUV(FACE* face, FacetMesh::FaceData& face_data) {
    INDEXED_MESH* mesh = GetIndexedMesh(face);
    if(nullptr == mesh) {
        // Application decision: do we throw for unfaceted faces?
        return;
    }
    SPAtransf tr = get_owner_transf(face);

    const int nv = mesh->number_of_vertices();
    int ntri = mesh->number_of_polygons();

    face_data.coords.resize(3 * nv);
    mesh->serialize_positions( face_data.coords.data());

    bool const has_normals = mesh->has_normals();
    if(has_normals) {
        face_data.normal_coords.resize(3 * nv);

        mesh->serialize_normals( face_data.normal_coords.data());
    }

    if(!tr.identity()) {
        for(int ii = 0; ii < nv; ii++) {
            int jj = 3 * ii;

            SPAposition pos(face_data.coords[jj], face_data.coords[jj + 1], face_data.coords[jj + 2]);
            SPAvector normal(face_data.normal_coords[jj], face_data.normal_coords[jj + 1], face_data.normal_coords[jj + 2]);
            pos *= tr;
            normal = normalise(normal * tr);
            face_data.coords[jj] = pos.x();
            face_data.coords[jj + 1] = pos.y();
            face_data.coords[jj + 2] = pos.z();
            face_data.normal_coords[jj] = normal.x();
            face_data.normal_coords[jj + 1] = normal.y();
            face_data.normal_coords[jj + 2] = normal.z();
        }
    }

    bool const has_uv = mesh->has_uv();
    if(has_uv) {
        face_data.uv_coords.resize(2 * nv);
        for(int i = 0; i < nv; i++) {
            SPApar_pos pos = mesh->get_uv_as_entered(i);
            face_data.uv_coords[2 * i] = pos.u;
            face_data.uv_coords[2 * i + 1] = pos.v;
        }
    }

    face_data.triangles.resize(3 * ntri);
    int ntri_actual = mesh->serialize_triangles( face_data.triangles.data());
    while(ntri_actual < ntri) {
        face_data.triangles.pop_back();
        ntri_actual = static_cast<int>(face_data.triangles.size());
    }
}

void WriteMeshInObjFormat(const std::string& filename, FACE* face) {
    std::ofstream file(filename);

    if(!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    FacetMesh::FaceData face_data;
    get_triangles_from_faceted_face_withUV(face, face_data);

    for(int i = 0; i < face_data.vertex_Nb; i++) {
        file << "v " << face_data.coords[3 * i] << " " << face_data.coords[3 * i + 1] << " " << face_data.coords[3 * i + 2] << std::endl;
    }
    for(int i = 0; i < face_data.vertex_Nb; i++) {
        file << "vt " << face_data.uv_coords[2 * i] << " " << face_data.uv_coords[2 * i + 1] << std::endl;
    }
    for(int i = 0; i < face_data.vertex_Nb; i++) {
        file << "vn " << face_data.normal_coords[3 * i] << " " << face_data.normal_coords[3 * i + 1] << " " << face_data.normal_coords[3 * i + 2] << std::endl;
    }
    for(int i = 0; i < face_data.triangle_Nb; i++) {
        file << "f " << face_data.triangles[3 * i] + 1 << "/" << face_data.triangles[3 * i] + 1 << "/" << face_data.triangles[3 * i] + 1 << " ";
        file << face_data.triangles[3 * i + 1] + 1 << "/" << face_data.triangles[3 * i + 1] + 1 << "/" << face_data.triangles[3 * i + 1] + 1 << " ";
        file << face_data.triangles[3 * i + 2] + 1 << "/" << face_data.triangles[3 * i + 2] + 1 << "/" << face_data.triangles[3 * i + 2] + 1 << std::endl;
    }

    file.close();
}

void WriteMeshInObjFormat(const std::string& filename, INDEXED_MESH* mesh) {
    std::ofstream file(filename);

    if(!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    int number_triangles = mesh->number_of_polygons();
    int number_vertices = mesh->number_of_vertices();

    for(int i = 0; i < number_vertices; i++) {
        SPAposition pos = mesh->get_position(i);
        file << "v " << pos.x() << " " << pos.y() << " " << pos.z() << std::endl;
    }
    for(int i = 0; i < number_vertices; i++) {
        SPApar_pos pos = mesh->get_uv_as_entered(i);
        file << "vt " << pos.u << " " << pos.v << std::endl;
    }
    for(int i = 0; i < number_vertices; i++) {
        SPAunit_vector normal = mesh->get_normal(i);
        file << "vn " << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;
    }
    for(int i = 0; i < number_triangles; i++) {
        indexed_polygon* poly = mesh->get_polygon(i);
        int node_count = poly->num_vertex();
        file << "f";
        for(int node_i = 0; node_i < node_count; node_i++) {
            polygon_vertex* vert = poly->get_vertex(node_i);
            int index = mesh->get_vertex_index(vert);
            file << " " << index + 1 << "/" << index + 1 << "/" << index + 1;
        }
        file << std::endl;
    }

    file.close();
}

void WriteMeshInObjFormat(const std::string& filename, ENTITY_LIST& faces, bool use_acis) {
    std::ofstream file(filename);

    if(!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    int nF = 0;
    int nV = 0;
    int nI = 0;
    int numFaces = faces.iteration_count();
    std::vector<FacetMesh::FaceData> entity_data;
    entity_data.resize(numFaces);

    // 记录各个面的数量
    std::vector<int> counter(5, 0);
    for(ENTITY* ent = faces.first(); ent; ent = faces.next()) {
        // assert(nF < numFaces);
        assert(is_FACE(ent));
        if(!is_FACE(ent)) {
            continue;
        }

        FacetMesh::FaceData face_data;
        FACE* face = dynamic_cast<FACE*>(ent);
        switch(face->geometry()->equation().type()) {
            case 1:
                file << "o plane" << counter[0] << std::endl;
                counter[0]++;
                break;
            case 2:
                file << "o cone" << counter[1] << std::endl;
                counter[1]++;
                break;
            case 3:
                file << "o sphere" << counter[2] << std::endl;
                counter[2]++;
                break;
            case 4:
                file << "o torus" << counter[3] << std::endl;
                counter[3]++;
                break;
            case 10:
                file << "o spline" << counter[4] << std::endl;
                counter[4]++;
                break;
            default:
                break;
        }
        if(use_acis)
            get_triangles_from_faceted_face_withUV(face, face_data);
        else
            get_triangles_from_faceted_faces_withUV(face, face_data);

        for(int i = 0; i < face_data.vertex_Nb; i++) {
            file << "v " << face_data.coords[3 * i] << " " << face_data.coords[3 * i + 1] << " " << face_data.coords[3 * i + 2] << std::endl;
        }
        for(int i = 0; i < face_data.vertex_Nb; i++) {
            file << "vt " << face_data.uv_coords[2 * i] << " " << face_data.uv_coords[2 * i + 1] << std::endl;
        }
        for(int i = 0; i < face_data.vertex_Nb; i++) {
            file << "vn " << face_data.normal_coords[3 * i] << " " << face_data.normal_coords[3 * i + 1] << " " << face_data.normal_coords[3 * i + 2] << std::endl;
        }
        for(int i = 0; i < face_data.triangle_Nb; i++) {
            file << "f " << face_data.triangles[3 * i] + 1 + nV << "/" << face_data.triangles[3 * i] + 1 + nV << "/" << face_data.triangles[3 * i] + 1 + nV << " ";
            file << face_data.triangles[3 * i + 1] + 1 + nV << "/" << face_data.triangles[3 * i + 1] + 1 + nV << "/" << face_data.triangles[3 * i + 1] + 1 + nV << " ";
            file << face_data.triangles[3 * i + 2] + 1 + nV << "/" << face_data.triangles[3 * i + 2] + 1 + nV << "/" << face_data.triangles[3 * i + 2] + 1 + nV << std::endl;
        }
        file << std::endl;

        nI += (unsigned int)face_data.triangle_Nb;
        nV += (unsigned int)face_data.vertex_Nb;
        nF++;
    }

    file.close();
}

void WriteEntityInObjFormat(const std::string& filename, ENTITY* entity, bool use_acis) {
    if(use_acis)
        api_facet_entity(entity);
    else
        api_facet_entity(entity);
    ENTITY_LIST faces;
    get_faces(entity, faces);
    WriteMeshInObjFormat(filename, faces, use_acis);
}

void WriteEntitiesInObjFormat(const std::string& filename, ENTITY_LIST& entities, bool use_acis) {
    ENTITY* owner;
    api_get_owner(entities.first(), owner);
    if(use_acis)
        api_facet_entities(owner, &entities);
    else {
        for(auto entity: entities) {
            api_facet_entity(entity);
        }
    }
    ENTITY_LIST faces;
    for(auto entity: entities) {
        api_get_faces(entity, faces);
    }
    WriteMeshInObjFormat(filename, faces, use_acis);
}

void WritePolylineInObjFormat(const std::string& filename, ENTITY_LIST& edges, bool use_acis) {
    std::ofstream file(filename);

    if(!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::vector<FacetMesh::EdgeData> edges_data;
    get_polylines_from_faceted_edges(edges, edges_data, use_acis);

    int cnt = 0;
    for(auto& edge_data: edges_data) {
        for(int i = 0; i < edge_data.coords.size(); i += 3) {
            file << "v " << edge_data.coords[i] << " " << edge_data.coords[i + 1] << " " << edge_data.coords[i + 2] << std::endl;
        }
        for(int i = 0; i < edge_data.coords.size() / 3; i++) {
            file << "l " << cnt + i + 1 << " " << cnt + i + 2 << std::endl;
        }
        cnt += edge_data.coords.size() / 3;
        file << std::endl;
    }

    file.close();
}

void WriteEntityLineInObjFormat(const std::string& filename, ENTITY* entity, bool use_acis) {
    if(use_acis)
        api_facet_entity(entity);
    else
        api_facet_entity(entity);
    ENTITY_LIST edges;
    api_get_edges(entity, edges);
    WritePolylineInObjFormat(filename, edges, use_acis);
}

void WriteEntitiesLineInObjFormat(const std::string& filename, ENTITY_LIST& entities, bool use_acis) {
    ENTITY* owner;
    api_get_owner(entities.first(), owner);
    if(use_acis)
        api_facet_entities(owner, &entities);
    else {
        for(auto entity: entities) {
            api_facet_entity(entity);
        }
    }
    ENTITY_LIST edges;
    for(auto entity: entities) {
        api_get_edges(entity, edges);
    }
    WritePolylineInObjFormat(filename, edges, use_acis);
}