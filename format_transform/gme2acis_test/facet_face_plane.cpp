#include "acis/gme/faceter/gme_facet_face_plane.hxx"

#include <array>
#include <vector>

#include "acis/include/curdef.hxx"
#include "acis/include/curve.hxx"

/**
 * @brief 将给定的平面使用约束的Delaunay三角剖分（CDT）离散成三角形网格
 *
 * @param f 一个指针指向要离散的平面对象
 * @param nt 法向容差
 * @param st 距离容差
 * @return 一个表示操作成功或失败的outcome对象
 */
outcome gme_facet_face_plane(FACE* f, double nt, double st) {
    // points变量用于存放边离散后的点
    // 格式为((u1,v1),(u2,v2),...)
    std::vector<CDT::V2d<double>> points;

    // edges变量用于存放连接点的边
    // 格式为((e1,e2),(e2,e3),...)
    // 其中(ei,ej)的ei为首顶点下标,ej为尾顶点下标
    std::vector<CDT::Edge> edges;

    // beg变量用于记录每个环的开始顶点下标
    int beg = 0;
    // loop_num变量用于记录当前处理的环的编号
    int loop_num = 0;

    // 获取面的法向
    REVBIT r = f->sense();

    // 开始遍历面的所有环
    for(LOOP* lp = f->loop(); lp != nullptr; lp = lp->next()) {
        // n_sample变量记录每个loop中得到的所有的点的个数
        int n_sample = 0;

// 在调试模式下输出环的编号
#if FACETER_DEBUG_MODE
        printf("loop_num: %d\n", loop_num);
#endif
        // loop_num++;

        // cedge变量是一个loop中的起始coedge指针
        COEDGE *cedge = lp->start(), *fcedge = cedge;

        // 使用for循环获取loop中的所有边
        for(int ic = 0; cedge != fcedge || ic == 0; cedge = cedge->next(), ic = 1) {
            // eg变量为coedge下的edge指针
            EDGE* eg = cedge->edge();

            // out变量用于记录coedge的方向是否和edge的方向相同
            // 相同为0，不同为1
            REVBIT out = cedge->sense();

            // flag变量是一个标识变量, 用于标识加入边离散后的点到points变量中是否顺序加入
            // flag=1表示顺序加入, flag=0表示逆序加入
            int flag = 0;

            // 如果coedge与基础edge方向相同，则顺序加入点COEDGE与基础EDGE方向相同
            if(!out) {
                flag = 1;
            }

            // 判断eg的形状, 如果eg只有一个点, 需要进行特殊处理
            if(eg->geometry() == nullptr) {
                // n_sample1变量记录边离散后点的个数
                int n_sample1;
                // pl1变量记录一组边离散后的顶点
                SPAposition* pl1;

                // 离散边获取点存入pl1
                api_get_facet_edge_points(eg, pl1, n_sample1);

                // 将n_sample置为1, 表示这个loop中只有一个点
                n_sample = 1;
                SPAposition pos = pl1[0];
                // 将得到的点的参数坐标加入points变量中, 并退出此次循环
                // 将得到的点加入points变量中, 并退出此次循环
                points.push_back({f->geometry()->equation().param(pos).u, f->geometry()->equation().param(pos).v});

                continue;
            }

// 调试模式下输出当前边的几何类型名称
#if FACETER_DEBUG_MODE
            printf(eg->geometry()->equation().type_name());
            printf("\n");

#endif

            // 如果edge不是只有一个点

            // n_sample1变量记录边离散后点的个数
            // pl1变量记录一组边离散后的顶点
            // 使用api_get_facet_edge_points函数获取边离散后的点
            int n_sample1;
            SPAposition* pl1;
            api_get_facet_edge_points(eg, pl1, n_sample1);

#if FACETER_DEBUG_MODE
            printf("n_sample: %d\n", n_sample1);
#endif

            // 考虑到平面中每个loop是一个首尾相连的环
            // 因此对于每条边, 只需各自加入一个首结点即可
            n_sample1--;

            // 记录一个loop中的点的个数
            n_sample += n_sample1;

            // 使用flag判断是否顺序加入点到points中
            if(flag == 1) {
// 在调试模式下打印当前处理的点的坐标信息
#if FACETER_DEBUG_MODE
                printf("It's same %f %f %f\n", pl1[n_sample1].coordinate(0), pl1[n_sample1].coordinate(1), pl1[n_sample1].coordinate(2));
#endif
                // 顺序加入点
                for(int i = 0; i < n_sample1; i++) {
                    SPAposition pos = pl1[i];
#if FACETER_DEBUG_MODE
                    printf("temp: %f %f %f\n", pos.coordinate(0), pos.coordinate(1), pos.coordinate(2));
#endif
                    // 将每个点的参数坐标加入points变量中

                    points.push_back({f->geometry()->equation().param(pos).u, f->geometry()->equation().param(pos).v});
                }
            } else {
#if FACETER_DEBUG_MODE
                printf("It's not same %f %f %f\n", pl1[n_sample1].coordinate(0), pl1[n_sample1].coordinate(1), pl1[n_sample1].coordinate(2));
#endif
                // 逆序加入点
                for(int i = n_sample1; i > 0; i--) {
                    SPAposition pos = pl1[i];
#if FACETER_DEBUG_MODE
                    printf("temp: %f %f %f\n", pos.coordinate(0), pos.coordinate(1), pos.coordinate(2));
#endif
                    // 将每个点的参数坐标加入points变量中
                    points.push_back({f->geometry()->equation().param(pos).u, f->geometry()->equation().param(pos).v});
                }
            }
        }

        // 对于loop中有大于1个点的情况, 增加边到edges变量中
        // 注意首尾结点的连接
        if(n_sample > 1) {
            // beg是当前环的起始点下标
            for(int i = beg; i < beg + n_sample - 1; i++) {
                edges.push_back(CDT::Edge(i, i + 1));
            }
            // 将环的最后一个点连接回环的第一个点
            edges.push_back(CDT::Edge(beg + n_sample - 1, beg));
        }

        // 下次loop的点的开始下标
        beg += n_sample;
    }

    // 生成CDT三角剖分对象
    CDT::Triangulation<double> cdt;

    // 去除重复点, 重映射边
    CDT::RemoveDuplicatesAndRemapEdges(points, edges);

    // 插入点和边
    cdt.insertVertices(points);
    cdt.insertEdges(edges);

    // 使用带限制的CDT函数进行离散
    cdt.eraseOuterTrianglesAndHoles();
    // 更新points为三角剖分后的顶点集
    points = cdt.vertices;

// 在调试模式下打印三角剖分后的顶点个三角形
#if FACETER_DEBUG_MODE
    for(const auto& v: cdt.vertices) {
        printf("%f %f\n", v.x, v.y);
    }
    for(const auto& t: cdt.triangles) {
        printf("%d, %d, %d\n", t.vertices[0], t.vertices[1], t.vertices[2]);
    }
#endif

    // 将points变量中的点加入mesh中
    // 需要判断点的法向和面的法向是否相同。若不同，需要对点的方向取逆向后加入
    // "gme"是网格名称
    INDEXED_MESH* mesh = new INDEXED_MESH( points.size(), cdt.triangles.size(), cdt.triangles.size() * 3);
    // 遍历所有点
    for(int i = 0; i < points.size(); i++) {
        // 参数空间中点的位置
        SPApar_pos par_pos(points[i].x, points[i].y);
        // 将参数位置转换为3D位置
        SPAposition pos = f->geometry()->equation().eval_position(par_pos);
        SPAunit_vector uv = f->geometry()->equation().eval_normal(par_pos);
        if(r) {
            uv = -uv;
        }
        mesh->add_vertex( pos, uv, par_pos);
    }

    // 根据CDT离散后得到的三角形，离散平面
    // m是三角形索引
    int m = 0;
    // 使用循环迭代cdt.triangles中的每个三角形t
    for(const auto& t: cdt.triangles) {
        // 向mesh对象中添加一个新的多边形
        mesh->add_polygon( m, 3);
        // 获取新添加的多边形对象poly0
        indexed_polygon* poly0 = mesh->get_polygon( m);
        // 将三角形的顶点分别设置为多边形的顶点
        poly0->set_vertex( 0, &mesh->get_vertex( t.vertices[0]));
        poly0->set_vertex( 1, &mesh->get_vertex( t.vertices[1]));
        poly0->set_vertex( 2, &mesh->get_vertex( t.vertices[2]));

// 用于调试
#if FACETER_DEBUG_MODE
        printf("%d %d %d\n", t.vertices[0], t.vertices[1], t.vertices[2]);
#endif
        m++;
    }

    // 将生成的三角网格与平面进行关联
    attach_indexed_mesh_to_face(f, mesh);

    return outcome();
}
