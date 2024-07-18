#include <iostream>
#include <fstream>
#include <vector>

void writeDataToFile(
    const std::vector<std::pair<double, double>>& points,
    const std::vector<int>& boundary_indices,
    const std::vector<std::pair<int, int>>& triangulation_edges,
    const std::vector<std::pair<int, int>>& constraint_edges
) {
    std::ofstream outfile("data.txt");

    // 输出二维点集数组
    outfile << "Points:\n";
    for (const auto& p : points) {
        outfile << p.first << " " << p.second << "\n";
    }

    // 输出边界点索引数组
    outfile << "Boundary Indices:\n";
    for (int idx : boundary_indices) {
        outfile << idx << "\n";
    }

    // 输出三角化后形成的边集合
    outfile << "Triangulation Edges:\n";
    for (const auto& e : triangulation_edges) {
        outfile << e.first << " " << e.second << "\n";
    }

    // 输出限制边集合
    outfile << "Constraint Edges:\n";
    for (const auto& e : constraint_edges) {
        outfile << e.first << " " << e.second << "\n";
    }

    outfile.close();
}

int main() {
    std::vector<std::pair<double, double>> points = {
        {0.0, 0.0}, {1.0, 0.0}, {0.5, 0.5}, {0.5, 1.0}
    };

    std::vector<int> boundary_indices = {0, 1, 3};

    std::vector<std::pair<int, int>> triangulation_edges = {
        {0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 0}
    };

    std::vector<std::pair<int, int>> constraint_edges = {
        {0, 3}, {1, 3}
    };

    writeDataToFile(points, boundary_indices, triangulation_edges, constraint_edges);

    return 0;
}