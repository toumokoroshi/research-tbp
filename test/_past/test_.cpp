#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

// #include <boost\numeric\odeint.hpp>
struct Point3D {
    double x, y, z;
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

std::vector<Point3D> create_PointsonAxis(const double ROI_radius,
                                      const int divisions,
                                      const Point3D& center) {
    std::vector<Point3D> meshPoints;
    if (divisions <= 0) return meshPoints;

    double step = (2.0 * ROI_radius) / (divisions - 1);
    double radiusSquared = ROI_radius * ROI_radius;

    for (int i = 0; i < divisions; ++i) {
        double x = center.x - ROI_radius + i * step;
        double y = center.y - ROI_radius + i * step;
        double z = center.z - ROI_radius + i * step;

        double x_ = std::sqrt(radiusSquared - y * y);
        double y_ = std::sqrt(radiusSquared - z * z);
        double z_ = std::sqrt(radiusSquared - x * x);

        meshPoints.emplace_back(x, 0.0, 0.0);
        meshPoints.emplace_back(0.0, y, 0.0);
        meshPoints.emplace_back(0.0, 0.0, z);

        meshPoints.emplace_back(x_, y, 0.0);
        meshPoints.emplace_back(x, 0.0, z_);
        meshPoints.emplace_back(0.0, y_, z);

        meshPoints.emplace_back(-x_, y, 0.0);
        meshPoints.emplace_back(x, 0.0, -z_);
        meshPoints.emplace_back(0.0, -y_, z);
    }

    std::sort(meshPoints.begin(), meshPoints.end(), [](const Point3D& a, const Point3D& b) {
        if (a.z != b.z) return a.z < b.z;
        if (a.y != b.y) return a.y < b.y;
        return a.x < b.x;
    });

    return meshPoints;
}


int main() {
    std::vector<Point3D> meshPoints = create_PointsonAxis(0.01, 10, Point3D(0.0, 0.0, 0.0));

    std::ofstream out("sphere_mesh.dat");
    for (const auto& point : meshPoints) {
        out << point.x << '\t' << point.y << '\t' << point.z << '\n';
    }

    return 0;
}