#include <fstream>
#include <iostream>
#include <string>
#include <utils.hpp>
#include <vector3d.hpp>
#include <vector>
#include <cmath>
#include <omp.h>

using std::cout;
using std::endl;

double potential(double x, double y, double z, double mu) {
    auto calc_r1 = [&x, &y, &z, &mu]() -> double {
        return std::sqrt(std::pow(x + mu, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
    };
    auto calc_r2 = [&x, &y, &z, &mu]() -> double {
        return std::sqrt(std::pow(x - 1. + mu, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
    };

    double r1 = calc_r1();
    double r2 = calc_r2();

    return 0.5 * (x * x + y * y) + (1. - mu) / r1 + mu / r2;
}

int main(void) {
    constexpr double mu = 3.003e-6;
    constexpr double eps = 1e-6;
    double C_j_start, C_j_end, C_j_step;

    cout << "Enter the start value of Jacobi constant: ";
    std::cin >> C_j_start;
    cout << endl;
    cout << "Enter the end value of Jacobi constant: ";
    std::cin >> C_j_end;
    cout << endl;
    cout << "Enter the step value of Jacobi constant: ";
    std::cin >> C_j_step;
    cout << endl;

    std::vector<std::array<double, 4>> zvcPoints;
    std::vector<Point3D> meshPoints = create_square_mesh(0.05, 200, Point3D(1 - mu, 0, 0));

    int num_steps = static_cast<int>((C_j_end - C_j_start) / C_j_step) + 1;

    #pragma omp parallel
    {
        std::vector<std::array<double, 4>> local_zvcPoints;

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < num_steps; ++i) {
            double C_j = C_j_start + i * C_j_step;
            for (const auto& point : meshPoints) {
                double zvp = 2 * potential(point.x, point.y, point.z, mu) - C_j;
                if (std::abs(zvp) < eps) {
                    local_zvcPoints.push_back({point.x, point.y, point.z, C_j});
                }
            }
        }

        #pragma omp critical
        {
            zvcPoints.insert(zvcPoints.end(), local_zvcPoints.begin(), local_zvcPoints.end());
        }
    }

    std::ofstream ofs("results/zvc.dat");

    if (!ofs) {
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }
    
    for (const auto& point : zvcPoints) {
        ofs << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
    }

    return 0;
}