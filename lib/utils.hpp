#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef __KEYWAIT__
#define __KEYWAIT__                                               \
    {                                                             \
        std::cout << "Press any key to continue..." << std::endl; \
        char i;                                                   \
        std::cin >> i;                                            \
    }
#endif

// テンプレート関数でキー入力待ちを汎用化
template <typename T>
void WaitForKey(const std::string& message) {
    std::cout << message << std::endl;
    T input;
    std::cin >> input;
}

// エンターキー入力まで待機する関数
void WaitForEnter(const std::string& message = "Press Enter to continue...") {
    std::cout << message << std::endl;
    while (std::cin.get() != '\n') continue;
}

struct Point3D {
    double x, y, z;
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

template <typename T>
std::vector<T> createSphereMesh(const double ROI_radius,
                                const int divisions,
                                const T& center) {
    std::vector<T> meshPoints;
    if (divisions <= 0) return meshPoints;

    double step = (2.0 * ROI_radius) / (divisions - 1);
    double radiusSquared = ROI_radius * ROI_radius;

    for (int i = 0; i < divisions; ++i) {
        for (int j = 0; j < divisions; ++j) {
            for (int k = 0; k < divisions; ++k) {
                double cx, cy, cz;
                if constexpr (std::is_same<T, Point3D>::value) {
                    cx = center.x;
                    cy = center.y;
                    cz = center.z;
                } else if constexpr (std::is_same<T, std::array<double, 3>>::
                                         value) {
                    cx = center[0];
                    cy = center[1];
                    cz = center[2];
                }
                double x = cx - ROI_radius + i * step;
                double y = cy - ROI_radius + j * step;
                double z = cz - ROI_radius + k * step;

                double distanceSquared = std::pow(x - cx, 2) +
                                         std::pow(y - cy, 2) +
                                         std::pow(z - cz, 2);

                if (distanceSquared <= radiusSquared) {
                    meshPoints.emplace_back(x, y, z);
                }
            }
        }
    }

    std::sort(meshPoints.begin(), meshPoints.end(), [](const T& a, const T& b) {
        if constexpr (std::is_same<T, Point3D>::value) {
            if (a.z != b.z) return a.z < b.z;
            if (a.y != b.y) return a.y < b.y;
            return a.x < b.x;
        } else if constexpr (std::is_same<T, std::array<double, 3>>::value) {
            if (a[2] != b[2]) return a[2] < b[2];
            if (a[1] != b[1]) return a[1] < b[1];
            return a[0] < b[0];
        }
    });

    return meshPoints;
}

template <typename T>
std::vector<T> create_cube_mesh(const double ROI_length,
                                const int divisions,
                                const T& center) {
    std::vector<T> meshPoints;
    if (divisions <= 0) return meshPoints;

    double step = (ROI_length) / (divisions - 1);

    for (int i = 0; i < divisions; ++i) {
        for (int j = 0; j < divisions; ++j) {
            for (int k = 0; k < divisions; ++k) {
                if constexpr (std::is_same<T, Point3D>::value) {
                    double x = center.x - ROI_length + i * step;
                    double y = center.y - ROI_length + j * step;
                    double z = center.z - ROI_length + k * step;
                    meshPoints.emplace_back(x, y, z);
                } else if constexpr (std::is_same<T, std::array<double, 3>>::
                                         value) {
                    double x = center[0] - ROI_length + i * step;
                    double y = center[1] - ROI_length + j * step;
                    double z = center[2] - ROI_length + k * step;
                    meshPoints.emplace_back(x, y, z);
                }
            }
        }
    }

    std::sort(meshPoints.begin(), meshPoints.end(), [](const T& a, const T& b) {
        if constexpr (std::is_same<T, Point3D>::value) {
            if (a.z != b.z) return a.z < b.z;
            if (a.y != b.y) return a.y < b.y;
            return a.x < b.x;
        } else if constexpr (std::is_same<T, std::array<double, 3>>::value) {
            if (a[2] != b[2]) return a[2] < b[2];
            if (a[1] != b[1]) return a[1] < b[1];
            return a[0] < b[0];
        }
    });

    return meshPoints;
}

std::vector<Point3D> create_square_mesh(const double ROI_length,
                                        const int divisions,
                                        const Point3D& center = Point3D(0,
                                                                        0,
                                                                        0),
                                        const double z = 0) {
    std::vector<Point3D> meshPoints;
    if (divisions <= 0) return meshPoints;

    double step = (2 * ROI_length) / (divisions - 1);

    for (int i = 0; i < divisions; ++i) {
        for (int j = 0; j < divisions; ++j) {
            for (int k = 0; k < divisions; ++k) {
                double x = center.x - ROI_length + i * step;
                double y = center.y - ROI_length + j * step;
                meshPoints.emplace_back(x, y, z);
            }
        }
    }

    std::sort(meshPoints.begin(),
              meshPoints.end(),
              [](const Point3D& a, const Point3D& b) {
                  if (a.z != b.z) return a.z < b.z;
                  if (a.y != b.y) return a.y < b.y;
                  return a.x < b.x;
              });

    return meshPoints;
}

void displayProgressBar(double progress, int barWidth = 40) {
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

std::string getcurrent_date() {
    // 現在の時刻を取得
    std::time_t now = std::time(nullptr);
    std::tm* ltm = std::localtime(&now);

    // 日付時刻をフォーマット
    std::ostringstream oss;
    oss << std::put_time(ltm, "%y_%m%d_%H%M");
    return oss.str();
}

#endif  // UTILS_HPP