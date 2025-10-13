/// @file           main_sali.cpp
/// @brief DE441を用いず、小惑星の速度を仮定し、ヤコビ積分を基準にCRTBPにおいてSALIを計算する
/// @author         tabata
/// @date           2024/10/22
/// @par            edittor/      date/ version/ description
///                 tabata/ 2024/11/20/     1.0/ 初版作成
///
////////////////////////////////////////////////////////////////////////////////

#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "crtbp.hpp"

#define __KEYWAIT__    \
    {                  \
        int i;         \
        std::cin >> i; \
    }

#define _DEBUG

struct Point3D {
    double x, y, z;
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

std::vector<Point3D> createSphereMesh(const double radius,
                                      const int divisions,
                                      const Point3D& center = Point3D(0, 0, 0)) {
    std::vector<Point3D> meshPoints;
    if (divisions <= 0) return meshPoints;

    double step = (2.0 * radius) / (divisions - 1);
    double radiusSquared = radius * radius;

    for (int i = 0; i < divisions; ++i) {
        for (int j = 0; j < divisions; ++j) {
            double x = center.x - radius + i * step;
            double y = center.y - radius + j * step;
            double z = center.z;

            double distanceSquared =
                std::pow(x - center.x, 2) + std::pow(y - center.y, 2) + std::pow(z - center.z, 2);

            if (distanceSquared <= radiusSquared) {
                meshPoints.emplace_back(x, y, z);
            }
        }
    }

    std::sort(meshPoints.begin(), meshPoints.end(), [](const Point3D& a, const Point3D& b) {
        if (a.y != b.y) return a.y < b.y;
        return a.x < b.x;
    });

    return meshPoints;
}

// 任意の2組の軌道要素からSALIを計算する関数
double calc_SALI_from2pair(const std::array<double, 6>& ref_state,
                           const std::array<double, 6>& perturbed_state1,
                           const std::array<double, 6>& perturbed_state2) {
    std::array<double, 3> q_ref{ref_state[0], ref_state[1], ref_state[2]};
    std::array<double, 3> q_pertubed1{
        perturbed_state1[0], perturbed_state1[1], perturbed_state1[2]};
    std::array<double, 3> q_pertubed2{
        perturbed_state2[0], perturbed_state2[1], perturbed_state2[2]};
    std::array<double, 3> p_ref{
        ref_state[3] - ref_state[0], ref_state[4] + ref_state[1], ref_state[5]};
    std::array<double, 3> p_pertubed1{perturbed_state1[3] - perturbed_state1[0],
                                      perturbed_state1[4] + perturbed_state1[1],
                                      perturbed_state1[5]};
    std::array<double, 3> p_pertubed2{perturbed_state2[3] - perturbed_state2[0],
                                      perturbed_state2[4] + perturbed_state2[1],
                                      perturbed_state2[5]};
    Vector3d deviation_vec1{perturbed_state1[0] - ref_state[0],
                            perturbed_state1[1] - ref_state[1],
                            perturbed_state1[2] - ref_state[2]};
    Vector3d deviation_vec2{perturbed_state2[0] - ref_state[0],
                            perturbed_state2[1] - ref_state[1],
                            perturbed_state2[2] - ref_state[2]};

    Vector3d normalized_dev_vec1 = deviation_vec1.normalise();
    Vector3d normalized_dev_vec2 = deviation_vec2.normalise();

    Vector3d sa = normalized_dev_vec1 - normalized_dev_vec2;
    Vector3d wa = normalized_dev_vec1 + normalized_dev_vec2;

    double sa_norm = sa.magnitude();
    double wa_norm = wa.magnitude();

    // SALIの計算
    return std::min(sa_norm, wa_norm);
}

double calc_r1(const Point3D& point, const double mu) {
    return std::sqrt(std::pow(point.x + mu, 2.) + std::pow(point.y, 2.) + std::pow(point.z, 2.));
}

double calc_r2(const Point3D& point, const double mu) {
    return std::sqrt(std::pow(point.x - 1. + mu, 2.) + std::pow(point.y, 2.) +
                     std::pow(point.z, 2.));
}

double calc_v_abs(const Point3D& point, const double mu, const double JACOBI_INTEGRAL) {
    double r1 = calc_r1(point, mu);
    double r2 = calc_r2(point, mu);
    // std::cout << "uiuoij" <<point.x * point.x + point.y * point.y + 2. * (1. - mu) / r1 + 2. * mu
    // / r2 +
    //                  mu * (1 - mu) - JACOBI_INTEGRAL << std::endl;
    return std::sqrt(point.x * point.x + point.y * point.y + 2. * (1. - mu) / r1 + 2. * mu / r2 +
                     mu * (1. - mu) - JACOBI_INTEGRAL);
}

double calc_jacobi_integral(const std::array<double, 6>& state, const double mu) {
    const double x = state[0];
    const double y = state[1];
    const double z = state[2];
    const double vx = state[3];
    const double vy = state[4];
    const double vz = state[5];

    // 第一質点からの距離
    const double r1 = std::sqrt(std::pow(x + mu, 2) + y * y + z * z);

    // 第二質点からの距離
    const double r2 = std::sqrt(std::pow(x - 1. + mu, 2) + y * y + z * z);

    // ヤコビ積分の計算
    return x * x + y * y - (vx * vx + vy * vy + vz * vz) + 2. * (1. - mu) / r1 + 2. * mu / r2 +
           mu * (1. - mu);
}

void displayProgressBar(double progress, int barWidth = 50) {
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}
int main() {
    //  実行時間の計測
    auto start = std::chrono::system_clock::now();

    constexpr double MU = 3.003e-6;  // 地球と太陽の質量比
    Point3D MeshCenter(1.0 - MU, 0, 0);
    std::cout << "<>****************************************************************" << std::endl;
    std::cout << "<>            CRTBP SALI Calculation based on Jacobi Integral" << std::endl;
    std::cout << "<>****************************************************************\n\n\n\n"
              << std::endl;
#ifndef _DEBUG
    std::cout << "<>****************************************************************" << std::endl;
    std::cout << "<>  [mode] : " << std::endl;
    std::cout << "<>        1. New simulation" << std::endl;
    std::cout << "<>        2. Detailed simulation for existing data" << std::endl;
    std::cout << "<>        else. Exit" << std::endl;
    std::cout << "<>      Input number : ";
    char input1;
    std::cin >> input1;
    std::cout << std::endl;
    if (input1 != '1' && input1 != '2') {
        return 0;
    }
    std::cout << "<>****************************************************************" << std::endl;
    std::cout << "<>  [simulation config]" << std::endl;
    int mesh_num = 0;
    if (input1 == '2') {
        ///////////////////////まだ作ってない
        // 既存データの詳細シミュレーション
        // ファイルを読み込んで、ターゲットのメッシュ番号を指定
        std::cout << "<>        input file name : ";
        std::string filename;
        std::cin >> filename;
        std::cout << std::endl;
        // メッシュ番号を指定
        std::cout << "<>        input mesh number to focus on : ";
        std::cin >> mesh_num;
        std::cout << std::endl;
    }
#endif
    std::cout << "<>    read config file\n\n" << std::endl;

    // 設定ファイル読み込み
    int MESH_SIZE = 0;
    double CALC_TIMESTEP = 0;
    double SALI_CALCTIME_THRESHOLD = 0;
    double SOI_RADIUS = 0;
    double FOREBIDDEN_AREA_RADIUS = 0;
    double DIRECTION = 0;
    double JACOBI_INTEGRAL = 0;
    std::ifstream ifs("configdata/SALI2Dconfig.txt");
    if (!ifs) {
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }
    std::string str;
    std::cout << std::setprecision(10);
    while (std::getline(ifs, str)) {
        if (str.find("MESH_SIZE") != std::string::npos) {
            MESH_SIZE = std::stoi(str.substr(str.find("=") + 1));
            std::cout << "<>        MESH_SIZE : " << MESH_SIZE << std::endl;
        } else if (str.find("CALC_TIMESTEP") != std::string::npos) {
            CALC_TIMESTEP = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        CALC_TIMESTEP : " << CALC_TIMESTEP << std::endl;
        } else if (str.find("SALI_CALCTIME_THRESHOLD") != std::string::npos) {
            SALI_CALCTIME_THRESHOLD = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        SALI_CALCTIME_THRESHOLD : " << SALI_CALCTIME_THRESHOLD
                      << std::endl;
        } else if (str.find("SOI_RADIUS") != std::string::npos) {
            SOI_RADIUS = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        SOI_RADIUS : " << SOI_RADIUS << std::endl;
        } else if (str.find("FOREBIDDEN_AREA_RADIUS") != std::string::npos) {
            FOREBIDDEN_AREA_RADIUS = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        FOREBIDDEN_AREA_RADIUS : " << FOREBIDDEN_AREA_RADIUS
                      << std::endl;
        } else if (str.find("DIRECTION") != std::string::npos) {
            DIRECTION = std::stod(str.substr(str.find("=") + 1));
        } else if (str.find("JACOBI_INTEGRAL") != std::string::npos) {
            JACOBI_INTEGRAL = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        JACOBI_INTEGRAL : " << JACOBI_INTEGRAL << std::endl;
        }
    }
    std::cout << "<>    config file read successfully" << std::endl;
    std::cout << "<>    Press any key to go ahead" << std::endl;
    __KEYWAIT__
#ifndef _DEBUG
    std::cout << "<>****************************************************************" << std::endl;

    std::cout << std::endl;
    std::cout << "<>        Check the param and Press keys to go ahead" << std::endl;
    std::cout << "<>        0: exit" << std::endl;
    std::cout << "<>        else: go ahead" << std::endl;
    char buf;
    std::cin >> buf;
    if (buf == '0') {
        return 0;
    }
#endif

    std::cout << std::endl;
    std::cout << "<>------------------------------------------------" << std::endl;
    std::cout << "<>-          --SALI caluculation --              -" << std::endl;
    std::cout << "<>------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "<>        Generating mesh ";
    // MeshCenter.x = 0;
    // MeshCenter.y = 0;
    // MeshCenter.z = 0;
    std::vector<Point3D> meshPoints = createSphereMesh(SOI_RADIUS, MESH_SIZE, MeshCenter);

    int countt = 0;
    countt = meshPoints.size();

    std::cout << std::endl;
    std::cout << "<>        " << countt << " mesh generated successfully\n" << std::endl;
    std::cout << std::endl;
    __KEYWAIT__

    std::cout << std::endl;

    std::cout << "<>        Start calclation" << std::endl;

    /* メッシュ上の点を用いてSALIを計算 */
    constexpr double perturbation = 1e-10;
    // [0]:count, [1]time , [2]:x, [3]:y, [4]:z, [5]:SALI(x-biased&y-biased),
    // [6]:SALI(x-biased&z-biased), [7]:SALI(y-biased&z-biased)[8]:SALI
    std::vector<std::array<double, 9>> SALI_data;
    SALI_data.reserve(countt);

    int totalIterations = meshPoints.size();
    double progress;
    for (const auto& point : meshPoints) {
        displayProgressBar(progress);
        // 0:計算終了, 1:計算継続
        bool calc_traj = 1;
        // 1:計算中断,途中でhill球から抜け出したり、禁止領域に入った場合
        bool abort_calc = 0;
        // 速度の絶対値が負の場合
        bool velo_err = 0;

        static int mesh_num = 1;  // メッシュ番号

        // non-biased velocity
        double v_abs = calc_v_abs(point, MU, JACOBI_INTEGRAL);
        double v_x, v_y, v_z;
        v_z = 0;
        if (v_abs > 0) {
            Point3D point_projected_to_xy_plane{point.x, point.y, 0};
            v_x = -v_abs * point.y * DIRECTION / calc_r2(point_projected_to_xy_plane, MU);
            v_y =
                v_abs * (point.x - 1. + MU) * DIRECTION / calc_r2(point_projected_to_xy_plane, MU);
        } else {
            calc_traj = 0;
            velo_err = 1;
            v_x = 0;
            v_y = 0;
        }

        // x-biased
        Point3D perturbed_point1(point.x + perturbation, point.y, point.z);
        Point3D point_projected_to_xy_plane{point.x + perturbation, point.y, 0};
        double v_abs_perturbed1 = calc_v_abs(perturbed_point1, MU, JACOBI_INTEGRAL);
        double v_x_perturbed1 = -v_abs_perturbed1 * perturbed_point1.y * DIRECTION /
                                calc_r2(point_projected_to_xy_plane, MU);
        double v_y_perturbed1 = v_abs_perturbed1 * (perturbed_point1.x - 1. + MU) * DIRECTION /
                                calc_r2(point_projected_to_xy_plane, MU);

        // y-biased
        Point3D perturbed_point2(point.x, point.y + perturbation, point.z);
        point_projected_to_xy_plane.x = point.x;
        point_projected_to_xy_plane.x = point.y + perturbation;
        double v_abs_perturbed2 = calc_v_abs(perturbed_point2, MU, JACOBI_INTEGRAL);
        double v_x_perturbed2 = -v_abs_perturbed2 * perturbed_point2.y * DIRECTION /
                                calc_r2(point_projected_to_xy_plane, MU);
        double v_y_perturbed2 = v_abs_perturbed2 * (perturbed_point2.x - 1. + MU) * DIRECTION /
                                calc_r2(point_projected_to_xy_plane, MU);

        // // z-biased
        // Point3D perturbed_point3(point.x, point.y, point.z + perturbation);
        // point_projected_to_xy_plane.x = point.x;
        // point_projected_to_xy_plane.x = point.y;
        // double v_abs_perturbed3 = calc_v_abs(perturbed_point3, MU, JACOBI_INTEGRAL);
        // double v_x_perturbed3 = -v_abs_perturbed3 * perturbed_point3.y * DIRECTION /
        //                         calc_r2(point_projected_to_xy_plane, MU);
        // double v_y_perturbed3 = v_abs_perturbed3 * (perturbed_point3.x - 1. + MU) * DIRECTION /
        //                         calc_r2(point_projected_to_xy_plane, MU);

        // perturbed stateとref stateを用いてCRTBPクラスを生成
        CRTBP ref_state(point.x, point.y, point.z, v_x, v_y, v_z, MU, CALC_TIMESTEP);
        CRTBP perturbed_state1(perturbed_point1.x,
                               perturbed_point1.y,
                               perturbed_point1.z,
                               v_x_perturbed1,
                               v_y_perturbed1,
                               v_z,
                               MU,
                               CALC_TIMESTEP);
        CRTBP perturbed_state2(perturbed_point2.x,
                               perturbed_point2.y,
                               perturbed_point2.z,
                               v_x_perturbed2,
                               v_y_perturbed2,
                               v_z,
                               MU,
                               CALC_TIMESTEP);
        // CRTBP perturbed_state3(perturbed_point3.x,
        //                        perturbed_point3.y,
        //                        perturbed_point3.z,
        //                        v_x_perturbed3,
        //                        v_y_perturbed3,
        //                        v_z,
        //                        MU,
        //                        CALC_TIMESTEP);

        double time = 0.0;
        double SALI1 = -1.0;
        double SALI2 = -1.0;
        double SALI3 = -1.0;
        double SALI = -1.0;
        if (!velo_err) {
            while (calc_traj) {
                if (time > SALI_CALCTIME_THRESHOLD) {
                    calc_traj = 0;
                    continue;
                }
                if (ref_state.calc_r2() > SOI_RADIUS ||
                    ref_state.calc_r2() < FOREBIDDEN_AREA_RADIUS) {
                    abort_calc = 1;
                    calc_traj = 0;
                    continue;
                }
                ref_state.RK4_step_noncanonical();
                perturbed_state1.RK4_step_noncanonical();
                perturbed_state2.RK4_step_noncanonical();
                // perturbed_state3.RK4_step_noncanonical();
                time += CALC_TIMESTEP;
            }
            if (!abort_calc) {
                SALI1 = calc_SALI_from2pair(ref_state.current_state(),
                                            perturbed_state1.current_state(),
                                            perturbed_state2.current_state());
                // SALI2 = calc_SALI_from2pair(ref_state.current_state(),
                //                             perturbed_state1.current_state(),
                //                             perturbed_state3.current_state());
                // SALI3 = calc_SALI_from2pair(ref_state.current_state(),
                //                             perturbed_state2.current_state(),
                //                             perturbed_state3.current_state());
                // SALI = std::min(SALI1, SALI2);
            }
        } else {
            SALI1 = -2.0;
            SALI2 = -2.0;
            SALI3 = -2.0;
            SALI = -2.0;
        }

        std::array<double, 9> SALI_data1{static_cast<double>(mesh_num),
                                         time,
                                         point.x,
                                         point.y,
                                         point.z,
                                         SALI1,
                                         SALI2,
                                         SALI3,
                                         SALI};
        SALI_data.emplace_back(SALI_data1);
        // 進捗を計算
        progress = (mesh_num + 1.0) / totalIterations;

        // プログレスバーを表示
        displayProgressBar(progress);

        mesh_num++;
    }

    //  日付と時刻をファイル名にしてファイル出力
    std::string filename =
        "results/SALI/SALI_2d_" +
        std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + ".txt";
    std::ofstream ofs1(filename);
    if (!ofs1) {
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }
    // 計算のコンフィグレーションをファイルに書き込む
    ofs1 << "MESH SIZE=" << MESH_SIZE << std::endl;
    ofs1 << "CALCULATION TIMESTEP=" << CALC_TIMESTEP << std::endl;
    ofs1 << "SIMULATION TIME=" << SALI_CALCTIME_THRESHOLD << std::endl;
    ofs1 << "RADIUSofSOI=" << SOI_RADIUS << std::endl;
    ofs1 << "FOREBIDDEN AREA RADIUS=" << FOREBIDDEN_AREA_RADIUS << std::endl;
    if (DIRECTION == 1) {
        ofs1 << "PROGRADE" << std::endl;
    } else {
        ofs1 << "RETROGRADE" << std::endl;
    }
    ofs1 << "INITIAL JACOBI INTEGRAL=" << JACOBI_INTEGRAL << std::endl;
    ofs1 << "mesh_num, time, x, y, z, SALI1, SALI2, SALI3, SALI" << std::endl;
    for (const auto& data : SALI_data) {
        ofs1 << std::setprecision(0) << std::fixed << data[0] << " " << std::setprecision(4)
             << data[1] << " " << std::setprecision(15) << std::fixed << data[2] << " " << data[3]
             << " " << data[4] << " " << data[5] << " " << data[6] << " " << data[7] << " "
             << data[8] << std::endl;
    }
    ofs1.close();

    // 実行時間の計測
    auto end = std::chrono::system_clock::now();

    // 実行時間の計算
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // 時、分、秒、ミリ秒に分解
    auto msec = duration.count() % 1000;
    auto sec = duration.count() / 1000 % 60;
    auto min = duration.count() / 1000 / 60 % 60;
    auto hour = duration.count() / 1000 / 60 / 60;
    displayProgressBar(1.0);
    std::cout << std::endl;
    std::cout << "<>        Calculation finished" << std::endl;
    std::cout << "<>        elapsed time : " << hour << "h " << min << "m " << sec << "s " << msec
              << "ms" << std::endl;

    return 0;
}
