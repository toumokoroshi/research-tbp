/**
 * @file main_sali_3d copy.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-01-28
 * 
 * @note rtbp
 * 
 */

#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <utils.hpp>
#include <vector3d.hpp>
#include <rtbp.hpp>

#define _DEBUG

Vector3d calc_velocity(const Point3D& point,
                       const double v_abs,
                       const double mu,
                       const double inclination,
                       const double OMEGA,
                       const double theta = 0.0);

std::vector<std::streampos> indexFile(const std::string& filename);

std::string readSpecificLine(const std::string& filename,
                             const std::vector<std::streampos>& linePositions,
                             int targetLine);

int main() {

    using namespace rtbp::crtbp;
    constexpr double MU = 3.003e-6;  // 地球と太陽の質量比
    constexpr double HEADER_SIZE = 9;
    Point3D MeshCenter(1.0 - MU, 0, 0);
    double ROI_length = 0;

    std::cout << "<>----------------------------------------------------------------" << std::endl;
    std::cout << "<>            CRTBP 3dSALI Calculation based on Jacobi Integral" << std::endl;
    std::cout << "<>----------------------------------------------------------------\n\n\n\n"
              << std::endl;
    // #ifndef _DEBUG
    std::cout << "<>****************************************************************" << std::endl;
    std::cout << "<>  [mode] : " << std::endl;
    std::cout << "<>        1. New simulation" << std::endl;
    std::cout << "<>        2. Detailed simulation for existing data" << std::endl;
    std::cout << "<>        else. Exit" << std::endl;
    std::cout << "<>      Input number : ";
    char mode;
    std::cin >> mode;
    std::cout << std::endl;
    if (mode == '1') {
        std::cout << "selected mode : New simulation\n" << std::endl;
    } else if (mode == '2') {
        std::cout << "selected mode : Detailed simulation for existing data\n" << std::endl;
    } else {
        std::cout << "selected mode : Exit\n" << std::endl;
        return 0;
    }

    if (mode == '2') {
        // ファイルを読み込んで、ターゲットのメッシュ番号を指定
        std::cout << "<>        [input file name to refer] : ";
        std::string filename_interested;
        std::cin >> filename_interested;
        std::cout << std::endl;
        // メッシュ番号を指定
        std::cout << "<>        [input mesh number to focus on] : ";
        std::string mesh_num_of_interest;
        std::cin >> mesh_num_of_interest;
        std::cout << std::endl;

        std::cout << "<>        [input the length of ROI] : ";
        std::string ROI_length_;
        std::cin >> ROI_length_;
        std::cout << std::endl;
        ROI_length = std::stod(ROI_length_);

        // ファイル読み込み
        std::ifstream ifs(filename_interested);
        if (!ifs) {
            std::cerr << "Can't open file : " << filename_interested << std::endl;
            return -1;
        }
        std::vector<std::streampos> linePositions = indexFile(filename_interested);

        // 指定した行を読み込む
        int targetLine = std::stoi(mesh_num_of_interest) + HEADER_SIZE;  // 読み込みたい行番号
        std::string line = readSpecificLine(filename_interested, linePositions, targetLine);
        std::cout << "<>        interested line : " << line << std::endl;
        std::stringstream ss(line);
        std::array<double, 7> data;
        for (int i = 0; i < 7; i++) {
            ss >> data[i];
        }

        MeshCenter.x = data[2];
        MeshCenter.y = data[3];
        MeshCenter.z = data[4];
    }
    // #endif
    std::cout << "<>****************************************************************" << std::endl;
    std::cout << "<>           --simulation config--\n" << std::endl;
    std::cout << "<>    reading config file\n" << std::endl;

    char mode2;
    std::cout << "<>  [single simulation or continuous simulation] : " << std::endl;
    std::cout << "<>        1. Scan Earth's vecinity (single config file)" << std::endl;
    std::cout << "<>        2. Scan Earth's vecinity (multi-config file)" << std::endl;
    std::cout << "<>        3. Scan the optional trajectory (single config file)" << std::endl;
    std::cout << "<>        else. Exit" << std::endl;
    std::cout << "<>      Input number : ";
    std::cin >> mode2;
    std::cout << std::endl;

    double is_continuous = 0;
    if (mode2 == '1') {
        std::cout << "<>    selected mode : single simulation\n" << std::endl;
        is_continuous = 0;
    } else if (mode2 == '2') {
        std::cout << "<>    selected mode : continuous simulation\n" << std::endl;
        is_continuous = 1;
    } else if (mode2 == '3')
        std::cout << "<>    selected mode : Scan the optional trajectory\n" << std::endl;
    else {
        std::cout << "selected mode : Exit\n" << std::endl;
        return 0;
    }
    

    // 設定ファイル読み込み
    int MESH_SIZE = 0;
    double CALC_TIMESTEP = 0;
    double SALI_CALCTIME_THRESHOLD = 0;
    double SOI_RADIUS = 0;
    double FOREBIDDEN_AREA_RADIUS = 0;
    double JACOBI_INTEGRAL = 0;
    double inclination = 0;
    double OMEGA = 0;
    double THETA = 0;
    std::ifstream ifs;
    std::string configfilename;
    if (mode2 == '1') {
        configfilename = "configdata/3D_crtbp_SALI/3DSALIconfig.txt";
    } else if (mode2 == '2') {
        configfilename = "configdata/3D_crtbp_SALI/3DSALIconfig_1.txt";
    }

    ifs.open(configfilename);

    if (!ifs) {
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }

    int configdata_num = 1;
    //  実行時間の計測
    auto start_ofall = std::chrono::system_clock::now();
    while (ifs) {
        double progress = 0;
        auto start = std::chrono::system_clock::now();
        std::string str;
        std::cout << std::setprecision(10);
        while (std::getline(ifs, str)) {
            if (str.find("MESH SIZE") != std::string::npos) {
                MESH_SIZE = std::stoi(str.substr(str.find("=") + 1));
                std::cout << "<>        MESH SIZE : " << MESH_SIZE << std::endl;
            } else if (str.find("CALC TIMESTEP") != std::string::npos) {
                CALC_TIMESTEP = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        CALC TIMESTEP : " << CALC_TIMESTEP << std::endl;
            } else if (str.find("SALI CALCTIME THRESHOLD") != std::string::npos) {
                SALI_CALCTIME_THRESHOLD = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        SALI CALCTIME THRESHOLD : " << SALI_CALCTIME_THRESHOLD
                          << std::endl;
            } else if (str.find("RADIUS OF SOI") != std::string::npos) {
                SOI_RADIUS = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        SOI RADIUS : " << SOI_RADIUS << std::endl;
            } else if (str.find("RADIUS OF FOREBIDDEN AREA") != std::string::npos) {
                FOREBIDDEN_AREA_RADIUS = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        RADIUS OF FOREBIDDEN AREA : " << FOREBIDDEN_AREA_RADIUS
                          << std::endl;
            } else if (str.find("JACOBI INTEGRAL") != std::string::npos) {
                JACOBI_INTEGRAL = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        JACOBI INTEGRAL : " << JACOBI_INTEGRAL << std::endl;
            } else if (str.find("INCLINATION AGAINST XY PLANE(deg)") != std::string::npos) {
                inclination = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        INCLINATION(deg) : " << inclination << std::endl;
                inclination = inclination * std::acos(-1) / 180.;  // deg to rad
            } else if (str.find("LONGTITUDE AGAINST X AXIS+(deg)") != std::string::npos) {
                OMEGA = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        LONGTITUDE(deg) : " << OMEGA << std::endl;
                OMEGA = OMEGA * std::acos(-1) / 180.;  // deg to rad
            } else if (str.find("DEGREE FROM TANGENT") != std::string::npos) {
                THETA = std::stod(str.substr(str.find("=") + 1));
                std::cout << "<>        DEGREE FROM TANGENT(deg) : " << THETA << std::endl;
                THETA = THETA * std::acos(-1) / 180.;  // deg to rad
            }
        }
        ifs.close();
        std::cout << "<>    config file read successfully\n" << std::endl;
        if (is_continuous == 0) {
            std::cout << "<>  [Check the param and Press any key to go ahead]" << std::endl;
            __KEYWAIT__
        }

        std::cout << std::endl;
        std::cout << "<>           --SALI caluculation --" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;

        std::cout << "<>        Generating mesh ";
        std::vector<std::array<double,3>> meshPoints;
        if (mode == '1') {
            std::cout << "based on SOI radius" << std::endl;
            meshPoints = createSphereMesh(SOI_RADIUS, MESH_SIZE, MeshCenter);
        } else if (mode == '2') {
            std::cout << "based on the specified point" << std::endl;
            meshPoints = create_cube_mesh(ROI_length, MESH_SIZE, MeshCenter);
        }

        int countt = 0;
        countt = meshPoints.size();

        std::cout << std::endl;
        std::cout << "<>        " << countt << " mesh generated successfully" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;

        std::cout << "<>        Start calclation" << std::endl;

        /* メッシュ上の点を用いてSALIを計算 */
        constexpr double perturbation = 1e-10;

        std::vector<std::array<double, 10>> SALI_data;

        // [0]:count, [1]time , [2]:x, [3]:y, [4]:z,
        // [5]:jacobi constant
        // [6]:SALI(x-biased&y-biased), [7]:SALI(x-biased&z-biased),
        // [8]:SALI(y-biased&z-biased)
        SALI_data.reserve(countt);

        int totalIterations = meshPoints.size();
        for (const auto& point : meshPoints) {
            displayProgressBar(progress);
            // 0:計算終了, 1:計算継続
            bool calc_traj = 1;
            // 0:計算継続, 1:計算中断,途中でhill球から抜け出したり、禁止領域に入った場合
            bool abort_calc = 0;
            // 速度が定義できない場合（速度の絶対値が負になる場合）
            bool velo_err = 0;

            static int mesh_num = 1;  // メッシュ番号

            // non-biased velocity
            double v_abs = calc_v_abs(point, MU, JACOBI_INTEGRAL);
            double vx, vy, vz = 0;

            if (v_abs > 0) {
                Vector3d velocity = calc_velocity(point, v_abs, MU, inclination, OMEGA, THETA);
                vx = velocity.x();
                vy = velocity.y();
                vz = velocity.z();
            } else {
                calc_traj = 0;
                velo_err = 1;
            }


            double time = 0.0;

            double SALIxy = -1.0;
            double SALIyz = -1.0;
            double SALIzx = -1.0;

            if (!velo_err) {
                while (calc_traj) {
     
                    if (ref_state.calc_r2() > SOI_RADIUS ||
                        ref_state.calc_r2() < FOREBIDDEN_AREA_RADIUS) {
                        abort_calc = 1;
                        calc_traj = 0;
                        continue;
                    }
                    ref_state.RK4_step_noncanonical();
                    perturbed_state1.RK4_step_noncanonical();
                    perturbed_state2.RK4_step_noncanonical();
#ifndef SALI_only_XY
                    perturbed_state3.RK4_step_noncanonical();
#endif
                    time += CALC_TIMESTEP;
                }
                if (!abort_calc) {
#ifdef SALI_only_XY
                    SALI = calc_SALI(ref_state.current_state(),
                                     perturbed_state1.current_state(),
                                     perturbed_state2.current_state());
#endif
#ifndef SALI_only_XY
                    SALIxy = calc_SALI(ref_state.current_state(),
                                       perturbed_state1.current_state(),
                                       perturbed_state2.current_state());
                    SALIyz = calc_SALI(ref_state.current_state(),
                                       perturbed_state2.current_state(),
                                       perturbed_state3.current_state());
                    SALIzx = calc_SALI(ref_state.current_state(),
                                       perturbed_state1.current_state(),
                                       perturbed_state3.current_state());
#endif
                }
            } else {
#ifdef SALI_only_XY
                SALI = -1.0;
#endif
#ifndef SALI_only_XY
                SALIxy = -1.0;
                SALIyz = -1.0;
                SALIzx = -1.0;
#endif
            }
            double jacobiii = calc_jacobi_integral(ref_state.current_state(), MU);
#ifdef SALI_only_XY
            std::array<double, 7> SALI_data1{
                static_cast<double>(mesh_num), time, point.x, point.y, point.z, jacobiii, SALI};
#endif
#ifndef SALI_only_XY
            std::array<double, 9> SALI_data1{static_cast<double>(mesh_num),
                                             time,
                                             point.x,
                                             point.y,
                                             point.z,
                                             jacobiii,
                                             SALIxy,
                                             SALIyz,
                                             SALIzx};
#endif
            SALI_data.emplace_back(SALI_data1);
            // 進捗を計算
            progress = (mesh_num + 1.0) / totalIterations;

            // プログレスバーを表示
            displayProgressBar(progress);

            mesh_num++;
        }

        //  日付と時刻をファイル名にしてファイル出力
        std::string filename = "results/SALI/3DSALI_" + getcurrent_date() + ".dat";
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
        ofs1 << "INITIAL JACOBI INTEGRAL=" << JACOBI_INTEGRAL << std::endl;
        ofs1 << "INCLINATION AGAINST XY PLANE=" << inclination / std::acos(-1) * 180.
             << std::endl;  // rad to deg
        ofs1 << "LONGTITUDE AGAINST X AXIS=" << OMEGA / std::acos(-1) * 180.
             << std::endl;  // rad to
                            // deg
        ofs1 << "DEGREE FROM TANGENT(deg)=" << THETA / std::acos(-1) * 180. << std::endl;
#ifdef SALI_only_XY
        ofs1 << "mesh_num, time, x, y, z, jacobi constant, SALI" << std::endl;
#endif
#ifndef SALI_only_XY
        ofs1 << "mesh_num, time, x, y, z, jacobi constant, SALIxy,SALIyz,SALIxz," << std::endl;
#endif

        for (const auto& data : SALI_data) {
#ifdef SALI_only_XY
            ofs1 << std::setprecision(0) << std::fixed << data[0] << " " << std::setprecision(4)
                 << data[1] << " " << std::setprecision(15) << std::fixed << data[2] << " "
                 << data[3] << " " << data[4] << " " << data[5] << " " << data[6] << " "
                 << std::endl;
#endif
#ifndef SALI_only_XY
            ofs1 << std::setprecision(0) << std::fixed << data[0] << " " << std::setprecision(4)
                 << data[1] << " " << std::setprecision(15) << std::fixed << data[2] << " "
                 << data[3] << " " << data[4] << " " << data[5] << " " << data[6] << " " << data[7]
                 << " " << data[8] << std::endl;
#endif
        }
        ofs1.close();
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
        if (mode2 == '1') {
            break;
        }
        std::cout << "<>        elapsed time : " << hour << "h " << min << "m " << sec << "s "
                  << msec << "ms" << std::endl;
        std::cout << "<>        Simulation for " << configfilename << " finished" << std::endl;
        configdata_num++;
        configfilename =
            "configdata/3D_crtbp_SALI/3DSALIconfig_" + std::to_string(configdata_num) + ".txt";
        std::cout << "<>        Next config file : " << configfilename << std::endl;
        ifs.open(configfilename);
    }


    // 実行時間の計測
    auto end = std::chrono::system_clock::now();

    // 実行時間の計算
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_ofall);

    // 時、分、秒、ミリ秒に分解
    auto msec = duration.count() % 1000;
    auto sec = duration.count() / 1000 % 60;
    auto min = duration.count() / 1000 / 60 % 60;
    auto hour = duration.count() / 1000 / 60 / 60;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "<>        Calculation finished" << std::endl;
    std::cout << "<>        Total elapsed time : " << hour << "h " << min << "m " << sec << "s "
              << msec << "ms" << std::endl;

    return 0;
}

std::vector<std::streampos> indexFile(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("ファイルを開けませんでした: " + filename);
    }

    std::vector<std::streampos> linePositions;
    std::string line;

    // 最初の行の開始位置を記録
    linePositions.push_back(file.tellg());

    // 各行の開始位置を記録
    while (std::getline(file, line)) {
        linePositions.push_back(file.tellg());
    }

    return linePositions;
}

std::string readSpecificLine(const std::string& filename,
                             const std::vector<std::streampos>& linePositions,
                             int targetLine) {
    if (targetLine < 1 || targetLine >= static_cast<int>(linePositions.size())) {
        throw std::out_of_range("指定した行が範囲外です: " + std::to_string(targetLine));
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("ファイルを開けませんでした: " + filename);
    }

    // 指定した行の位置にシーク
    file.seekg(linePositions[targetLine - 1]);

    std::string line;
    std::getline(file, line);
    return line;
}