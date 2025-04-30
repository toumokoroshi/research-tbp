////////////////////////////////////////////////////////////////////////////////
/// @file           .cpp
/// @brief          受け取った小惑星のデータを回転座標系に変換し、SALIを計算したい
/// @author         tabata
/// @date           2024/07/25
/// @par            edittor/ date/ version/ description
///                 tabata/ 2024/7/25/ 0.0/ 初版作成
///
////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "chaosmap_generator.hpp"
#include "crtbp.hpp"
#include "vector3d.hpp"

#define __KEYWAIT__    \
    {                  \
        std::cout << "Press any key to continue..." << std::endl; \
        char i;        \
        std::cin >> i; \
    }


std::vector<double> readdata_fromfile(const std::string& filename);
std::string getcurrent_date();

int main() {
    constexpr double MU = 3.003e-6;
    constexpr double SoI_radiusUS = 0.01;
    constexpr double FOREBIDDEN_AREA_RADIUS = 0.00007;
    constexpr double TIMESTEP = 0.001;
    constexpr double SALI_CALCTIME_THRESHOLD = 10;

    //////////ファイルからデータを読み込む
    std::string datafilename = "configdata/DE441_crtbp/input/asteroid.dat";

    std::vector<double> red_databuf =
        readdata_fromfile(datafilename);  // スペース区切りのデータを読み込む
    if (red_databuf.empty()) {
        std::cerr << "Can't open file : " << datafilename << std::endl;
        return -1;
    }
    std::cout << "<>               Data read successfully" << std::endl;

    Vector3d init_ast_G_pos(red_databuf[0], red_databuf[1], red_databuf[2]);
    Vector3d init_ast_G_vel(red_databuf[3], red_databuf[4], red_databuf[5]);
    std::array<Vector3d, 2> init_ast_G{init_ast_G_pos, init_ast_G_vel};
    // check
    std::cout << "<> Initial asteroid state -> " << std::endl;
    std::cout << "<>               Position : " << init_ast_G[0].x() << " "
              << init_ast_G[0].y() << " " << init_ast_G[0].z() << std::endl;
    std::cout << "<>               Velocity : " << init_ast_G[1].x() << " "
              << init_ast_G[1].y() << " " << init_ast_G[1].z() << std::endl;

    Vector3d init_e_G_pos(red_databuf[6], red_databuf[7], red_databuf[8]);
    Vector3d init_e_G_vel(red_databuf[9], red_databuf[10], red_databuf[11]);
    std::array<Vector3d, 2> init_e_G_state{init_e_G_pos, init_e_G_vel};

    std::cout << "<>                    Data read successfully " << std::endl;
    std::cout << std::endl;

    __KEYWAIT__
    std::cout << "<>           ------------------------------------- " << std::endl;
    std::cout << "<>                        Simulation Start " << std::endl;
    std::cout << "<>           ------------------------------------- " << std::endl;
    std::cout << "<>" << std::endl;
    std::cout << "<>           --trajectory propagation start--" << std::endl;
    std::cout << std::endl;

    CRTBPChaosMap chaos_map(
        init_e_G_state, init_ast_G, SoI_radiusUS, FOREBIDDEN_AREA_RADIUS, MU);
    // まずは軌道計算

    chaos_map.frame_transformation_();

    return 0;
}

std::vector<double> readdata_fromfile(const std::string& filename) {
    std::vector<double> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return data;  // 空のvectorを返す
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double value;
        while (iss >> value) {
            data.push_back(value);
        }
    }

    file.close();
    return data;
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
