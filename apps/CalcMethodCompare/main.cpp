

#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector3d.hpp>
#include <vector>

#include "rtbp.hpp"

/**
 * @brief enumの概要説明
 */
enum class MeshCenter {
  kSUN = 0,
  kEARTH = 1,
};

int main() {
  using namespace param;
  using namespace crtbp;
  using namespace utils;
  namespace fs = std::filesystem;
  // CMakeから渡されたCONFIG_DIRマクロを使用
  const std::string kConfigFilePath = CONFIG_DIR;
  const std::string kCalcConfigPath = kConfigFilePath + "/3D_crtbp_SALI/";
  constexpr std::string kCalcConfigPrefix = "3DSALIconfig_";
  std::string astro_param_file = kConfigFilePath + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);

  const double kAU = astro_params.au;                 // astronomical unit in meters
  const double kGMSUN = astro_params.gm_sun;          // heliocentric gravitational constant m3 s-2
  const double kGMEARTH = astro_params.gm_earth;      // geocentric gravitational constant m3 s-2
  const double kMU = kGMEARTH / (kGMEARTH + kGMSUN);  // mu parameter of Earth-Sun
  std::cout << "-" << std::endl;

  constexpr double HEADER_SIZE = 9;
  State3d<double> MeshCenter{1.0 - kMU, 0, 0};
  double ROI_length = 0;

  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            Integration Method Compareing" << std::endl;
  std::cout << "<>-------------------------------------------------------------"
               "---\n"
            << std::endl;
  std::cout << "<>****************************************************************" << std::endl;

  double is_continuous = 0;
  std::string configfilename;
  std::vector<std::string> config_file_list;

  const std::regex pattern("^" + kCalcConfigPrefix + "\\d+\\.txt$");
  try {
    for (const auto& entry : fs::directory_iterator(kCalcConfigPath)) {
      if (entry.is_regular_file()) {
        std::string filename = entry.path().filename().string();
        if (std::regex_match(filename, pattern)) {
          config_file_list.push_back(fs::absolute(entry.path()).string());
        }
      }
    }
  } catch (fs::filesystem_error& e) {
    std::cerr << "Error accessing directory: " << e.what() << std::endl;
  }
  std::sort(config_file_list.begin(), config_file_list.end(),
            [](const std::string& a, const std::string& b) {
              auto getNumber = [](const std::string& path_str) -> int {
                std::string stem = std::filesystem::path(path_str).stem().string();
                size_t lastUnderscore = stem.find_last_of('_');
                return std::stoi(stem.substr(lastUnderscore + 1));
              };
              return getNumber(a) < getNumber(b);
            });
  std::cout << "<> Loaded config file list:" << std::endl;
  for (const auto& filename : config_file_list) {
    std::cout << "<>    - " << filename << std::endl;
  }

  WaitForEnter();
  std::cout << "<> " << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  WaitForEnter();
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  std::ifstream ifs;
  // 積分器
  class EquationOfMotion<double> rtbp_eom(astro_params);
  auto integrator_rk4 = [&](const my_type::State<double>& state_ptr, double time,
                            double h) -> my_type::State<double> {
    return RungeKutta4Step(rtbp_eom, state_ptr, time, h);
  };
  auto integrator_symplectic = [&](const my_type::State<double>& state_ptr, double time,
                                   double h) -> my_type::State<double> {
    return SymplecticStep4thOrder(kMU, state_ptr, h);
  };
  auto integrator_symplectic6 = [&](const my_type::State<double>& state_ptr, double time,
                                    double h) -> my_type::State<double> {
    return SymplecticStep6thOrder(kMU, state_ptr, h);
  };

  //  実行時間の計測
  auto start_ofall = std::chrono::system_clock::now();

  // while (ifs) {
  for (const auto& configfilepath : config_file_list) {
    std::cout << "<>        Next config file : " << configfilepath << std::endl;
    ifs.open(configfilepath);
    if (!ifs) {
      std::cerr << "<> !err! config : " << configfilepath << " does NOT EXIST" << std::endl;
      continue;
    }
    double progress = 0;
    auto start = std::chrono::system_clock::now();
    std::string str;
    std::cout << std::setprecision(10);
    std::cout << "<>    loaded config >>>" << std::endl;

    //--------- 設定ファイル読み込み部分---------
    // ConfigParserを使用
    utils::ConfigParser config(configfilepath);

    int MESH_CENTER = config.GetInt("MESH CENTER", 0);
    State3d<double> MeshCenter_local{1.0 - kMU, 0, 0};
    if (MESH_CENTER == static_cast<int>(MeshCenter::kSUN)) {
      MeshCenter_local.x = 1.0 - kMU;
    } else if (MESH_CENTER == static_cast<int>(MeshCenter::kEARTH)) {
      MeshCenter_local.x = -kMU;
    }

    State3d<int> MESH_DIVISION = config.GetState3d<int>("MESH DIVISION", {50, 50, 50});
    State3d<double> MESH_HALF_WIDTH = config.GetState3d<double>("MESH SIZE", {0.01, 0.01, 0.01});
    double CALC_TIMESTEP = config.GetDouble("CALC TIMESTEP", 0.001);
    double SALI_CALCTIME_THRESHOLD = config.GetDouble("SALI CALCTIME THRESHOLD", 100.0);
    double SOI_RADIUS = config.GetDouble("RADIUS OF SOI", 0.01);
    double FOREBIDDEN_AREA_RADIUS = config.GetDouble("RADIUS OF FOREBIDDEN AREA", 0.001);
    double JACOBI_INTEGRAL = config.GetDouble("JACOBI INTEGRAL", 3.0);
    double inclination =
        config.GetDouble("INCLINATION AGAINST XY PLANE(deg)", 0.0) * std::acos(-1) / 180.;
    double OMEGA = config.GetDouble("LONGTITUDE AGAINST X AXIS+(deg)", 0.0) * std::acos(-1) / 180.;
    double THETA = config.GetDouble("DEGREE FROM TANGENT", 0.0) * std::acos(-1) / 180.;

    ifs.close();
    std::string center_str = (MESH_CENTER == static_cast<int>(MeshCenter::kSUN)) ? "SUN" : "EARTH";
    std::cout << "<>        MESH CENTER : " << center_str << std::endl;
    std::cout << "<>        MESH DIVISION : " << MESH_DIVISION.x << ", " << MESH_DIVISION.y << ", "
              << MESH_DIVISION.z << std::endl;
    std::cout << "<>        MESH HALF WIDTH : " << MESH_HALF_WIDTH.x << ", " << MESH_HALF_WIDTH.y
              << ", " << MESH_HALF_WIDTH.z << std::endl;
    std::cout << "<>        CALC TIMESTEP : " << CALC_TIMESTEP << std::endl;
    std::cout << "<>        SALI CALCTIME THRESHOLD : " << SALI_CALCTIME_THRESHOLD << std::endl;
    std::cout << "<>        SOI RADIUS : " << SOI_RADIUS << std::endl;
    std::cout << "<>        RADIUS OF FOREBIDDEN AREA : " << FOREBIDDEN_AREA_RADIUS << std::endl;
    std::cout << "<>        JACOBI INTEGRAL : " << JACOBI_INTEGRAL << std::endl;
    std::cout << "<>        INCLINATION(deg) : " << inclination << std::endl;
    std::cout << "<>        LONGTITUDE(deg) : " << OMEGA << std::endl;
    std::cout << "<>        DEGREE FROM TANGENT(deg) : " << THETA << std::endl;
    std::cout << "<>    config file read successfully\n" << std::endl;
    if (is_continuous == 0) {
      WaitForEnter();
    }

    std::cout << "<>    -- Start SALI caluculation for " << configfilepath << std::endl;
    std::cout << "<>        Generating mesh ";
    std::vector<State3d<double>> meshPoints;

    std::cout << "based on SOI radius" << std::endl;

    meshPoints = createDimensionlessCartesianMesh(MeshCenter, MESH_HALF_WIDTH, MESH_DIVISION);

    int countt = meshPoints.size();
    std::cout << "<>        " << countt << " mesh generated successfully" << std::endl;
    std::cout << "<>        Start calclation" << std::endl;

    // 計算のステップ数
    int num_step = static_cast<int>(SALI_CALCTIME_THRESHOLD / CALC_TIMESTEP);
    int totalIterations = meshPoints.size();
    // 進捗カウンタ
    int completed_count = 0;

    const auto& point = meshPoints[1];
    int mesh_num = 1 + 1;
    bool velo_err = 0;
    double final_sali = -1;
    double vx = 0.0, vy = 0.0, vz = 0.0;

    double v_abs = calc_v_abs(point, JACOBI_INTEGRAL, kMU);

    if (v_abs > 0) {
      State3d<double> velocity = calc_velocity(point, v_abs, kMU, inclination, OMEGA, THETA);
      vx = velocity.x;
      vy = velocity.y;
      vz = velocity.z;
    } else {
      std::cout << "<>        velocity error at mesh point " << mesh_num << std::endl;
      velo_err = 1;
    }
    State<double> initial_state_rk4 = {point.x, point.y, point.z, vx, vy, vz};
    State<double> initial_state_symp4 = {point.x, point.y, point.z, vx, vy, vz};
    State<double> initial_state_symp6 = {point.x, point.y, point.z, vx, vy, vz};
    State<double> initial_state_dopri5 = {point.x, point.y, point.z, vx, vy, vz};
    std::vector<std::array<double, 8>> state_history_rk4;
    std::vector<std::array<double, 8>> state_history_symp4;
    std::vector<std::array<double, 8>> state_history_symp6;
    std::vector<std::array<double, 8>> state_history_dopri5;
    StateBufferObserver<double> observer_rk4(state_history_rk4, kMU, JACOBI_INTEGRAL);
    StateBufferObserver<double> observer_symp4(state_history_symp4, kMU, JACOBI_INTEGRAL);
    StateBufferObserver<double> observer_symp6(state_history_symp6, kMU, JACOBI_INTEGRAL);
    StateBufferObserver<double> observer_dopri5(state_history_dopri5, kMU, JACOBI_INTEGRAL);
    // Integrate(initial_state, integrator_symplectic, observer, 0.0, CALC_TIMESTEP, 2);
    // Integrate(initial_state1, integrator_rk4, observer_rk4, 0.0, CALC_TIMESTEP, 2);
    // 積分時間計測
    auto start_integration = std::chrono::high_resolution_clock::now();
    Integrate(initial_state_rk4, integrator_rk4, observer_rk4, 0.0, CALC_TIMESTEP,
              static_cast<int>(SALI_CALCTIME_THRESHOLD / CALC_TIMESTEP));
    auto end_integration = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_integration = end_integration - start_integration;
    std::cout << "<>        RK4 Integration Time: " << elapsed_integration.count() << " seconds"
              << std::endl;
    start_integration = std::chrono::high_resolution_clock::now();
    Integrate(initial_state_symp4, integrator_symplectic, observer_symp4, 0.0, CALC_TIMESTEP,
              static_cast<int>(SALI_CALCTIME_THRESHOLD / CALC_TIMESTEP));
    end_integration = std::chrono::high_resolution_clock::now();
    elapsed_integration = end_integration - start_integration;
    std::cout << "<>        Symplectic 4th Order Integration Time: " << elapsed_integration.count()
              << " seconds" << std::endl;
    start_integration = std::chrono::high_resolution_clock::now();
    Integrate(initial_state_symp6, integrator_symplectic6, observer_symp6, 0.0, CALC_TIMESTEP,
              static_cast<int>(SALI_CALCTIME_THRESHOLD / CALC_TIMESTEP));
    end_integration = std::chrono::high_resolution_clock::now();
    elapsed_integration = end_integration - start_integration;
    std::cout << "<>        Symplectic 6th Order Integration Time: " << elapsed_integration.count()
              << " seconds" << std::endl;
    start_integration = std::chrono::high_resolution_clock::now();
    IntegrateDopri5Orbit(astro_params, initial_state_dopri5, 0.0, SALI_CALCTIME_THRESHOLD,
                         CALC_TIMESTEP, observer_dopri5);
    end_integration = std::chrono::high_resolution_clock::now();
    elapsed_integration = end_integration - start_integration;
    std::cout << "<>        Dopri5 Integration Time: " << elapsed_integration.count() << " seconds"
              << std::endl;
    // // ---------出力ファイル設定---------
    // std::string output_base_path = OUTPUT_DIR;
    // // シミュレーション終了時刻が同じでもファイル名が被らないようにする
    // std::string filename = output_base_path + "/rk4.csv";
    // std::ofstream ofs1(filename);
    // if (!ofs1) {
    //   std::filesystem::path filepath(filename);
    //   std::cerr << "<> !err! output file : " << filepath << " cannot be opened. EXITING..."
    //             << std::endl;
    //   return -1;
    // }
    // ofs1 << std::fixed << std::setprecision(15);
    // for (const auto& state : state_history_rk4) {
    //   ofs1 << state[0] << "," << state[1] << "," << state[2] << "," << state[3] << "," <<
    //   state[4]
    //        << "," << state[5] << "," << state[6] << "," << state[7] << "\n";
    // }
    // ofs1.close();
    // filename = output_base_path + "/symp4.csv";
    // std::ofstream ofs2(filename);
    // if (!ofs2) {
    //   std::filesystem::path filepath(filename);
    //   std::cerr << "<> !err! output file : " << filepath << " cannot be opened. EXITING..."
    //             << std::endl;
    //   return -1;
    // }
    // ofs2 << std::fixed << std::setprecision(15);
    // for (const auto& state : state_history_symp4) {
    //   ofs2 << state[0] << "," << state[1] << "," << state[2] << "," << state[3] << "," <<
    //   state[4]
    //        << "," << state[5] << "," << state[6] << "," << state[7] << "\n";
    // }
    // ofs2.close();
    // filename = output_base_path + "/symp6.csv";
    // std::ofstream ofs3(filename);
    // if (!ofs3) {
    //   std::filesystem::path filepath(filename);
    //   std::cerr << "<> !err! output file : " << filepath << " cannot be opened. EXITING..."
    //             << std::endl;
    //   return -1;
    // }
    // ofs3 << std::fixed << std::setprecision(15);
    // for (const auto& state : state_history_symp6) {
    //   ofs3 << state[0] << "," << state[1] << "," << state[2] << "," << state[3] << "," <<
    //   state[4]
    //        << "," << state[5] << "," << state[6] << "," << state[7] << "\n";
    // }
    // ofs3.close();
    // filename = output_base_path + "/dopri5.csv";
    // std::ofstream ofs4(filename);
    // if (!ofs4) {
    //   std::filesystem::path filepath(filename);
    //   std::cerr << "<> !err! output file : " << filepath << " cannot be opened. EXITING..."
    //             << std::endl;
    //   return -1;
    // }
    // ofs4 << std::fixed << std::setprecision(15);
    // for (const auto& state : state_history_dopri5) {
    //   ofs4 << state[0] << "," << state[1] << "," << state[2] << "," << state[3] << "," <<
    //   state[4]
    //        << "," << state[5] << "," << state[6] << "," << state[7] << "\n";
    // }
    // ofs4.close();

    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    auto msec = duration.count() % 1000;
    auto sec = duration.count() / 1000 % 60;
    auto min = duration.count() / 1000 / 60 % 60;
    auto hour = duration.count() / 1000 / 60 / 60;

    // if (mode2 == '1') {
    //   break;
    // }
    std::cout << "<>        elapsed time : " << hour << "h " << min << "m " << sec << "s " << msec
              << "ms" << std::endl;
    // std::cout << "<>        Simulation for " << configfilename << " finished" << std::endl;
    std::cout << "<>        Simulation for " << configfilepath << " finished" << std::endl;
    // configdata_num++;
    // configfilename =
    //     config_base_path + "/3D_crtbp_SALI/3DSALIconfig_" + std::to_string(configdata_num) +
    //     ".txt";
    // std::cout << "<>        Next config file : " << configfilename << std::endl;
    // ifs.open(configfilename);
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

  std::cout << "<>" << std::endl;
  std::cout << "<>        Calculation finished" << std::endl;
  std::cout << "<>        Total elapsed time : " << hour << "h " << min << "m " << sec << "s "
            << msec << "ms" << std::endl;

  return 0;
}
