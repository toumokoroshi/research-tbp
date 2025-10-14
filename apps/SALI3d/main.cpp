/** @file           main_sali.cpp
    @brief
   DE441を用いず、小惑星の速度を仮定し、ヤコビ積分を基準にCRTBPにおいてSALIを計算する
    @author         tabata
    @date           2024/10/22
    @par            edittor/      date/ version/ description
                    tabata/ 2024/11/20/     1.0/ 初版作成
                    tabata/ 2025/1/28/     1.1/
 */

#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <crtbp.hpp>
#include <utils.hpp>
#include <vector3d.hpp>

#define _DEBUG
// #define SALI_only_XY

// CMakeから渡されるマクロを文字列に変換するためのヘルパーマクロ
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

// 任意の2組の軌道要素からSALIを計算する関数
double calc_SALI(const std::array<double, 6> &ref_state,
                 const std::array<double, 6> &perturbed_state1,
                 const std::array<double, 6> &perturbed_state2, int mode = 6);
double calc_r1(const Point3D &point, const double mu);

double calc_r2(const Point3D &point, const double mu);

double calc_v_abs(const Point3D &point, const double mu,
                  const double JACOBI_INTEGRAL);

double calc_jacobi_integral(const std::array<double, 6> &state,
                            const double mu);

Vector3d calc_velocity(const Point3D &point, const double v_abs,
                       const double mu, const double inclination,
                       const double OMEGA, const double theta = 0.0);

std::vector<std::streampos> indexFile(const std::string &filename);

std::string readSpecificLine(const std::string &filename,
                             const std::vector<std::streampos> &linePositions,
                             int targetLine);

int main() {
  constexpr double MU = 3.003e-6; // 地球と太陽の質量比
  constexpr double HEADER_SIZE = 9;
  Point3D MeshCenter(1.0 - MU, 0, 0);
  double ROI_length = 0;

  std::cout
      << "<>----------------------------------------------------------------"
      << std::endl;
  std::cout << "<>            CRTBP 3dSALI Calculation based on Jacobi Integral"
            << std::endl;
  std::cout << "<>-------------------------------------------------------------"
               "---\n\n\n\n"
            << std::endl;
  // #ifndef _DEBUG
  std::cout
      << "<>****************************************************************"
      << std::endl;
  std::cout << "<>  [mode] : " << std::endl;
  std::cout << "<>        1. New simulation" << std::endl;
  std::cout << "<>        2. Detailed simulation for existing data"
            << std::endl;
  std::cout << "<>        else. Exit" << std::endl;
  std::cout << "<>      Input number : ";
  char mode;
  std::cin >> mode;
  std::cout << std::endl;
  if (mode == '1') {
    std::cout << "selected mode : New simulation\n" << std::endl;
  } else if (mode == '2') {
    std::cout << "selected mode : Detailed simulation for existing data\n"
              << std::endl;
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
    int targetLine =
        std::stoi(mesh_num_of_interest) + HEADER_SIZE; // 読み込みたい行番号
    std::string line =
        readSpecificLine(filename_interested, linePositions, targetLine);
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
  std::cout
      << "<>****************************************************************"
      << std::endl;
  std::cout << "<>           --simulation config--\n" << std::endl;
  std::cout << "<>    reading config file\n" << std::endl;

  char mode2;
  std::cout << "<>  [single simulation or continuous simulation] : "
            << std::endl;
  std::cout << "<>        1. single simulation" << std::endl;
  std::cout << "<>        2. continuous simulation" << std::endl;
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
  } else {
    std::cout << "selected mode : Exit\n" << std::endl;
    return 0;
  }

  constexpr double SOI = 0.03;
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

  // CMakeから渡されたCONFIG_DIRマクロを使用
  std::string config_base_path = CONFIG_DIR;

  if (mode2 == '1') {
    configfilename = config_base_path + "/3D_crtbp_SALI/3DSALIconfig.txt";
  } else if (mode2 == '2') {
    configfilename = config_base_path + "/3D_crtbp_SALI/3DSALIconfig_1.txt";
  }

  std::cout << "<>        config file : " << configfilename << std::endl;
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
        std::cout << "<>        SALI CALCTIME THRESHOLD : "
                  << SALI_CALCTIME_THRESHOLD << std::endl;
      } else if (str.find("RADIUS OF SOI") != std::string::npos) {
        SOI_RADIUS = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        SOI RADIUS : " << SOI_RADIUS << std::endl;
      } else if (str.find("RADIUS OF FOREBIDDEN AREA") != std::string::npos) {
        FOREBIDDEN_AREA_RADIUS = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        RADIUS OF FOREBIDDEN AREA : "
                  << FOREBIDDEN_AREA_RADIUS << std::endl;
      } else if (str.find("JACOBI INTEGRAL") != std::string::npos) {
        JACOBI_INTEGRAL = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        JACOBI INTEGRAL : " << JACOBI_INTEGRAL
                  << std::endl;
      } else if (str.find("INCLINATION AGAINST XY PLANE(deg)") !=
                 std::string::npos) {
        inclination = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        INCLINATION(deg) : " << inclination
                  << std::endl;
        inclination = inclination * std::acos(-1) / 180.; // deg to rad
      } else if (str.find("LONGTITUDE AGAINST X AXIS+(deg)") !=
                 std::string::npos) {
        OMEGA = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        LONGTITUDE(deg) : " << OMEGA << std::endl;
        OMEGA = OMEGA * std::acos(-1) / 180.; // deg to rad
      } else if (str.find("DEGREE FROM TANGENT") != std::string::npos) {
        THETA = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        DEGREE FROM TANGENT(deg) : " << THETA
                  << std::endl;
        THETA = THETA * std::acos(-1) / 180.; // deg to rad
      }
    }
    ifs.close();
    std::cout << "<>    config file read successfully\n" << std::endl;
    if (is_continuous == 0) {
      std::cout << "<>  [Check the param and Press any key to go ahead]"
                << std::endl;
      __KEYWAIT__
    }

    std::cout << std::endl;
    std::cout << "<>           --SALI caluculation --" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "<>        Generating mesh ";
    std::vector<Point3D> meshPoints;
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
    std::cout << "<>        " << countt << " mesh generated successfully"
              << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "<>        Start calclation" << std::endl;

    /* メッシュ上の点を用いてSALIを計算 */
    constexpr double perturbation = 1e-10;
#ifdef SALI_only_XY
    std::vector<std::array<double, 7>> SALI_data;
#else
    std::vector<std::array<double, 10>> SALI_data;
#endif
    // [0]:count, [1]time , [2]:x, [3]:y, [4]:z,
    // [5]:jacobi constant
    // [6]:SALI(x-biased&y-biased), [7]:SALI(x-biased&z-biased),
    // [8]:SALI(y-biased&z-biased)
    SALI_data.reserve(countt);

    int totalIterations = meshPoints.size();
    for (const auto &point : meshPoints) {
      displayProgressBar(progress);
      // 0:計算終了, 1:計算継続
      bool calc_traj = 1;
      // 0:計算継続,
      // 1:計算中断,途中でhill球から抜け出したり、禁止領域に入った場合
      bool abort_calc = 0;
      // 速度が定義できない場合（速度の絶対値が負になる場合）
      bool velo_err = 0;

      static int mesh_num = 1; // メッシュ番号

      // non-biased velocity
      double v_abs = calc_v_abs(point, MU, JACOBI_INTEGRAL);
      double vx = 0.0, vy = 0.0, vz = 0.0;

      if (v_abs > 0) {
        Vector3d velocity =
            calc_velocity(point, v_abs, MU, inclination, OMEGA, THETA);
        vx = velocity.x();
        vy = velocity.y();
        vz = velocity.z();
      } else {
        calc_traj = 0;
        velo_err = 1;
      }

      CRTBP ref_state(point.x, point.y, point.z, vx, vy, vz, MU, CALC_TIMESTEP);
      CRTBP perturbed_state1(point.x + perturbation, point.y, point.z, vx, vy,
                             vz, MU, CALC_TIMESTEP);
      CRTBP perturbed_state2(point.x, point.y + perturbation, point.z, vx, vy,
                             vz, MU, CALC_TIMESTEP);

#ifndef SALI_only_XY
      CRTBP perturbed_state3(point.x, point.y, point.z + perturbation, vx, vy,
                             vz, MU, CALC_TIMESTEP);
#endif
      CRTBP perturbed_state4(point.x, point.y, point.z, vx, vy,
                             vz + perturbation, MU, CALC_TIMESTEP);

      double time = 0.0;
#ifndef SALI_only_XY
      double SALIxy = -1.0;
      double SALIyz = -1.0;
      double SALIzx = -1.0;
#endif
#ifdef SALI_only_XY
      double SALI = -1.0;
#endif
      if (!velo_err) {
        while (calc_traj) {
          if (time > SALI_CALCTIME_THRESHOLD) {
            calc_traj = 0;
            continue;
          }
          if (ref_state.calc_r2() > SOI ||
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
      std::array<double, 7> SALI_data1{static_cast<double>(mesh_num),
                                       time,
                                       point.x,
                                       point.y,
                                       point.z,
                                       jacobiii,
                                       SALI};
#endif
#ifndef SALI_only_XY
      // 最小値
      double SALI = std::min({SALIxy, SALIyz, SALIzx});
      std::array<double, 10> SALI_data1{static_cast<double>(mesh_num),
                                        time,
                                        point.x,
                                        point.y,
                                        point.z,
                                        jacobiii,
                                        SALIxy,
                                        SALIyz,
                                        SALIzx,
                                        SALI};
#endif
      SALI_data.emplace_back(SALI_data1);
      // 進捗を計算
      progress = (mesh_num + 1.0) / totalIterations;

      // プログレスバーを表示
      displayProgressBar(progress);

      mesh_num++;
    }

    //  日付と時刻をファイル名にしてファイル出力
    // CMakeから渡されたOUTPUT_DIRマクロを使用
    std::string output_base_path = OUTPUT_DIR;
    std::string filename =
        output_base_path + "/SALI/3DSALI_" + getcurrent_date() + ".dat";
    std::ofstream ofs1(filename);
    if (!ofs1) {
      // ファイルが開けなかった場合、ディレクトリが存在しない可能性があるため作成を試みる
      std::filesystem::path filepath(filename);
      std::filesystem::path dir = filepath.parent_path();
      if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
      }
      // 再度ファイルを開く
      ofs1.open(filename);
      if (!ofs1) {
        std::cerr << "Failed to open file even after creating directory: "
                  << filename << std::endl;
        return -1;
      }
    }
    // 計算のコンフィグレーションをファイルに書き込む
    ofs1 << "MESH SIZE=" << MESH_SIZE << std::endl;
    ofs1 << "CALCULATION TIMESTEP=" << CALC_TIMESTEP << std::endl;
    ofs1 << "SIMULATION TIME=" << SALI_CALCTIME_THRESHOLD << std::endl;
    ofs1 << "RADIUSofSOI=" << SOI_RADIUS << std::endl;
    ofs1 << "FOREBIDDEN AREA RADIUS=" << FOREBIDDEN_AREA_RADIUS << std::endl;
    ofs1 << "INITIAL JACOBI INTEGRAL=" << JACOBI_INTEGRAL << std::endl;
    ofs1 << "INCLINATION AGAINST XY PLANE="
         << inclination / std::acos(-1) * 180. << std::endl; // rad to deg
    ofs1 << "LONGTITUDE AGAINST X AXIS=" << OMEGA / std::acos(-1) * 180.
         << std::endl; // rad to
                       // deg
    ofs1 << "DEGREE FROM TANGENT(deg)=" << THETA / std::acos(-1) * 180.
         << std::endl;
#ifdef SALI_only_XY
    ofs1 << "mesh_num, time, x, y, z, jacobi constant, SALI" << std::endl;
#endif
#ifndef SALI_only_XY
    ofs1 << "mesh_num, time, x, y, z, jacobi constant, SALIxy,SALIyz,SALIxz,"
         << std::endl;
#endif

    for (const auto &data : SALI_data) {
#ifdef SALI_only_XY
      ofs1 << std::setprecision(0) << std::fixed << data[0] << " "
           << std::setprecision(4) << data[1] << " " << std::setprecision(15)
           << std::fixed << data[2] << " " << data[3] << " " << data[4] << " "
           << data[5] << " " << data[6] << " " << std::endl;
#endif
#ifndef SALI_only_XY
      ofs1 << std::setprecision(0) << std::fixed << data[0] << " "
           << std::setprecision(4) << data[1] << " " << std::setprecision(15)
           << std::fixed << data[2] << " " << data[3] << " " << data[4] << " "
           << data[5] << " " << data[6] << " " << data[7] << " " << data[8]
           << " " << data[9] << std::endl;
#endif
    }
    ofs1.close();
    auto end = std::chrono::system_clock::now();

    // 実行時間の計算
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

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
    std::cout << "<>        elapsed time : " << hour << "h " << min << "m "
              << sec << "s " << msec << "ms" << std::endl;
    std::cout << "<>        Simulation for " << configfilename << " finished"
              << std::endl;
    configdata_num++;
    configfilename = config_base_path + "/3D_crtbp_SALI/3DSALIconfig_" +
                     std::to_string(configdata_num) + ".txt";
    std::cout << "<>        Next config file : " << configfilename << std::endl;
    ifs.open(configfilename);
  }

  // 実行時間の計測
  auto end = std::chrono::system_clock::now();

  // 実行時間の計算
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start_ofall);

  // 時、分、秒、ミリ秒に分解
  auto msec = duration.count() % 1000;
  auto sec = duration.count() / 1000 % 60;
  auto min = duration.count() / 1000 / 60 % 60;
  auto hour = duration.count() / 1000 / 60 / 60;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "<>        Calculation finished" << std::endl;
  std::cout << "<>        Total elapsed time : " << hour << "h " << min << "m "
            << sec << "s " << msec << "ms" << std::endl;

  return 0;
}

double calc_SALI(const std::array<double, 6> &ref_state,
                 const std::array<double, 6> &perturbed_state1,
                 const std::array<double, 6> &perturbed_state2, int mode) {
  double SALI = 5.0;
  if (mode == 3) {
    std::array<double, 3> q_ref{ref_state[0], ref_state[1], ref_state[2]};
    std::array<double, 3> q_pertubed1{perturbed_state1[0], perturbed_state1[1],
                                      perturbed_state1[2]};
    std::array<double, 3> q_pertubed2{perturbed_state2[0], perturbed_state2[1],
                                      perturbed_state2[2]};
    std::array<double, 3> p_ref{ref_state[3] - ref_state[0],
                                ref_state[4] + ref_state[1], ref_state[5]};
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

    SALI = std::min(sa_norm, wa_norm);
  } else if (mode == 6) {
    // refとの差分を計算
    std::array<double, 6> diff1{
        perturbed_state1[0] - ref_state[0], perturbed_state1[1] - ref_state[1],
        perturbed_state1[2] - ref_state[2], perturbed_state1[3] - ref_state[3],
        perturbed_state1[4] - ref_state[4], perturbed_state1[5] - ref_state[5]};
    std::array<double, 6> diff2{
        perturbed_state2[0] - ref_state[0], perturbed_state2[1] - ref_state[1],
        perturbed_state2[2] - ref_state[2], perturbed_state2[3] - ref_state[3],
        perturbed_state2[4] - ref_state[4], perturbed_state2[5] - ref_state[5]};
    // 差分ベクトルを正規化
    double norm1 = std::sqrt(diff1[0] * diff1[0] + diff1[1] * diff1[1] +
                             diff1[2] * diff1[2] + diff1[3] * diff1[3] +
                             diff1[4] * diff1[4] + diff1[5] * diff1[5]);
    double norm2 = std::sqrt(diff2[0] * diff2[0] + diff2[1] * diff2[1] +
                             diff2[2] * diff2[2] + diff2[3] * diff2[3] +
                             diff2[4] * diff2[4] + diff2[5] * diff2[5]);
    if (norm1 == 0 || norm2 == 0) {
      // 差分ベクトルがゼロの場合、SALIは定義できない
      return -1.0;
    }
    std::array<double, 6> normalized_diff1{diff1[0] / norm1, diff1[1] / norm1,
                                           diff1[2] / norm1, diff1[3] / norm1,
                                           diff1[4] / norm1, diff1[5] / norm1};
    std::array<double, 6> normalized_diff2{diff2[0] / norm2, diff2[1] / norm2,
                                           diff2[2] / norm2, diff2[3] / norm2,
                                           diff2[4] / norm2, diff2[5] / norm2};

    std::array<double, 6> d_plus{normalized_diff1[0] + normalized_diff2[0],
                                 normalized_diff1[1] + normalized_diff2[1],
                                 normalized_diff1[2] + normalized_diff2[2],
                                 normalized_diff1[3] + normalized_diff2[3],
                                 normalized_diff1[4] + normalized_diff2[4],
                                 normalized_diff1[5] + normalized_diff2[5]};
    std::array<double, 6> d_minus{normalized_diff1[0] - normalized_diff2[0],
                                  normalized_diff1[1] - normalized_diff2[1],
                                  normalized_diff1[2] - normalized_diff2[2],
                                  normalized_diff1[3] - normalized_diff2[3],
                                  normalized_diff1[4] - normalized_diff2[4],
                                  normalized_diff1[5] - normalized_diff2[5]};

    // d_plusとd_minusのノルムを計算
    double d_plus_norm = std::sqrt(
        d_plus[0] * d_plus[0] + d_plus[1] * d_plus[1] + d_plus[2] * d_plus[2] +
        d_plus[3] * d_plus[3] + d_plus[4] * d_plus[4] + d_plus[5] * d_plus[5]);
    double d_minus_norm =
        std::sqrt(d_minus[0] * d_minus[0] + d_minus[1] * d_minus[1] +
                  d_minus[2] * d_minus[2] + d_minus[3] * d_minus[3] +
                  d_minus[4] * d_minus[4] + d_minus[5] * d_minus[5]);
    // d_plusとd_minusのノルムを比較して最小値を返す
    // SALIの計算
    SALI = std::min(d_plus_norm, d_minus_norm);
  }
  return SALI;
}

double calc_r1(const Point3D &point, const double mu) {
  return std::sqrt(std::pow(point.x + mu, 2.) + std::pow(point.y, 2.) +
                   std::pow(point.z, 2.));
}

double calc_r2(const Point3D &point, const double mu) {
  return std::sqrt(std::pow(point.x - 1. + mu, 2.) + std::pow(point.y, 2.) +
                   std::pow(point.z, 2.));
}

double calc_v_abs(const Point3D &point, const double mu,
                  const double JACOBI_INTEGRAL) {
  double r1 = calc_r1(point, mu);
  double r2 = calc_r2(point, mu);
  // std::cout << "uiuoij" <<point.x * point.x + point.y * point.y + 2. * (1. -
  // mu) / r1 + 2.
  // * mu / r2 +
  //                  mu * (1 - mu) - JACOBI_INTEGRAL << std::endl;
  return std::sqrt(point.x * point.x + point.y * point.y + 2. * (1. - mu) / r1 +
                   2. * mu / r2 + mu * (1. - mu) - JACOBI_INTEGRAL);
}

double calc_jacobi_integral(const std::array<double, 6> &state,
                            const double mu) {
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
  return x * x + y * y - (vx * vx + vy * vy + vz * vz) + 2. * (1. - mu) / r1 +
         2. * mu / r2 + mu * (1. - mu);
}

Vector3d calc_velocity(const Point3D &point, const double v_abs,
                       const double mu, const double inclination,
                       const double OMEGA, const double theta) {
  // 回転軸と回転角からクオータニオン経由で回転行列を生成
  auto create_rot_matrix =
      [](const Vector3d &unit_n,
         double theta_) -> std::array<std::array<double, 3>, 3> {
    double half_theta = theta_ / 2.0;
    double q0 = std::cos(half_theta);
    double sin_half_theta = std::sin(half_theta);
    double q1 = unit_n.x() * sin_half_theta;
    double q2 = unit_n.y() * sin_half_theta;
    double q3 = unit_n.z() * sin_half_theta;

    double q0q0 = q0 * q0;
    double q1q1 = q1 * q1;
    double q2q2 = q2 * q2;
    double q3q3 = q3 * q3;
    double q0q1 = q0 * q1;
    double q0q2 = q0 * q2;
    double q0q3 = q0 * q3;
    double q1q2 = q1 * q2;
    double q1q3 = q1 * q3;
    double q2q3 = q2 * q3;

    std::array<std::array<double, 3>, 3> rot_matrix = {
        {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 - q0q3), 2.0 * (q1q3 + q0q2)},
         {2.0 * (q1q2 + q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 - q0q1)},
         {2.0 * (q1q3 - q0q2), 2.0 * (q2q3 + q0q1),
          q0q0 - q1q1 - q2q2 + q3q3}}};
    return rot_matrix;
  };

  // 入力されたベクトルを入力された行列で変換
  auto convert = [](std::array<std::array<double, 3>, 3> convert_matrix,
                    const Vector3d &v) -> Vector3d {
    return {convert_matrix[0][0] * v.x() + convert_matrix[0][1] * v.y() +
                convert_matrix[0][2] * v.z(),
            convert_matrix[1][0] * v.x() + convert_matrix[1][1] * v.y() +
                convert_matrix[1][2] * v.z(),
            convert_matrix[2][0] * v.x() + convert_matrix[2][1] * v.y() +
                convert_matrix[2][2] * v.z()};
  };
  double vx_, vy_, vz_ = 0;

  // inclinationとOMEGAを用いて軌道面の法線ベクトルを計算
  Vector3d normal_vector{std::sin(inclination) * std::cos(OMEGA),
                         std::sin(inclination) * std::sin(OMEGA),
                         std::cos(inclination)};
  // 法線ベクトルと位置ベクトルの外積を計算
  Vector3d r2_vector{point.x - 1. + mu, point.y, point.z};
  Vector3d h_vector = r2_vector.gaiseki(normal_vector);
  Vector3d normalized_h_vector = h_vector.normalise(); // theta = 0は逆行回転

  // 速度ベクトルを計算
  vx_ = v_abs * normalized_h_vector.x();
  vy_ = v_abs * normalized_h_vector.y();
  vz_ = v_abs * normalized_h_vector.z();

  if (theta == 0.0)
    return Vector3d(vx_, vy_, vz_);
  else {
    // クオータニオンを用いて速度ベクトルをnormal_vector周りにthetaだけ回転させる
    std::array<std::array<double, 3>, 3> rot_matrix =
        create_rot_matrix(normal_vector, theta);
    Vector3d velocity{vx_, vy_, vz_};
    Vector3d rotated_velocity = convert(rot_matrix, velocity);

    return rotated_velocity;
  }
}

std::vector<std::streampos> indexFile(const std::string &filename) {
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

std::string readSpecificLine(const std::string &filename,
                             const std::vector<std::streampos> &linePositions,
                             int targetLine) {
  if (targetLine < 1 || targetLine >= static_cast<int>(linePositions.size())) {
    throw std::out_of_range("指定した行が範囲外です: " +
                            std::to_string(targetLine));
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