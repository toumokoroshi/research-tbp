

#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <crtbp.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector3d.hpp>
#include <vector>

#include "rtbp.hpp"

std::vector<std::streampos> indexFile(const std::string& filename);

std::string readSpecificLine(const std::string& filename,
                             const std::vector<std::streampos>& linePositions, int targetLine);

int main() {
  using namespace param;
  using namespace crtbp;
  using namespace utils;
  // CMakeから渡されたCONFIG_DIRマクロを使用
  std::string config_base_path = CONFIG_DIR;
  std::string astro_param_file = config_base_path + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);
  const double kAU = astro_params.au;                 // astronomical unit in meters
  const double kGMSUN = astro_params.gm_sun;          // heliocentric gravitational constant m3 s-2
  const double kGMEARTH = astro_params.gm_earth;      // geocentric gravitational constant m3 s-2
  const double kMU = kGMEARTH / (kGMEARTH + kGMSUN);  // mu parameter of Earth-Sun
  constexpr double kPERTUBATION = 1e-10;

  constexpr double HEADER_SIZE = 9;
  my_type::State3d<double> MeshCenter{1.0 - kMU, 0, 0};
  double ROI_length = 0;

  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP 3dSALI Calculation based on Jacobi Integral" << std::endl;
  std::cout << "<>-------------------------------------------------------------"
               "---\n\n"
            << std::endl;
  std::cout << "<>****************************************************************" << std::endl;
  std::cout << "<>  [mode selection] : " << std::endl;
  std::cout << "<>        1. New simulation" << std::endl;
  std::cout << "<>        2. Detailed simulation for existing data" << std::endl;
  std::cout << "<>        else. Exit" << std::endl;
  std::cout << "<>   enter number " << std::endl;
  std::cout << "<> >>>";
  char mode;
  std::cin >> mode;
  std::cout << std::endl;
  if (mode == '1') {
    std::cout << "<> > selected mode : New simulation\n" << std::endl;
  } else if (mode == '2') {
    std::cout << "<> > selected mode : Detailed simulation for existing data\n" << std::endl;
  } else {
    std::cout << "<> > selected mode : Exit\n" << std::endl;
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
  char mode2;
  std::cout << "<>  [single simulation or continuous simulation] : " << std::endl;
  std::cout << "<>        1. single simulation" << std::endl;
  std::cout << "<>        2. continuous simulation" << std::endl;
  std::cout << "<>        else. Exit" << std::endl;
  std::cout << "<>      enter number " << std::endl;
  std::cout << "<> >>> ";
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

  if (mode2 == '1') {
    configfilename = config_base_path + "/3D_crtbp_SALI/3DSALIconfig.txt";
  } else if (mode2 == '2') {
    // configファイルの数だけ計算する
    configfilename = config_base_path + "/3D_crtbp_SALI/3DSALIconfig_1.txt";
  }

  std::cout << "<>        config file : " << configfilename << std::endl;
  ifs.open(configfilename);

  if (!ifs) {
    std::cerr << "Failed to open file." << std::endl;
    return -1;
  }
  int Core_Max = omp_get_max_threads();
  int OMP_Fmax{};
  std::cout << "<>  [OpenMP preparation]" << std::endl;
  std::cout << "<> On your PC, " << Core_Max
            << " threads can be used for parallel computing employing OMP." << std::endl;
  std::cout << "<> >  " << std::endl;
  std::cout << "<>   * How many threads do you want use for simulation? "
            << "(input an integer)" << std::endl;
  std::cout << "<>   * （※最大コア数を指定すると計算が終わるまでPCが" << std::endl;
  std::cout << "<>   *    激重になるので，最大値-1くらいが良いかも？）> " << std::endl;
  std::cout << "<> >>> ";
  int getcore = 0;
  std::cin >> getcore;
  if (getcore <= Core_Max) {
    OMP_Fmax = getcore;
    std::cout << "  <>     >> Number of OMP threads is " << OMP_Fmax << std::endl;
  } else {
    OMP_Fmax = Core_Max;
    std::cout << "  <>     >> Your input is INVALID. OMP threads is "
              << "automatically determined as " << OMP_Fmax << std::endl;
  }
  WaitForEnter();
  omp_set_num_threads(OMP_Fmax);
  // 積分器
  auto integrator = [&](SaliState<double>* state_ptr, double h) {
    SymplecticStep4thOrderSALI(kMU, state_ptr, h);
  };
  int configdata_num = 1;
  //  実行時間の計測
  auto start_ofall = std::chrono::system_clock::now();

  // --------  configファイルの数だけSALI計算全体を繰り返す-----------------------
  while (ifs) {
    double progress = 0;
    auto start = std::chrono::system_clock::now();
    std::string str;
    std::cout << std::setprecision(10);

    //--------- 設定ファイル読み込み部分---------
    while (std::getline(ifs, str)) {
      if (str.find("MESH SIZE") != std::string::npos) {
        MESH_SIZE = std::stoi(str.substr(str.find("=") + 1));
        std::cout << "<>        MESH SIZE : " << MESH_SIZE << std::endl;
      } else if (str.find("CALC TIMESTEP") != std::string::npos) {
        CALC_TIMESTEP = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        CALC TIMESTEP : " << CALC_TIMESTEP << std::endl;
      } else if (str.find("SALI CALCTIME THRESHOLD") != std::string::npos) {
        SALI_CALCTIME_THRESHOLD = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        SALI CALCTIME THRESHOLD : " << SALI_CALCTIME_THRESHOLD << std::endl;
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
        inclination = inclination * std::acos(-1) / 180.;
      } else if (str.find("LONGTITUDE AGAINST X AXIS+(deg)") != std::string::npos) {
        OMEGA = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        LONGTITUDE(deg) : " << OMEGA << std::endl;
        OMEGA = OMEGA * std::acos(-1) / 180.;
      } else if (str.find("DEGREE FROM TANGENT") != std::string::npos) {
        THETA = std::stod(str.substr(str.find("=") + 1));
        std::cout << "<>        DEGREE FROM TANGENT(deg) : " << THETA << std::endl;
        THETA = THETA * std::acos(-1) / 180.;
      }
    }
    ifs.close();
    std::cout << "<>    config file read successfully\n" << std::endl;
    if (is_continuous == 0) {
      std::cout << "<>  [read config validation]" << std::endl;
      WaitForEnter();
      std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "<>    Start SALI caluculation --" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "<>        Generating mesh ";
    std::vector<State3d<double>> meshPoints;
    if (mode == '1') {
      std::cout << "based on SOI radius" << std::endl;
      meshPoints = createSphereMesh(SOI_RADIUS, MESH_SIZE, MeshCenter);
    } else if (mode == '2') {
      // std::cout << "based on the specified point" << std::endl;
      // meshPoints = create_cube_mesh(ROI_length, MESH_SIZE, MeshCenter);
    }
    int countt = meshPoints.size();
    std::cout << std::endl;
    std::cout << "<>        " << countt << " mesh generated successfully" << std::endl;
    std::cout << std::endl;
    std::cout << "<>        Start calclation" << std::endl;

    // ---------出力ファイル設定---------
    std::string output_base_path = OUTPUT_DIR;
    // シミュレーション終了時刻が同じでもファイル名が被らないようにする
    std::string filename = output_base_path + "/3D_crtbp_SALI/configdata_" +
                           std::to_string(configdata_num) + "_" + getcurrent_date() + ".txt";
    std::ofstream ofs1(filename);
    if (!ofs1) {
      std::filesystem::path filepath(filename);
      std::filesystem::path dir = filepath.parent_path();
      if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
      }
      ofs1.open(filename);
      if (!ofs1) {
        std::cerr << "Failed to open file even after creating directory: " << filename << std::endl;
        return -1;
      }
    }
    // ヘッダーを書き込む
    //     ofs1 << "MESH SIZE=" << MESH_SIZE << std::endl;
    ofs1 << "CALCULATION TIMESTEP=" << CALC_TIMESTEP << std::endl;
    ofs1 << "SIMULATION TIME=" << SALI_CALCTIME_THRESHOLD << std::endl;
    ofs1 << "RADIUSofSOI=" << SOI_RADIUS << std::endl;
    ofs1 << "FOREBIDDEN AREA RADIUS=" << FOREBIDDEN_AREA_RADIUS << std::endl;
    ofs1 << "INITIAL JACOBI INTEGRAL=" << JACOBI_INTEGRAL << std::endl;
    ofs1 << "INCLINATION AGAINST XY PLANE=" << inclination / std::acos(-1) * 180. << std::endl;
    ofs1 << "LONGTITUDE AGAINST X AXIS=" << OMEGA / std::acos(-1) * 180. << std::endl;
    ofs1 << "DEGREE FROM TANGENT(deg)=" << THETA / std::acos(-1) * 180. << std::endl;
    ofs1 << "Time,SALI,x,y,z,px,py,pz\n";

    // 積分クラスに渡すためのファイル出力用オブザーバー
    SaliFileObserver<double> observer(ofs1);
    // 計算のステップ数
    int num_step = static_cast<int>(SALI_CALCTIME_THRESHOLD / CALC_TIMESTEP);
    int totalIterations = meshPoints.size();
    // 進捗カウンタ
    int completed_count = 0;
    std::vector<std::array<double, 8>> ref_state;
    std::vector<std::array<double, 8>> pertubed_state1;
    std::vector<std::array<double, 8>> pertubed_state2;

// OpenMP並列化ループ
#pragma omp parallel shared(meshPoints, completed_count, totalIterations, progress, ofs1)
    {
      // 各スレッドがSALIの結果を一時的に保存する文字列
      std::stringstream local_output_buffer;
      local_output_buffer << std::fixed << std::setprecision(15);
#pragma omp for schedule(dynamic)
      for (int idx = 0; idx < static_cast<int>(meshPoints.size()); ++idx) {
        const auto& point = meshPoints[idx];

        int mesh_num = idx + 1;
        bool velo_err = 0;

        // non-biased velocity
        double v_abs = calc_v_abs(point, JACOBI_INTEGRAL, kMU);
        double vx = 0.0, vy = 0.0, vz = 0.0;

        if (v_abs > 0) {
          State3d<double> velocity = calc_velocity(point, v_abs, kMU, inclination, OMEGA, THETA);
          vx = velocity.x;
          vy = velocity.y;
          vz = velocity.z;
        } else {
          velo_err = 1;
        }

        if (velo_err) {
          completed_count++;
          continue;
        }
        State<double> initial_state = {point.x, point.y, point.z, vx, vy, vz};
        SaliState<double> sali_state;
        sali_state.state = ConvertToCanonical(initial_state);
        // 偏差ベクトル w1, w2 の初期化
        sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        // 積分ループ (オブザーバー無し)
        for (int step = 0; step < num_step; ++step) {
          // 1. 積分
          integrator(&sali_state, CALC_TIMESTEP);
          // 2. 正規化
          sali_state.w1.Normalize();
          sali_state.w2.Normalize();
        }
        const double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
        const double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
        const double final_sali = std::min(norm_plus, norm_minus);
        local_output_buffer << mesh_num << "," << final_sali << "," << point.x << "," << point.y
                            << "," << point.z << "," << vx << "," << vy << "," << vz << "\n";
        // 5. 一定件数ごとにバッファをファイルに書き込む (排他制御)
        // (またはループ終了後にまとめて書き込む)
        if (idx % 100 == 0 || idx == totalIterations - 1) {
#pragma omp critical
          {
            ofs1 << local_output_buffer.str();
            local_output_buffer.str("");  // バッファをクリア
          }
        }
#pragma omp atomic
        completed_count++;

#pragma omp critical
        {
          int display_interval = std::max(totalIterations / 100, 1);
          if (completed_count % display_interval == 0 || completed_count == totalIterations) {
            double current_progress = static_cast<double>(completed_count) / totalIterations;

            // (注: この関数が内部で std::cout を使う前提)
            displayProgressBarThreadSafe(current_progress);
          }
        }
      }
    }  // end of parallel region

    // 最終的な進捗バー表示
    displayProgressBarThreadSafe(1.0);
    std::cout << std::endl;

    ofs1.close();

    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    auto msec = duration.count() % 1000;
    auto sec = duration.count() / 1000 % 60;
    auto min = duration.count() / 1000 / 60 % 60;
    auto hour = duration.count() / 1000 / 60 / 60;

    if (mode2 == '1') {
      break;
    }
    std::cout << "<>        elapsed time : " << hour << "h " << min << "m " << sec << "s " << msec
              << "ms" << std::endl;
    std::cout << "<>        Simulation for " << configfilename << " finished" << std::endl;
    configdata_num++;
    configfilename =
        config_base_path + "/3D_crtbp_SALI/3DSALIconfig_" + std::to_string(configdata_num) + ".txt";
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
                             const std::vector<std::streampos>& linePositions, int targetLine) {
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
