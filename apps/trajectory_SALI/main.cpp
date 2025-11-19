

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
 * @brief クラスの概要説明
 *
 * クラスの詳細説明
 */
template <typename ScalarType>
class TrajectoryObserver {
 public:
  /**
   * @brief コンストラクタ
   */
  TrajectoryObserver(std::vector<std::array<ScalarType, 8>>& history, ScalarType jacobi_constant,
                     ScalarType mu)
      : history_(history), jacobi_constant_(jacobi_constant), mu_(mu) {};

  /**
   * @brief デストラクタ
   */
  virtual ~TrajectoryObserver();

  TrajectoryObserver& operator=(const TrajectoryObserver&) = delete;

  void operator()(const State<ScalarType>& state, ScalarType t) {
    history_.push_back({t, jacobi_constant_ - calc_jacobi_integral(state, mu_), state.x, state.y,
                        state.z, state.vx, state.vy, state.vz});
  }

 private:
  std::vector<std::array<ScalarType, 8>>& history_;
  ScalarType jacobi_constant_;
  ScalarType mu_;
};

int main() {
  using namespace param;
  using namespace crtbp;
  using namespace utils;
  namespace fs = std::filesystem;
  // CMakeから渡されたCONFIG_DIRマクロを使用
  const std::string kConfigFilePath = CONFIG_DIR;
  const std::string kCalcConfigPath = kConfigFilePath + "/3D_crtbp_SALI/";
  std::string calc_config_prefix = "3DSALIconfig";
  std::string astro_param_file = kConfigFilePath + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);

  const double kAU = astro_params.au;                 // astronomical unit in meters
  const double kGMSUN = astro_params.gm_sun;          // heliocentric gravitational constant m3 s-2
  const double kGMEARTH = astro_params.gm_earth;      // geocentric gravitational constant m3 s-2
  const double kMU = kGMEARTH / (kGMEARTH + kGMSUN);  // mu parameter of Earth-Sun
  std::cout << "-" << std::endl;

  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            SALI calc on spesific trajectory" << std::endl;
  std::cout << "<>-------------------------------------------------------------"
               "---\n\n"
            << std::endl;
  std::cout << "<>****************************************************************" << std::endl;
  State<double> ast_state_ssb{};
  State<double> sun_state_ssb{};
  State<double> earth_state_ssb{};
  State<double> init_ast_state_crtbp{};
  /* ファイルから軌道要素と計算コンフィグをよんで表示
  →座標変換
  →座標変換後の値の表示 */

  WaitForEnter();
  /* 読んだ要素から軌道計算 */
  State<double> ast_state_crtbp = init_ast_state_crtbp;  // 読んだ要素をここにコピー
  std::vector<std::array<double, 8>> history;
  TrajectoryObserver<double> observer(history, 0.0, kMU);
  auto trajectory_integrator = [&](const State<double>& state_ptr, double time,
                                   double h) -> State<double> {
    return SymplecticStep6thOrder(kMU, state_ptr, h);
  };

  Integrate(ast_state_crtbp, trajectory_integrator, observer, 0.0, 0.01, 1000);

  /* 軌道上の各点でSALI計算 */
  /* どの向きに速度ベクトルを与えるのがいいか調べる */
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>           SALI calculation on the trajectory" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n" << std::endl;
  double jacobi_range_min;
  double SALI_range_max;
  double SALI_range_step;
  std::vector<double> sali_initial_range;
  int SALI_point_interval;
  int SALI_calc_num = history.size() / SALI_point_interval;
  // 一番外側がヤコビ積分を変化させる次元、
  // その内側が時刻歴の次元
  // さらにその内側が各時刻のデータ配列
  std::vector<std::vector<std::array<double, 8>>> sali_history_list;
  auto sali_integrator = [&](SaliState<double>* state_ptr, double h) {
    SymplecticStep6thOrderSALI(kMU, state_ptr, h);
  };

  for
    /* 呼んだ要素から軌道計算 */

    std::cout << "<>  [mode selection] : " << std::endl;
  std::cout << "<> " << std::endl;
  std::cout << "<>        1. New simulation" << std::endl;
  std::cout << "<>        2. Detailed simulation for existing data" << std::endl;
  std::cout << "<>        else. Exit" << std::endl;
  std::cout << "<>   enter number " << std::endl;
  std::cout << "<> >>>";
  char mode;
  std::cin >> mode;
  std::cout << "<> " << std::endl;
  if (mode == '1') {
    std::cout << "<> > selected mode : New simulation" << std::endl;
  } else if (mode == '2') {
    std::cout << "<> > selected mode : Detailed simulation for existing data\n" << std::endl;
  } else {
    std::cout << "<> > selected mode : Exit\n" << std::endl;
    return 0;
  }
  std::cout << "<> " << std::endl;

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
  std::cout << "<> " << std::endl;
  std::cout << "<>        1. single simulation" << std::endl;
  std::cout << "<>        2. continuous simulation" << std::endl;
  std::cout << "<>        else. Exit" << std::endl;
  std::cout << "<>      enter number " << std::endl;
  std::cout << "<> >>> ";
  std::cin >> mode2;
  std::cout << "<> " << std::endl;

  double is_continuous = 0;
  std::string configfilename;
  std::vector<std::string> config_file_list;

  if (mode2 == '1') {
    std::cout << "<> >  selected mode : single simulation" << std::endl;
    is_continuous = 0;
  } else if (mode2 == '2') {
    std::cout << "<>    selected mode : continuous simulation\n" << std::endl;
    is_continuous = 1;
  } else {
    std::cout << "selected mode : Exit\n" << std::endl;
    return 0;
  }
  const std::string calc_config_pattern_str =
      "^" + calc_config_prefix + (is_continuous ? "_\\d+\\" : "") + ".txt$";
  const std::regex pattern("^" + calc_config_pattern_str);
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
  std::cout << "<> Loaded config file list:" <<　std::endl;
  for (const auto& filename : config_file_list) {
    std::cout << "<>    - " << filename << std::endl;
  }

  WaitForEnter();
  std::cout << "<> " << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  int Core_Max = omp_get_max_threads();
  int OMP_Fmax{};
  std::cout << "<>  [OpenMP preparation]" << std::endl;
  std::cout << "<> " << std::endl;
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
    std::cout << "<>" << std::endl;
    std::cout << "<>     >> Number of OMP threads is " << OMP_Fmax << std::endl;
    std::cout << "<>" << std::endl;
  } else {
    OMP_Fmax = Core_Max;
    std::cout << "  <>     >> Your input is INVALID. OMP threads is "
              << "automatically determined as " << OMP_Fmax << std::endl;
  }
  WaitForEnter();
  omp_set_num_threads(OMP_Fmax);
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  std::ifstream ifs;
  // 積分器
  auto integrator = [&](SaliState<double>* state_ptr, double h) {
    SymplecticStep4thOrderSALI(kMU, state_ptr, h);
  };
  int configdata_num = 0;
  //  実行時間の計測
  auto start_ofall = std::chrono::system_clock::now();

  // --------  configファイルの数だけSALI計算全体を繰り返す-----------------------
  // while (ifs) {
  for (const auto& configfilepath : config_file_list) {
    std::cout << "<>        Next config file : " << configfilepath << std::endl;
    configdata_num++;
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
    //   // 設定ファイル読み込み
    int MESH_CENTER = 0;
    State3d<int> MESH_DIVISION = {50, 50, 50};
    State3d<double> MESH_HALF_WIDTH{0.01, 0.01, 0.01};
    double CALC_TIMESTEP = 0;
    double SALI_CALCTIME_THRESHOLD = 0;
    double SOI_RADIUS = 0;
    double FOREBIDDEN_AREA_RADIUS = 0;
    double JACOBI_INTEGRAL = 0;
    double inclination = 0;
    double OMEGA = 0;
    double THETA = 0;
    while (std::getline(ifs, str)) {
      if (str.find("MESH CENTER") != std::string::npos) {
        MESH_CENTER = std::stoi(str.substr(str.find("=") + 1));
        if (MESH_CENTER == static_cast<int>(MeshCenter::kSUN)) {
          MeshCenter.x = -kMU;
          MeshCenter.y = 0.0;
          MeshCenter.z = 0.0;
        } else if (MESH_CENTER == static_cast<int>(MeshCenter::kEARTH)) {
          MeshCenter.x = 1.0 - kMU;
          MeshCenter.y = 0.0;
          MeshCenter.z = 0.0;
        } else {
          std::cerr << "<> !err! MESH CENTER is INVALID. EXITING..." << std::endl;
          return -1;
        }
      } else if (str.find("MESH DIVISION") != std::string::npos) {
        MESH_DIVISION.x = std::stod(str.substr(str.find("=") + 1));
        size_t first_space = str.find(" ", str.find("=") + 1);
        size_t second_space = str.find(" ", first_space + 1);
        MESH_DIVISION.y = std::stod(str.substr(first_space + 1, second_space - first_space - 1));
        MESH_DIVISION.z = std::stod(str.substr(second_space + 1));
      } else if (str.find("MESH SIZE") != std::string::npos) {
        MESH_HALF_WIDTH.x = std::stod(str.substr(str.find("=") + 1));
        size_t first_space = str.find(" ", str.find("=") + 1);
        size_t second_space = str.find(" ", first_space + 1);
        MESH_HALF_WIDTH.y = std::stod(str.substr(first_space + 1, second_space - first_space - 1));
        MESH_HALF_WIDTH.z = std::stod(str.substr(second_space + 1));
      } else if (str.find("CALC TIMESTEP") != std::string::npos) {
        CALC_TIMESTEP = std::stod(str.substr(str.find("=") + 1));
      } else if (str.find("SALI CALCTIME THRESHOLD") != std::string::npos) {
        SALI_CALCTIME_THRESHOLD = std::stod(str.substr(str.find("=") + 1));
      } else if (str.find("RADIUS OF SOI") != std::string::npos) {
        SOI_RADIUS = std::stod(str.substr(str.find("=") + 1));
      } else if (str.find("RADIUS OF FOREBIDDEN AREA") != std::string::npos) {
        FOREBIDDEN_AREA_RADIUS = std::stod(str.substr(str.find("=") + 1));
      } else if (str.find("JACOBI INTEGRAL") != std::string::npos) {
        JACOBI_INTEGRAL = std::stod(str.substr(str.find("=") + 1));
      } else if (str.find("INCLINATION AGAINST XY PLANE(deg)") != std::string::npos) {
        inclination = std::stod(str.substr(str.find("=") + 1));
        inclination = inclination * std::acos(-1) / 180.;
      } else if (str.find("LONGTITUDE AGAINST X AXIS+(deg)") != std::string::npos) {
        OMEGA = std::stod(str.substr(str.find("=") + 1));
        OMEGA = OMEGA * std::acos(-1) / 180.;
      } else if (str.find("DEGREE FROM TANGENT") != std::string::npos) {
        THETA = std::stod(str.substr(str.find("=") + 1));
        THETA = THETA * std::acos(-1) / 180.;
      }
    }
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
    if (mode == '1') {
      std::cout << "based on SOI radius" << std::endl;
      // meshPoints = CreateCircleMesh(SOI_RADIUS, MESH_SIZE, MeshCenter);
      meshPoints = createDimensionlessCartesianMesh(MeshCenter, MESH_HALF_WIDTH, MESH_DIVISION);
    } else if (mode == '2') {
      // std::cout << "based on the specified point" << std::endl;
      // meshPoints = create_cube_mesh(ROI_length, MESH_SIZE, MeshCenter);
    }
    int countt = meshPoints.size();
    std::cout << "<>        " << countt << " mesh generated successfully" << std::endl;
    std::cout << "<>        Start calclation" << std::endl;

    // ---------出力ファイル設定---------
    std::string output_base_path = OUTPUT_DIR;
    // シミュレーション終了時刻が同じでもファイル名が被らないようにする
    std::string filename = output_base_path + "/3D_crtbp_SALI/configdata_" +
                           std::to_string(configdata_num) + "_" + getcurrent_date() + ".dat";
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
    ofs1 << "# MESH CENTER=" << center_str << std::endl;
    ofs1 << "# MESH DIVISION=" << MESH_DIVISION.x << " " << MESH_DIVISION.y << " "
         << MESH_DIVISION.z << std::endl;
    ofs1 << "# MESH HALF WIDTH=" << MESH_HALF_WIDTH.x << " " << MESH_HALF_WIDTH.y << " "
         << MESH_HALF_WIDTH.z << std::endl;
    ofs1 << "# CALCULATION TIMESTEP=" << CALC_TIMESTEP << std::endl;
    ofs1 << "# SIMULATION TIME=" << SALI_CALCTIME_THRESHOLD << std::endl;
    ofs1 << "# RADIUSofSOI=" << SOI_RADIUS << std::endl;
    ofs1 << "# FOREBIDDEN AREA RADIUS=" << FOREBIDDEN_AREA_RADIUS << std::endl;
    ofs1 << "# INITIAL JACOBI INTEGRAL=" << JACOBI_INTEGRAL << std::endl;
    ofs1 << "# INCLINATION AGAINST XY PLANE=" << inclination / std::acos(-1) * 180. << std::endl;
    ofs1 << "# LONGTITUDE AGAINST X AXIS=" << OMEGA / std::acos(-1) * 180. << std::endl;
    ofs1 << "# DEGREE FROM TANGENT(deg)=" << THETA / std::acos(-1) * 180. << std::endl;
    ofs1 << "Time,SALI,x,y,z,px,py,pz\n";

    // 計算のステップ数
    int num_step = static_cast<int>(SALI_CALCTIME_THRESHOLD / CALC_TIMESTEP);
    int totalIterations = meshPoints.size();
    // 進捗カウンタ
    int completed_count = 0;

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
        double final_sali = -1;
        double vx = 0.0, vy = 0.0, vz = 0.0;

        double v_abs = calc_v_abs(point, JACOBI_INTEGRAL, kMU);

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
          // 1. 積分s
          integrator(&sali_state, CALC_TIMESTEP);
          // 2. 正規化
          sali_state.w1.Normalize();
          sali_state.w2.Normalize();
        }

        State<double> final_state = ConvertToPhysical(sali_state.state);
        if (calc_r2(final_state.x, final_state.y, final_state.z, kMU) < SOI_RADIUS &&
            calc_r2(final_state.x, final_state.y, final_state.z, kMU) > FOREBIDDEN_AREA_RADIUS) {
          const double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
          const double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
          final_sali = std::min(norm_plus, norm_minus);
        }
        local_output_buffer << mesh_num << "," << final_sali << "," << point.x << "," << point.y
                            << "," << point.z << "," << vx << "," << vy << "," << vz << "\n";
        // 一定件数ごとにバッファをファイルに書き込む (排他制御)
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
    std::cout << "<>    Output File:" << filename << std::endl;
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
