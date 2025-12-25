

#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector3d.hpp>
#include <vector>

#include "rtbp.hpp"

struct SaliOutputRow {
  int mesh_num;
  double sali;
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  int collision;
  int escape;
  double lower_limit_reach_time;
};

bool ParseSaliDataLine(const std::string& line, SaliOutputRow* output_row);

bool WriteSaliOutputsSortedByValue(const std::string& input_filename,
                                   std::string* sorted_output_filename);

/**
 * @brief enumの概要説明
 */
enum class MeshCenter {
  kSUN = 0,
  kEARTH = 1,
};

// ChaosIndexType は utils.hpp で定義されているため、ローカル定義は不要

int main(int argc, char* argv[]) {
  using namespace param;
  using namespace crtbp;
  using namespace utils;
  namespace fs = std::filesystem;

  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv);
  bool is_continuous = args.is_continuous;
  bool skip_wait = args.skip_wait;
  std::string output_tag = args.output_tag;

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

  State3d<double> MeshCenter{1.0 - kMU, 0, 0};

  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP 3dSALI  Calculation ver3.0" << std::endl;
  std::cout << "<>-------------------------------------------------------------"
               "---\n\n"
            << std::endl;

  // 選択されたモードを表示
  if (is_continuous) {
    std::cout << "<>    selected mode : continuous simulation (--continuous)\n" << std::endl;
  } else {
    std::cout << "<> >  selected mode : single simulation (default)" << std::endl;
  }

  if (skip_wait) {
    std::cout << "<>    WaitForEnter : skipped (--no-wait)" << std::endl;
  }

  // 設定ファイルを検索（_sample付きは除外）
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = is_continuous;
  std::vector<std::string> config_file_list =
      DiscoverConfigFilesToml(kCalcConfigPath, calc_config_prefix, discovery_opts);

  if (config_file_list.empty()) {
    std::cerr << "<> No config files found in " << kCalcConfigPath << std::endl;
    return -1;
  }

  std::cout << "<> Loaded config file list:" << std::endl;
  for (const auto& filename : config_file_list) {
    std::cout << "<>    - " << filename << std::endl;
  }

  if (!skip_wait) {
    WaitForEnter();
  }
  std::cout << "<> " << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  int Core_Max = omp_get_max_threads();
  int OMP_Fmax = Core_Max - 1;
  if (!skip_wait) {
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
    }
    std::cout << "  <>     >> Your input is INVALID. OMP threads is "
              << "automatically determined as " << OMP_Fmax << std::endl;
    WaitForEnter();
  } else {
    std::cout << "  <>     >> OMP threads is automatically determined as " << OMP_Fmax << std::endl;
  }
  omp_set_num_threads(OMP_Fmax);
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  std::ifstream ifs;
  // SALI積分器
  auto sali_integrator = [&](SaliState<double>* state_ptr, double h) {
    SymplecticStep6thOrderSALI(kMU, state_ptr, h);
  };
  // GALI4積分器
  auto gali4_integrator = [&](GaliState<double, 4>* state_ptr, double h) {
    SymplecticStep6thOrderGALI(kMU, state_ptr, h);
  };
  // GALI6積分器
  auto gali6_integrator = [&](GaliState<double, 6>* state_ptr, double h) {
    SymplecticStep6thOrderGALI(kMU, state_ptr, h);
  };
  //  実行時間の計測
  auto start_ofall = std::chrono::system_clock::now();

  // 実行ごとのセッション出力ディレクトリを作成
  std::string output_base_path = OUTPUT_DIR;
  std::string session_dir_name = getcurrent_date();
  if (!output_tag.empty()) {
    session_dir_name += "_" + output_tag;
  }
  std::string session_output_dir = output_base_path + "/3D_crtbp_SALI_v3/" + session_dir_name;
  if (!fs::exists(session_output_dir)) {
    fs::create_directories(session_output_dir);
  }
  std::cout << "<>" << std::endl;
  std::cout << "<>    Session output directory: " << session_output_dir << std::endl;
  std::cout << "<>" << std::endl;

  // --------  configファイルの数だけSALI計算全体を繰り返す-----------------------
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
    //   // 設定ファイル読み込み
    // TomlConfigParserを使用した設定読み込み
    std::unique_ptr<utils::TomlConfigParser> config_ptr;
    try {
      config_ptr = std::make_unique<utils::TomlConfigParser>(configfilepath);
    } catch (const std::exception& e) {
      std::cerr << "\n<> !ERROR! Failed to parse config file: " << configfilepath << std::endl;
      std::cerr << "<> !ERROR! Reason: " << e.what() << std::endl;
      std::cerr << "<> !ERROR! Please check if all required values are properly set." << std::endl;
      std::cerr << "<> !ERROR! Common issues:" << std::endl;
      std::cerr << "<>         - Empty values (e.g., 'jacobi_integral = ')" << std::endl;
      std::cerr << "<>         - Missing quotes around string values" << std::endl;
      std::cerr << "<>         - Invalid number format" << std::endl;
      std::cerr << "<> Skipping this config file...\n" << std::endl;
      ifs.close();
      continue;
    }
    auto& config = *config_ptr;

    // 必須パラメータのバリデーション
    std::vector<std::string> validation_errors;

    // jacobi_integral は必須（特別な値-999.0をセンチネルとして使用）
    const double kSentinelValue = -999.0;
    double jacobi_check = config.GetDouble("boundary.jacobi_integral", kSentinelValue);
    if (jacobi_check == kSentinelValue || !config.HasKey("boundary.jacobi_integral")) {
      validation_errors.push_back("boundary.jacobi_integral is missing or invalid");
    }

    // 配列パラメータのチェック
    auto division_check = config.GetDoubleArray("mesh.division");
    if (division_check.size() < 3) {
      validation_errors.push_back("mesh.division must have 3 values [x, y, z]");
    }

    auto halfwidth_check = config.GetDoubleArray("mesh.half_width");
    if (halfwidth_check.size() < 3) {
      validation_errors.push_back("mesh.half_width must have 3 values [x, y, z]");
    }

    // バリデーションエラーがある場合
    if (!validation_errors.empty()) {
      std::cerr << "\n<> !ERROR! Config validation failed for: " << configfilepath << std::endl;
      for (const auto& err : validation_errors) {
        std::cerr << "<> !ERROR!   - " << err << std::endl;
      }
      std::cerr << "<> Skipping this config file...\n" << std::endl;
      ifs.close();
      continue;
    }

    int MESH_CENTER = config.GetInt("mesh.center", 1);
    if (MESH_CENTER == static_cast<int>(MeshCenter::kSUN)) {
      MeshCenter.x = -kMU;
      MeshCenter.y = 0.0;
      MeshCenter.z = 0.0;
    } else if (MESH_CENTER == static_cast<int>(MeshCenter::kEARTH)) {
      MeshCenter.x = 1.0 - kMU;
      MeshCenter.y = 0.0;
      MeshCenter.z = 0.0;
    }

    auto division_arr = config.GetDoubleArray("mesh.division");
    State3d<int> MESH_DIVISION{50, 50, 50};
    if (division_arr.size() >= 3) {
      MESH_DIVISION.x = static_cast<int>(division_arr[0]);
      MESH_DIVISION.y = static_cast<int>(division_arr[1]);
      MESH_DIVISION.z = static_cast<int>(division_arr[2]);
    }

    auto halfwidth_arr = config.GetDoubleArray("mesh.half_width");
    State3d<double> MESH_HALF_WIDTH{0.01, 0.01, 0.01};
    if (halfwidth_arr.size() >= 3) {
      MESH_HALF_WIDTH.x = halfwidth_arr[0];
      MESH_HALF_WIDTH.y = halfwidth_arr[1];
      MESH_HALF_WIDTH.z = halfwidth_arr[2];
    }

    double CALC_TIMESTEP = config.GetDouble("simulation.calc_timestep", 0.000001);
    double SALI_CALCTIME_THRESHOLD = config.GetDouble("simulation.time_threshold", 10.0);
    double SOI_RADIUS = config.GetDouble("boundary.soi_radius", 0.03);
    double FOREBIDDEN_AREA_RADIUS = config.GetDouble("boundary.forbidden_area_radius", 1e-10);
    double JACOBI_INTEGRAL = config.GetDouble("boundary.jacobi_integral", 3.0008);

    double inclination =
        config.GetDouble("orientation.inclination_deg", 0.0) * std::acos(-1) / 180.;
    double OMEGA = config.GetDouble("orientation.longitude_deg", 0.0) * std::acos(-1) / 180.;
    double THETA = config.GetDouble("orientation.tangent_deg", 0.0) * std::acos(-1) / 180.;

    double SALI_LOWER_LIMIT = config.GetDouble("chaos.sali_lower_limit", 1e-8);

    std::string chaos_str = config.GetString("chaos.index_type", "SALI");
    ChaosIndexType chaos_index_type = ChaosIndexType::SALI;
    int gali_k = 2;
    if (chaos_str == "SALI" || chaos_str == "sali") {
      chaos_index_type = ChaosIndexType::SALI;
      gali_k = 2;
    } else if (chaos_str == "GALI2" || chaos_str == "gali2") {
      chaos_index_type = ChaosIndexType::GALI;
      gali_k = 2;
    } else if (chaos_str == "GALI4" || chaos_str == "gali4") {
      chaos_index_type = ChaosIndexType::GALI;
      gali_k = 4;
    } else if (chaos_str == "GALI6" || chaos_str == "gali6") {
      chaos_index_type = ChaosIndexType::GALI;
      gali_k = 6;
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
    std::cout << "<>        SALI LOWER LIMIT : " << SALI_LOWER_LIMIT << std::endl;
    std::string chaos_index_str;
    switch (chaos_index_type) {
      case ChaosIndexType::SALI:
        chaos_index_str = "SALI";
        break;
      case ChaosIndexType::GALI:
        chaos_index_str = "GALI" + std::to_string(gali_k);
        break;
    }
    std::cout << "<>        CHAOS INDEX : " << chaos_index_str << std::endl;
    std::cout << "<>    config file read successfully\n" << std::endl;
    if (!is_continuous && !skip_wait) {
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

    // ---------出力ファイル設定---------
    // configファイル名からベース名を取得
    std::string config_basename = fs::path(configfilepath).stem().string();
    std::string filename = session_output_dir + "/" + config_basename + ".dat";
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
    ofs1 << "# SALI LOWER LIMIT=" << SALI_LOWER_LIMIT << std::endl;
    ofs1 << "Time,SALI,x,y,z,px,py,pz,collision,escape,lower_limit_reach_time\n";

    // 計算のステップ数
    int num_step = static_cast<int>(SALI_CALCTIME_THRESHOLD / CALC_TIMESTEP);
    int totalIterations = meshPoints.size();
    // 進捗カウンタ
    int completed_count = 0;
    // 進捗表示間隔（ループ外で計算してオーバーヘッド削減）
    int display_interval = std::max(totalIterations / 100, 1);
    // ファイル書き込み間隔（大きくしてcritical削減）
    constexpr int kWriteInterval = 1000;

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
        double sali = -1;
        double vx = 0.0, vy = 0.0, vz = 0.0;
        int collision_flag = 0;
        int escape_flag = 0;
        double lower_limit_reach_time = SALI_CALCTIME_THRESHOLD;
        int lower_limit_reach_flag = 0;

        double v_abs = calc_v_abs(point, JACOBI_INTEGRAL, kMU);

        if (v_abs > 0) {
          State3d<double> velocity = calc_velocity(point, v_abs, kMU, inclination, OMEGA, THETA);
          vx = velocity.x;
          vy = velocity.y;
          vz = velocity.z;
          velo_err = 0;
        } else {
          velo_err = 1;
          lower_limit_reach_time = 0.0;
        }

        if (velo_err == 0) {
          State<double> initial_state = {point.x, point.y, point.z, vx, vy, vz};
          CanonicalState<double> canonical_state = ConvertToCanonical(initial_state);

          // 使用するstateポインタ（衝突/脱出判定用）
          CanonicalState<double>* current_state_ptr = nullptr;

          // SALI/GALI状態の初期化
          SaliState<double> sali_state;
          GaliState<double, 4> gali4_state;
          GaliState<double, 6> gali6_state;

          if (chaos_index_type == ChaosIndexType::SALI ||
              (chaos_index_type == ChaosIndexType::GALI && gali_k == 2)) {
            sali_state.state = canonical_state;
            sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
            current_state_ptr = &sali_state.state;
          } else if (chaos_index_type == ChaosIndexType::GALI && gali_k == 4) {
            gali4_state.state = canonical_state;
            gali4_state.w[0] = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            gali4_state.w[1] = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
            gali4_state.w[2] = CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
            gali4_state.w[3] = CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
            current_state_ptr = &gali4_state.state;
          } else if (chaos_index_type == ChaosIndexType::GALI && gali_k == 6) {
            gali6_state.state = canonical_state;
            gali6_state.w[0] = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            gali6_state.w[1] = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
            gali6_state.w[2] = CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
            gali6_state.w[3] = CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
            gali6_state.w[4] = CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
            gali6_state.w[5] = CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
            current_state_ptr = &gali6_state.state;
          }

          // 積分ループ (オブザーバー無し)
          for (int step = 0; step < num_step; ++step) {
            if (chaos_index_type == ChaosIndexType::SALI ||
                (chaos_index_type == ChaosIndexType::GALI && gali_k == 2)) {
              sali_integrator(&sali_state, CALC_TIMESTEP);
              sali_state.w1.Normalize();
              sali_state.w2.Normalize();
              double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
              double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
              sali = std::min(norm_plus, norm_minus);
              current_state_ptr = &sali_state.state;
            } else if (chaos_index_type == ChaosIndexType::GALI && gali_k == 4) {
              gali4_integrator(&gali4_state, CALC_TIMESTEP);
              gali4_state.NormalizeDeviationVectors();
              sali = gali4_state.ComputeGALI();
              current_state_ptr = &gali4_state.state;
            } else if (chaos_index_type == ChaosIndexType::GALI && gali_k == 6) {
              gali6_integrator(&gali6_state, CALC_TIMESTEP);
              gali6_state.NormalizeDeviationVectors();
              sali = gali6_state.ComputeGALI();
              current_state_ptr = &gali6_state.state;
            }

            if ((sali < SALI_LOWER_LIMIT) && lower_limit_reach_flag == 0) {
              lower_limit_reach_time = step * CALC_TIMESTEP;
              lower_limit_reach_flag = 1;
            }
            // Check for collision or escape
            if (current_state_ptr != nullptr) {
              double r2 =
                  calc_r2(current_state_ptr->qx, current_state_ptr->qy, current_state_ptr->qz, kMU);
              if (r2 < FOREBIDDEN_AREA_RADIUS) {
                collision_flag = 1;
              }
              if (r2 > SOI_RADIUS) {
                escape_flag = 1;
              }
            }
          }
        }
        local_output_buffer << mesh_num << "," << sali << "," << point.x << "," << point.y << ","
                            << point.z << "," << vx << "," << vy << "," << vz << ","
                            << collision_flag << "," << escape_flag << "," << lower_limit_reach_time
                            << "\n";
        // 一定件数ごとにバッファをファイルに書き込む (排他制御)
        if (idx % kWriteInterval == 0) {
#pragma omp critical(file_write)
          {
            ofs1 << local_output_buffer.str();
            local_output_buffer.str("");  // バッファをクリア
            local_output_buffer.clear();  // ストリーム状態もクリア
          }
        }
#pragma omp atomic
        completed_count++;

        // スレッド0のみが進捗表示（criticalを回避）
        if (omp_get_thread_num() == 0) {
          if (completed_count % display_interval == 0 ||
              completed_count >= totalIterations - OMP_Fmax) {
            double current_progress = static_cast<double>(completed_count) / totalIterations;
            displayProgressBarThreadSafe(current_progress);
          }
        }
      }
      // ループ終了後に残りのバッファを必ず書き込む
#pragma omp critical
      {
        ofs1 << local_output_buffer.str();
      }
    }  // end of parallel region

    // 最終的な進捗バー表示
    displayProgressBarThreadSafe(1.0);
    std::cout << std::endl;

    ofs1.close();
    std::string sorted_output_filename;
    if (WriteSaliOutputsSortedByValue(filename, &sorted_output_filename)) {
      std::cout << "<>    SALI sorted (desc) file: " << sorted_output_filename << std::endl;
    } else {
      std::cout << "<>    SALI sorting was skipped due to parse error." << std::endl;
    }
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
    std::cout << "<>        Simulation for " << configfilepath << " finished" << std::endl;
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

bool ParseSaliDataLine(const std::string& line, SaliOutputRow* output_row) {
  if (output_row == nullptr) {
    return false;
  }

  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string token;
  while (std::getline(ss, token, ',')) {
    fields.push_back(token);
  }

  if (fields.size() != 11) {
    return false;
  }

  try {
    output_row->mesh_num = static_cast<int>(std::stod(fields[0]));
    output_row->sali = std::stod(fields[1]);
    output_row->x = std::stod(fields[2]);
    output_row->y = std::stod(fields[3]);
    output_row->z = std::stod(fields[4]);
    output_row->vx = std::stod(fields[5]);
    output_row->vy = std::stod(fields[6]);
    output_row->vz = std::stod(fields[7]);
    output_row->collision = std::stoi(fields[8]);
    output_row->escape = std::stoi(fields[9]);
    output_row->lower_limit_reach_time = std::stod(fields[10]);
  } catch (const std::exception&) {
    return false;
  }

  return true;
}

bool WriteSaliOutputsSortedByValue(const std::string& input_filename,
                                   std::string* sorted_output_filename) {
  if (sorted_output_filename == nullptr) {
    return false;
  }

  std::ifstream input(input_filename);
  if (!input.is_open()) {
    std::cerr << "Failed to open file for SALI sorting: " << input_filename << std::endl;
    return false;
  }

  std::vector<std::string> headers;
  std::vector<SaliOutputRow> rows;
  std::string line;
  bool header_finished = false;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    if (!header_finished && (line[0] == '#' || line.rfind("Time", 0) == 0)) {
      headers.push_back(line);
      if (line[0] != '#') {
        header_finished = true;
      }
      continue;
    }
    SaliOutputRow parsed;
    if (ParseSaliDataLine(line, &parsed)) {
      rows.push_back(parsed);
    }
  }

  if (rows.empty()) {
    std::cerr << "No SALI data lines to sort in " << input_filename << std::endl;
    return false;
  }

  std::sort(rows.begin(), rows.end(), [](const SaliOutputRow& lhs, const SaliOutputRow& rhs) {
    // まずlower_limit_reach_timeで昇順ソート（短い方が上位）
    if (lhs.lower_limit_reach_time != rhs.lower_limit_reach_time) {
      return lhs.lower_limit_reach_time < rhs.lower_limit_reach_time;
    }
    // lower_limit_reach_timeが同じ場合はsaliで昇順ソート（小さい方が上位）
    return lhs.sali < rhs.sali;
  });

  std::filesystem::path input_path(input_filename);
  std::filesystem::path sorted_path =
      input_path.parent_path() /
      (input_path.stem().string() + "_sorted" + input_path.extension().string());

  std::ofstream output(sorted_path);
  if (!output.is_open()) {
    std::cerr << "Failed to create SALI sorted file: " << sorted_path << std::endl;
    return false;
  }
  output << std::fixed << std::setprecision(15);
  for (const auto& header : headers) {
    output << header << '\n';
  }
  for (const auto& row : rows) {
    output << row.mesh_num << "," << row.sali << "," << row.x << "," << row.y << "," << row.z << ","
           << row.vx << "," << row.vy << "," << row.vz << "," << row.collision << "," << row.escape
           << "," << row.lower_limit_reach_time << "\n";
  }

  *sorted_output_filename = sorted_path.string();
  return true;
}
