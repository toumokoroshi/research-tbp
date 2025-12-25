/**
 * @file main.cpp
 * @brief 円制限三体問題のSALI計算アプリケーション (J2000→回転座標変換対応)
 *
 * J2000系の小惑星・地球状態量を読み込み、CR3BP回転座標に変換後、
 * 速度減速を加えてSALI計算を行う。
 *
 * @note Google C++ Style Guide準拠、Doxygenフォーマット
 */

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <rtbp.hpp>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector3d.hpp>
#include <vector>

namespace fs = std::filesystem;

using namespace my_type;
using namespace crtbp;
using namespace param;
using namespace utils;

/**
 * @brief 積分器の種類
 */
enum class IntegratorType {
  kSymplectic6th,  ///< 6次吉田法（デフォルト）
  kSymplectic4th,  ///< 4次吉田法
};

/**
 * @brief カオス指標の種類
 */
enum class ChaosIndexType {
  SALI,  ///< SALI (K=2) - デフォルト
  GALI   ///< GALI (K可変)
};

/**
 * @brief 設定パラメータ構造体
 */
struct TrajectorySaliConfig {
  std::string data_dir;            ///< 軌道データディレクトリ
  double calc_duration_nd = 10.0;  ///< 無次元計算時間
  double calc_timestep_nd = 0.01;  ///< 無次元タイムステップ
  double forbidden_radius = 1e-6;  ///< 衝突判定半径
  double escape_judge_hill = 1.0;  ///< エスケープ判定倍率
  double dv_max = 0.01;            ///< 最大減速量（無次元）
  double dv_step = 0.002;          ///< 減速ステップ（無次元）
  int phase_step = 1;              ///< 軌道位相の間引き（1=全て使用）
  IntegratorType integrator = IntegratorType::kSymplectic6th;
  double sali_lower_limit = 1e-10;                         ///< SALI閾値
  bool output_sali_timeseries = false;                     ///< SALI時系列出力フラグ
  ChaosIndexType chaos_index_type = ChaosIndexType::SALI;  ///< カオス指標の種類
  int gali_k = 2;                                          ///< GALIの偏差ベクトル数 (2, 4, 6)
  bool output_converted_trajectory = false;                ///< 座標変換後軌道の時刻歴出力フラグ
  bool output_impulse_trajectory = false;  ///< インパルス付与後軌道の時刻歴出力フラグ
};

/**
 * @brief 軌道データ1行分
 */
struct OrbitDataRow {
  double time_j2000;       ///< J2000時刻
  State<double> asteroid;  ///< 小惑星状態量（AU, AU/day）
  State<double> earth;     ///< 地球状態量（AU, AU/day）
};

/**
 * @brief TOML設定ファイルを読み込む
 * @param[in] filepath 設定ファイルパス
 * @param[out] config 設定構造体
 * @return 成功時true
 */
bool LoadTrajectorySaliConfig(const std::string& filepath, TrajectorySaliConfig* config) {
  try {
    utils::TomlConfigParser toml_config(filepath);

    config->data_dir = toml_config.GetString("data.data_dir", "");
    config->calc_duration_nd = toml_config.GetDouble("simulation.calc_duration_nd", 10.0);
    config->calc_timestep_nd = toml_config.GetDouble("simulation.calc_timestep_nd", 0.01);
    config->forbidden_radius = toml_config.GetDouble("simulation.forbidden_radius", 1e-6);
    config->escape_judge_hill = toml_config.GetDouble("simulation.escape_judge_hill", 1.0);
    config->dv_max = toml_config.GetDouble("deceleration.dv_max", 0.01);
    config->dv_step = toml_config.GetDouble("deceleration.dv_step", 0.002);
    config->phase_step = toml_config.GetInt("deceleration.phase_step", 1);

    std::string integrator_str = toml_config.GetString("integrator.type", "SYMPLECTIC_6TH");
    if (integrator_str == "SYMPLECTIC_6TH") {
      config->integrator = IntegratorType::kSymplectic6th;
    } else if (integrator_str == "SYMPLECTIC_4TH") {
      config->integrator = IntegratorType::kSymplectic4th;
    }

    config->sali_lower_limit = toml_config.GetDouble("chaos.sali_lower_limit", 1e-8);
    config->output_sali_timeseries = toml_config.GetBool("chaos.output_sali_timeseries", false);

    std::string chaos_str = toml_config.GetString("chaos.index_type", "SALI");
    if (chaos_str == "SALI" || chaos_str == "sali") {
      config->chaos_index_type = ChaosIndexType::SALI;
      config->gali_k = 2;
    } else if (chaos_str == "GALI2" || chaos_str == "gali2") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 2;
    } else if (chaos_str == "GALI4" || chaos_str == "gali4") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 4;
    } else if (chaos_str == "GALI6" || chaos_str == "gali6") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 6;
    }

    // 出力オプション
    config->output_converted_trajectory =
        toml_config.GetBool("output.output_converted_trajectory", false);
    config->output_impulse_trajectory =
        toml_config.GetBool("output.output_impulse_trajectory", false);

  } catch (const std::exception& e) {
    std::cerr << "<> !err! Cannot load config: " << e.what() << std::endl;
    return false;
  }
  return true;
}

/**
 * @brief 積分器名を文字列で返す
 */
std::string GetIntegratorName(IntegratorType type) {
  switch (type) {
    case IntegratorType::kSymplectic6th:
      return "Symplectic 6th Order (Yoshida)";
    case IntegratorType::kSymplectic4th:
      return "Symplectic 4th Order (Yoshida)";
    default:
      return "Unknown";
  }
}

/**
 * @brief 軌道データファイルを読み込む
 * @param[in] filepath ファイルパス
 * @param[in] phase_step 間引きステップ
 * @param[out] data 読み込んだデータ
 * @return 成功時true
 */
bool LoadOrbitData(const std::string& filepath, int phase_step, std::vector<OrbitDataRow>* data) {
  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "<> !err! Cannot open orbit data file: " << filepath << std::endl;
    return false;
  }

  std::string line;
  bool data_started = false;
  int line_count = 0;

  while (std::getline(ifs, line)) {
    // 空行でデータ開始
    if (!data_started) {
      if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
        data_started = true;
      }
      continue;
    }

    // 空行スキップ
    if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
      continue;
    }

    // 間引き
    if (line_count % phase_step != 0) {
      line_count++;
      continue;
    }
    line_count++;

    // データパース: 時刻 小惑星(6) 地球(6) = 13個
    std::stringstream ss(line);
    OrbitDataRow row;
    ss >> row.time_j2000;
    ss >> row.asteroid.x >> row.asteroid.y >> row.asteroid.z;
    ss >> row.asteroid.vx >> row.asteroid.vy >> row.asteroid.vz;
    ss >> row.earth.x >> row.earth.y >> row.earth.z;
    ss >> row.earth.vx >> row.earth.vy >> row.earth.vz;

    if (!ss.fail()) {
      data->push_back(row);
    }
  }
  return !data->empty();
}

/**
 * @brief DATA_DIRからEarth軌道ファイルを発見
 */
std::vector<std::string> DiscoverOrbitFiles(const std::string& data_dir) {
  std::vector<std::string> files;
  if (!fs::exists(data_dir)) {
    std::cerr << "<> !err! DATA_DIR does not exist: " << data_dir << std::endl;
    return files;
  }

  const std::regex pattern(R"(OBT_\d+_Earth\.txt)");
  for (const auto& entry : fs::directory_iterator(data_dir)) {
    if (!entry.is_regular_file()) continue;
    std::string filename = entry.path().filename().string();
    if (std::regex_match(filename, pattern)) {
      files.push_back(entry.path().string());
    }
  }

  std::sort(files.begin(), files.end());
  return files;
}

/**
 * @brief gnuplotでカラーコンタ図を生成（sali_final, sali_threshold_time）
 * @param[in] csv_path 入力CSVファイルパス
 * @param[in] output_dir 出力ディレクトリ
 * @param[in] file_basename ファイル名のベース
 */
void GenerateSaliContourPlots(const std::string& csv_path, const std::string& output_dir,
                              const std::string& file_basename) {
  // gnuplotスクリプトファイル
  std::string gp_script_path = output_dir + "/plot_sali_contour.gp";
  std::ofstream gp(gp_script_path);
  if (!gp) {
    std::cerr << "<> !warn! Cannot create gnuplot script: " << gp_script_path << std::endl;
    return;
  }

  // sali_final用プロット
  std::string eps_sali = output_dir + "/" + file_basename + "_sali_final.eps";
  std::string png_sali = output_dir + "/" + file_basename + "_sali_final.png";

  // sali_threshold_time用プロット
  std::string eps_time = output_dir + "/" + file_basename + "_sali_threshold_time.eps";
  std::string png_time = output_dir + "/" + file_basename + "_sali_threshold_time.png";

  gp << "set datafile separator ','\n";
  gp << "set datafile commentschars '#'\n";
  gp << "\n";

  // sali_final カラーコンタ (EPS)
  gp << "set terminal postscript eps enhanced color font 'Arial,14'\n";
  gp << "set output '" << eps_sali << "'\n";
  gp << "set xlabel 'Phase Index'\n";
  gp << "set ylabel 'Delta-V (non-dim)'\n";
  gp << "set title 'SALI Final Value - " << file_basename << "'\n";
  gp << "set cblabel 'log_{10}(SALI)'\n";
  gp << "set palette defined (0 'dark-blue', 1 'blue', 2 'cyan', 3 'green', 4 'yellow', 5 "
        "'orange', 6 'red')\n";
  gp << "set pm3d map corners2color c1\n";
  gp << "set dgrid3d 50,50\n";
  gp << "splot '" << csv_path << "' using 1:3:(log10($4 > 0 ? $4 : 1e-16)) notitle\n";
  gp << "\n";

  // sali_final カラーコンタ (PNG)
  gp << "set terminal pngcairo enhanced font 'Arial,14' size 1200,900\n";
  gp << "set output '" << png_sali << "'\n";
  gp << "replot\n";
  gp << "\n";

  // sali_threshold_time カラーコンタ (EPS)
  gp << "set terminal postscript eps enhanced color font 'Arial,14'\n";
  gp << "set output '" << eps_time << "'\n";
  gp << "set xlabel 'Phase Index'\n";
  gp << "set ylabel 'Delta-V (non-dim)'\n";
  gp << "set title 'SALI Threshold Time - " << file_basename << "'\n";
  gp << "set cblabel 'Threshold Time (non-dim)'\n";
  gp << "set palette defined (0 'dark-blue', 1 'blue', 2 'cyan', 3 'green', 4 'yellow', 5 "
        "'orange', 6 'red')\n";
  gp << "set pm3d map corners2color c1\n";
  gp << "splot '" << csv_path << "' using 1:3:5 notitle\n";
  gp << "\n";

  // sali_threshold_time カラーコンタ (PNG)
  gp << "set terminal pngcairo enhanced font 'Arial,14' size 1200,900\n";
  gp << "set output '" << png_time << "'\n";
  gp << "replot\n";
  gp << "\n";

  gp.close();

  // gnuplot実行
  std::string cmd = "gnuplot \"" + gp_script_path + "\"";
  std::cout << "<>    Generating plots..." << std::endl;
  int ret = std::system(cmd.c_str());
  if (ret != 0) {
    std::cerr << "<> !warn! gnuplot returned non-zero: " << ret << std::endl;
  } else {
    std::cout << "<>    Plots generated:" << std::endl;
    std::cout << "<>      - " << eps_sali << std::endl;
    std::cout << "<>      - " << png_sali << std::endl;
    std::cout << "<>      - " << eps_time << std::endl;
    std::cout << "<>      - " << png_time << std::endl;
  }
}

/**
 * @brief gnuplotでSALI時系列プロットを生成
 * @param[in] csv_path 入力CSVファイルパス（時系列データ）
 * @param[in] output_dir 出力ディレクトリ
 * @param[in] file_basename ファイル名のベース
 * @param[in] phase_idx 位相インデックス
 * @param[in] dv_mag 速度変化量
 */
void GenerateSaliTimeSeriesPlot(const std::string& csv_path, const std::string& output_dir,
                                const std::string& file_basename, int phase_idx, double dv_mag) {
  // gnuplotスクリプトファイル
  std::ostringstream script_name;
  script_name << "/plot_sali_ts_p" << phase_idx << "_dv" << std::fixed << std::setprecision(4)
              << dv_mag << ".gp";
  std::string gp_script_path = output_dir + script_name.str();
  std::ofstream gp(gp_script_path);
  if (!gp) {
    return;  // スクリプト作成失敗時は無視
  }

  // 出力ファイル名生成
  std::ostringstream base_name;
  base_name << file_basename << "_sali_ts_p" << phase_idx << "_dv" << std::fixed
            << std::setprecision(4) << dv_mag;
  std::string eps_path = output_dir + "/" + base_name.str() + ".eps";
  std::string png_path = output_dir + "/" + base_name.str() + ".png";

  gp << "set datafile separator ','\n";
  gp << "set datafile commentschars '#'\n";
  gp << "\n";

  // EPS出力
  gp << "set terminal postscript eps enhanced color font 'Arial,14'\n";
  gp << "set output '" << eps_path << "'\n";
  gp << "set xlabel 'Time (non-dim)'\n";
  gp << "set ylabel 'SALI'\n";
  gp << "set title 'SALI Time Series - Phase " << phase_idx << ", dv=" << std::fixed
     << std::setprecision(4) << dv_mag << "'\n";
  gp << "set grid\n";
  gp << "set logscale y\n";
  gp << "set format y '10^{%L}'\n";
  gp << "plot '" << csv_path
     << "' using 1:($2 > 0 ? $2 : 1e-16) with lines lw 2 lc rgb 'blue' title 'SALI'\n";
  gp << "\n";

  // PNG出力
  gp << "set terminal pngcairo enhanced font 'Arial,14' size 1200,600\n";
  gp << "set output '" << png_path << "'\n";
  gp << "replot\n";

  gp.close();

  // gnuplot実行
  std::string cmd = "gnuplot \"" + gp_script_path + "\"";
  std::system(cmd.c_str());
}

int main(int argc, char* argv[]) {
  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv);
  bool skip_wait = args.skip_wait;
  std::string output_tag = args.output_tag;

  // ヘッダー出力
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP Trajectory-SALI Calculator ver1.0" << std::endl;
  std::cout << "<>            (J2000 -> Rotating Frame Conversion)" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // パス設定
  const std::string kConfigBasePath = CONFIG_DIR;
  const std::string kOutputBasePath = OUTPUT_DIR;
  const std::string kConfigFilePath =
      kConfigBasePath + "/trajectory_SALI/trajectory_sali_config_sample.toml";
  const std::string kAstroParamFile = kConfigBasePath + "/astro_param/astro_param.txt";

  // 天文定数読み込み
  AstroConstants<double> astro_params = loadConstants<double>(kAstroParamFile);
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  std::cout << "<>    mu parameter: " << std::setprecision(15) << kMU << std::endl;
  std::cout << "<>" << std::endl;

  // 設定ファイル読み込み
  std::cout << "<>    Loading config file: " << kConfigFilePath << std::endl;
  TrajectorySaliConfig config;
  if (!LoadTrajectorySaliConfig(kConfigFilePath, &config)) {
    return -1;
  }

  // 軌道データファイル発見
  std::vector<std::string> orbit_files = DiscoverOrbitFiles(config.data_dir);
  if (orbit_files.empty()) {
    std::cerr << "<> !err! No orbit files found in: " << config.data_dir << std::endl;
    return -1;
  }

  // 設定確認表示
  std::cout << "<>" << std::endl;
  std::cout << "<>****************************************************************" << std::endl;
  std::cout << "<>  [Configuration Summary]" << std::endl;
  std::cout << "<>" << std::endl;
  std::cout << "<>    DATA_DIR: " << config.data_dir << std::endl;
  std::cout << "<>    Orbit files found: " << orbit_files.size() << std::endl;
  std::cout << "<>    CALC_DURATION_ND: " << config.calc_duration_nd << std::endl;
  std::cout << "<>    CALC_TIMESTEP_ND: " << config.calc_timestep_nd << std::endl;
  std::cout << "<>    FORBIDDEN_RADIUS: " << config.forbidden_radius << std::endl;
  std::cout << "<>    ESCAPE_JUDGE_HILL: " << config.escape_judge_hill << std::endl;
  std::cout << "<>    DV_MAX: " << config.dv_max << std::endl;
  std::cout << "<>    DV_STEP: " << config.dv_step << std::endl;
  std::cout << "<>    PHASE_STEP: " << config.phase_step << std::endl;
  std::cout << "<>    INTEGRATOR: " << GetIntegratorName(config.integrator) << std::endl;
  std::cout << "<>    SALI_LOWER_LIMIT: " << config.sali_lower_limit << std::endl;
  std::cout << "<>    OUTPUT_SALI_TIMESERIES: " << (config.output_sali_timeseries ? "ON" : "OFF")
            << std::endl;
  std::string chaos_index_str;
  switch (config.chaos_index_type) {
    case ChaosIndexType::SALI:
      chaos_index_str = "SALI";
      break;
    case ChaosIndexType::GALI:
      chaos_index_str = "GALI" + std::to_string(config.gali_k);
      break;
  }
  std::cout << "<>    CHAOS_INDEX: " << chaos_index_str << std::endl;
  std::cout << "<>    OUTPUT_CONVERTED_TRAJ: "
            << (config.output_converted_trajectory ? "ON" : "OFF") << std::endl;
  std::cout << "<>    OUTPUT_IMPULSE_TRAJ: " << (config.output_impulse_trajectory ? "ON" : "OFF")
            << std::endl;
  std::cout << "<>" << std::endl;
  std::cout << "<>    Velocity deceleration steps: "
            << static_cast<int>(config.dv_max / config.dv_step) + 1 << std::endl;
  std::cout << "<>" << std::endl;
  std::cout << "<>****************************************************************" << std::endl;
  std::cout << "<>" << std::endl;

  // ユーザー確認
  int core_max = omp_get_max_threads();
  if (!skip_wait) {
    WaitForEnter("<>    Press Enter to start simulation...");

    // OpenMP設定
    int omp_threads = 1;
    std::cout << "<>" << std::endl;
    std::cout << "<>  [OpenMP Configuration]" << std::endl;
    std::cout << "<>    Available threads: " << core_max << std::endl;
    std::cout << "<>    Enter number of threads to use: ";
    std::cin >> omp_threads;
    if (omp_threads > core_max) omp_threads = core_max;
    if (omp_threads < 1) omp_threads = 1;
    omp_set_num_threads(omp_threads);
    std::cout << "<>    Using " << omp_threads << " threads" << std::endl;
    std::cout << "<>" << std::endl;
  } else {
    // skip_wait時は最大スレッド数-1を自動選択
    int omp_threads = std::max(1, core_max - 1);
    omp_set_num_threads(omp_threads);
    std::cout << "<>    Using " << omp_threads << " threads (auto-selected)" << std::endl;
  }

  // SALI積分器ラムダ
  auto sali_integrator = [&](SaliState<double>* state, double h) {
    switch (config.integrator) {
      case IntegratorType::kSymplectic6th:
        SymplecticStep6thOrderSALI(kMU, state, h);
        break;
      case IntegratorType::kSymplectic4th:
        SymplecticStep4thOrderSALI(kMU, state, h);
        break;
    }
  };

  // GALI4積分器ラムダ
  auto gali4_integrator = [&](GaliState<double, 4>* state, double h) {
    switch (config.integrator) {
      case IntegratorType::kSymplectic6th:
        SymplecticStep6thOrderGALI(kMU, state, h);
        break;
      case IntegratorType::kSymplectic4th:
        SymplecticStep4thOrderGALI(kMU, state, h);
        break;
    }
  };

  // GALI6積分器ラムダ
  auto gali6_integrator = [&](GaliState<double, 6>* state, double h) {
    switch (config.integrator) {
      case IntegratorType::kSymplectic6th:
        SymplecticStep6thOrderGALI(kMU, state, h);
        break;
      case IntegratorType::kSymplectic4th:
        SymplecticStep4thOrderGALI(kMU, state, h);
        break;
    }
  };

  // ベース出力ディレクトリ
  std::string output_base_dir = kOutputBasePath + "/trajectory_SALI";
  if (!fs::exists(output_base_dir)) {
    fs::create_directories(output_base_dir);
  }

  // 全体時間計測
  auto start_total = std::chrono::system_clock::now();

  // シミュレーション開始時刻のタイムスタンプを生成（全軌道ファイル共通）
  std::time_t start_time_t = std::chrono::system_clock::to_time_t(start_total);
  std::tm start_local_tm;
#ifdef _WIN32
  localtime_s(&start_local_tm, &start_time_t);
#else
  localtime_r(&start_time_t, &start_local_tm);
#endif
  std::ostringstream sim_timestamp_ss;
  sim_timestamp_ss << std::put_time(&start_local_tm, "%y_%m%d_%H%M") << "_run";
  if (!output_tag.empty()) {
    sim_timestamp_ss << "_" << output_tag;
  }
  std::string sim_output_dir = output_base_dir + "/" + sim_timestamp_ss.str();
  fs::create_directories(sim_output_dir);
  std::cout << "<>    Simulation output folder: " << sim_output_dir << std::endl;
  std::cout << "<>" << std::endl;

  // 各軌道ファイルをループ
  for (size_t file_idx = 0; file_idx < orbit_files.size(); ++file_idx) {
    const std::string& orbit_file = orbit_files[file_idx];
    std::string file_basename = fs::path(orbit_file).stem().string();

    // 軌道ファイルごとのサブフォルダをシミュレーションフォルダ内に作成
    std::string output_dir = sim_output_dir + "/" + file_basename;
    fs::create_directories(output_dir);

    // SALI時系列データ用サブディレクトリ
    std::string timeseries_dir = output_dir + "/sali_timeseries";
    fs::create_directories(timeseries_dir);

    // インパルス付与後軌道データ用サブディレクトリ
    std::string impulse_traj_dir = output_dir + "/impulse_trajectory";
    if (config.output_impulse_trajectory) {
      fs::create_directories(impulse_traj_dir);
    }

    // 座標変換後軌道データ用CSVファイル（全位相点の変換結果を並べたもの）
    std::string converted_csv_path = output_dir + "/" + file_basename + "_converted_states.csv";
    std::ofstream converted_ofs;
    if (config.output_converted_trajectory) {
      converted_ofs.open(converted_csv_path);
      if (converted_ofs) {
        converted_ofs << std::setprecision(15) << std::fixed;
        converted_ofs << "# Converted States (J2000 -> Rotating Frame)\n";
        converted_ofs << "# Source: " << orbit_file << "\n";
        converted_ofs << "phase_idx,time_j2000,x,y,z,vx,vy,vz,jacobi\n";
      }
    }

    std::cout << "<>----------------------------------------------------------------" << std::endl;
    std::cout << "<>    Processing file " << (file_idx + 1) << " / " << orbit_files.size() << ": "
              << file_basename << std::endl;
    std::cout << "<>    Output folder: " << output_dir << std::endl;

    // 軌道データ読み込み
    std::vector<OrbitDataRow> orbit_data;
    if (!LoadOrbitData(orbit_file, config.phase_step, &orbit_data)) {
      std::cerr << "<> !warn! Skipping file due to load error" << std::endl;
      continue;
    }
    std::cout << "<>    Loaded " << orbit_data.size() << " phase points" << std::endl;

    // 出力ファイル
    std::string csv_path = output_dir + "/" + file_basename + "_sali_map.csv";
    std::ofstream ofs(csv_path);
    if (!ofs) {
      std::cerr << "<> !err! Cannot create output file: " << csv_path << std::endl;
      continue;
    }

    // CSVヘッダー
    ofs << std::setprecision(15) << std::fixed;
    ofs << "# TrajectorySALI Output\n";
    ofs << "# Source: " << orbit_file << "\n";
    ofs << "# CALC_DURATION_ND=" << config.calc_duration_nd << "\n";
    ofs << "# CALC_TIMESTEP_ND=" << config.calc_timestep_nd << "\n";
    ofs << "# DV_MAX=" << config.dv_max << "\n";
    ofs << "# DV_STEP=" << config.dv_step << "\n";
    ofs << "# SALI_LOWER_LIMIT=" << config.sali_lower_limit << "\n";
    ofs << "phase_idx,time_j2000,dv_mag,sali_final,sali_threshold_time,"
        << "escape,collision,x0,y0,z0,vx0,vy0,vz0,jacobi\n";

    // 計算ステップ数
    int num_steps = static_cast<int>(config.calc_duration_nd / config.calc_timestep_nd);
    int dv_count = static_cast<int>(config.dv_max / config.dv_step) + 1;
    int total_calcs = orbit_data.size() * dv_count;
    int completed = 0;

    auto start_file = std::chrono::system_clock::now();

// 各位相点×各速度変化量でループ
#pragma omp parallel for schedule(dynamic) shared(completed, ofs)
    for (int phase_idx = 0; phase_idx < static_cast<int>(orbit_data.size()); ++phase_idx) {
      const OrbitDataRow& row = orbit_data[phase_idx];

      // J2000→CR3BP回転座標変換
      State<double> asteroid_rot =
          ConvertInertial2RotatingV2(row.asteroid, row.earth, astro_params);

      // 速度変化ループ
      for (int dv_idx = 0; dv_idx < dv_count; ++dv_idx) {
        double dv_mag = dv_idx * config.dv_step;

        // 速度逆方向に減速
        double v_mag =
            std::sqrt(asteroid_rot.vx * asteroid_rot.vx + asteroid_rot.vy * asteroid_rot.vy +
                      asteroid_rot.vz * asteroid_rot.vz);
        double dv_x = 0.0, dv_y = 0.0, dv_z = 0.0;
        if (v_mag > 1e-15) {
          dv_x = -dv_mag * asteroid_rot.vx / v_mag;
          dv_y = -dv_mag * asteroid_rot.vy / v_mag;
          dv_z = -dv_mag * asteroid_rot.vz / v_mag;
        }

        // 減速適用後の状態
        State<double> modified_state = asteroid_rot;
        modified_state.vx += dv_x;
        modified_state.vy += dv_y;
        modified_state.vz += dv_z;

        // カオス指標計算用の状態初期化
        SaliState<double> sali_state;
        GaliState<double, 4> gali4_state;
        GaliState<double, 6> gali6_state;
        CanonicalState<double> canonical_state = ConvertToCanonical(modified_state);

        if (config.chaos_index_type == ChaosIndexType::SALI ||
            (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 2)) {
          sali_state.state = canonical_state;
          sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
          sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        } else if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 4) {
          gali4_state.state = canonical_state;
          gali4_state.w[0] = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
          gali4_state.w[1] = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
          gali4_state.w[2] = CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
          gali4_state.w[3] = CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        } else if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 6) {
          gali6_state.state = canonical_state;
          gali6_state.w[0] = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
          gali6_state.w[1] = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
          gali6_state.w[2] = CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
          gali6_state.w[3] = CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
          gali6_state.w[4] = CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
          gali6_state.w[5] = CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
        }

        double chaos_final = -1.0;
        double chaos_threshold_time = config.calc_duration_nd;
        int escape_flag = 0;
        int collision_flag = 0;
        bool threshold_reached = false;

        // 初期ヤコビ積分
        double jacobi = calc_jacobi_integral(modified_state, kMU);

        // カオス指標時系列データを記録するベクトル
        std::vector<std::pair<double, double>> chaos_timeseries;
        chaos_timeseries.reserve(num_steps);

        // 軌道時刻歴データを記録するベクトル (time, x, y, z, vx, vy, vz, jacobi)
        // インパルス付与後軌道用（時系列出力が必要な場合のみ）
        std::vector<std::tuple<double, double, double, double, double, double, double, double>>
            trajectory_timeseries;
        if (config.output_impulse_trajectory && dv_idx > 0) {
          trajectory_timeseries.reserve(num_steps + 1);
          // 初期状態を記録
          trajectory_timeseries.emplace_back(0.0, modified_state.x, modified_state.y,
                                             modified_state.z, modified_state.vx, modified_state.vy,
                                             modified_state.vz, jacobi);
        }

        // 積分ループ
        for (int step = 0; step < num_steps; ++step) {
          CanonicalState<double>* current_state = nullptr;

          if (config.chaos_index_type == ChaosIndexType::SALI ||
              (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 2)) {
            sali_integrator(&sali_state, config.calc_timestep_nd);
            sali_state.w1.Normalize();
            sali_state.w2.Normalize();
            double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
            double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
            chaos_final = std::min(norm_plus, norm_minus);
            current_state = &sali_state.state;
          } else if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 4) {
            gali4_integrator(&gali4_state, config.calc_timestep_nd);
            gali4_state.NormalizeDeviationVectors();
            chaos_final = gali4_state.ComputeGALI();
            current_state = &gali4_state.state;
          } else if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 6) {
            gali6_integrator(&gali6_state, config.calc_timestep_nd);
            gali6_state.NormalizeDeviationVectors();
            chaos_final = gali6_state.ComputeGALI();
            current_state = &gali6_state.state;
          }

          // 時系列データを記録
          double current_time = (step + 1) * config.calc_timestep_nd;
          chaos_timeseries.emplace_back(current_time, chaos_final);

          // 閾値到達判定
          if (!threshold_reached && chaos_final < config.sali_lower_limit) {
            chaos_threshold_time = current_time;
            threshold_reached = true;
          }

          // 衝突・エスケープ判定
          if (current_state != nullptr) {
            double r2 = calc_r2(current_state->qx, current_state->qy, current_state->qz, kMU);
            if (r2 < config.forbidden_radius) {
              collision_flag = 1;
            }
            if (r2 > config.escape_judge_hill) {
              escape_flag = 1;
            }

            // 軌道時刻歴データを収集（インパルス付与後軌道用）
            if (config.output_impulse_trajectory && dv_idx > 0) {
              State<double> current_cartesian = ConvertToPhysical(*current_state);
              double current_jacobi = calc_jacobi_integral(current_cartesian, kMU);
              trajectory_timeseries.emplace_back(
                  current_time, current_cartesian.x, current_cartesian.y, current_cartesian.z,
                  current_cartesian.vx, current_cartesian.vy, current_cartesian.vz, current_jacobi);
            }
          }
        }

        // 後方互換のためsali_final, sali_threshold_timeに設定
        double sali_final = chaos_final;
        double sali_threshold_time = chaos_threshold_time;

        // SALI時系列CSVファイル出力（設定で有効な場合のみ）
        if (config.output_sali_timeseries) {
          std::ostringstream ts_filename;
          ts_filename << timeseries_dir << "/" << file_basename << "_sali_ts_p" << phase_idx
                      << "_dv" << std::fixed << std::setprecision(4) << dv_mag << ".csv";
          std::string ts_csv_path = ts_filename.str();
          {
            std::ofstream ts_ofs(ts_csv_path);
            if (ts_ofs) {
              ts_ofs << std::setprecision(15) << std::fixed;
              ts_ofs << "# SALI Time Series\n";
              ts_ofs << "# phase_idx=" << phase_idx << ", dv_mag=" << dv_mag << "\n";
              ts_ofs << "# Source: " << orbit_file << "\n";
              ts_ofs << "time_nd,sali\n";
              for (const auto& ts_point : chaos_timeseries) {
                ts_ofs << ts_point.first << "," << ts_point.second << "\n";
              }
              ts_ofs.close();

              // SALI時系列プロット生成
              GenerateSaliTimeSeriesPlot(ts_csv_path, timeseries_dir, file_basename, phase_idx,
                                         dv_mag);
            }
          }
        }

        // 座標変換後軌道の状態量をCSVに出力（dv_mag=0の場合のみ、各位相点を1行）
        if (config.output_converted_trajectory && dv_idx == 0 && converted_ofs) {
#pragma omp critical(converted_csv)
          {
            double jacobi_rot = calc_jacobi_integral(asteroid_rot, kMU);
            converted_ofs << phase_idx << "," << row.time_j2000 << "," << asteroid_rot.x << ","
                          << asteroid_rot.y << "," << asteroid_rot.z << "," << asteroid_rot.vx
                          << "," << asteroid_rot.vy << "," << asteroid_rot.vz << "," << jacobi_rot
                          << "\n";
          }
        }

        // インパルス付与後軌道の時刻歴CSVファイル出力（dv_mag > 0の場合のみ）
        if (config.output_impulse_trajectory && dv_idx > 0) {
          std::ostringstream traj_filename;
          traj_filename << impulse_traj_dir << "/" << file_basename << "_impulse_traj_p"
                        << phase_idx << "_dv" << std::fixed << std::setprecision(4) << dv_mag
                        << ".csv";
          std::string traj_csv_path = traj_filename.str();
          std::ofstream traj_ofs(traj_csv_path);
          if (traj_ofs) {
            traj_ofs << std::setprecision(15) << std::fixed;
            traj_ofs << "# Impulse Trajectory (After Delta-V Application)\n";
            traj_ofs << "# phase_idx=" << phase_idx << ", dv_mag=" << dv_mag << "\n";
            traj_ofs << "# Source: " << orbit_file << "\n";
            traj_ofs << "# Initial state after impulse: x=" << modified_state.x
                     << ", y=" << modified_state.y << ", z=" << modified_state.z << "\n";
            traj_ofs << "# Delta-V: dvx=" << dv_x << ", dvy=" << dv_y << ", dvz=" << dv_z << "\n";
            traj_ofs << "time_nd,x,y,z,vx,vy,vz,jacobi\n";
            for (const auto& traj_point : trajectory_timeseries) {
              traj_ofs << std::get<0>(traj_point) << "," << std::get<1>(traj_point) << ","
                       << std::get<2>(traj_point) << "," << std::get<3>(traj_point) << ","
                       << std::get<4>(traj_point) << "," << std::get<5>(traj_point) << ","
                       << std::get<6>(traj_point) << "," << std::get<7>(traj_point) << "\n";
            }
            traj_ofs.close();
          }
        }

// CSV出力（排他制御）
#pragma omp critical
        {
          ofs << phase_idx << "," << row.time_j2000 << "," << dv_mag << "," << sali_final << ","
              << sali_threshold_time << "," << escape_flag << "," << collision_flag << ","
              << modified_state.x << "," << modified_state.y << "," << modified_state.z << ","
              << modified_state.vx << "," << modified_state.vy << "," << modified_state.vz << ","
              << jacobi << "\n";

          completed++;
          if (completed % 100 == 0 || completed == total_calcs) {
            double progress = static_cast<double>(completed) / total_calcs;
            displayProgressBarThreadSafe(progress);
          }
        }
      }
    }

    displayProgressBarThreadSafe(1.0);
    std::cout << std::endl;

    ofs.close();
    std::cout << "<>    Output: " << csv_path << std::endl;

    // gnuplotでカラーコンタ図生成
    GenerateSaliContourPlots(csv_path, output_dir, file_basename);

    auto end_file = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_file - start_file);
    std::cout << "<>    Elapsed: " << duration.count() / 1000 << "s " << duration.count() % 1000
              << "ms" << std::endl;
  }

  // 全体時間
  auto end_total = std::chrono::system_clock::now();
  auto duration_total =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
  auto sec = duration_total.count() / 1000;
  auto min = sec / 60;
  sec = sec % 60;

  std::cout << "<>" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    All calculations finished" << std::endl;
  std::cout << "<>    Total elapsed time: " << min << "m " << sec << "s" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  return 0;
}
