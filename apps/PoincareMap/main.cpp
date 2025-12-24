/**
 * @file main.cpp
 * @brief ポアンカレマップ生成アプリケーション
 * @details CR3BP（円制限三体問題）のポアンカレマップを生成する
 * @date 2025-12-23
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
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector>

#include "rtbp.hpp"

namespace fs = std::filesystem;
using namespace crtbp;
using namespace utils;
using namespace my_type;

// ---------------------------------------------------------------------------
// 断面交差方向の列挙型
// ---------------------------------------------------------------------------
enum class CrossingDirection {
  kPositive,  // 正方向のみ（変数が負→正）
  kNegative,  // 負方向のみ（変数が正→負）
  kBoth       // 両方向
};

// ---------------------------------------------------------------------------
// 積分器タイプの列挙型
// ---------------------------------------------------------------------------
enum class IntegratorType { kSymplectic4, kSymplectic6, kRK4, kRK45 };

// ---------------------------------------------------------------------------
// コンフィグ構造体
// ---------------------------------------------------------------------------
struct PoincareMapConfig {
  // シミュレーション設定
  double calc_time = 100.0;
  double timestep = 0.0001;

  // ポアンカレ断面設定
  std::string section_var = "y";
  double section_value = 0.0;
  CrossingDirection crossing_direction = CrossingDirection::kPositive;

  // 初期条件設定
  double jacobi_integral = 3.0008;
  std::string var1_name = "x";
  double var1_min = 0.98;
  double var1_max = 1.02;
  int var1_division = 20;
  std::string var2_name = "z";
  double var2_min = -0.01;
  double var2_max = 0.01;
  int var2_division = 20;
  std::string fixed_var1_name = "vx";
  double fixed_var1_value = 0.0;
  std::string fixed_var2_name = "vy";
  double fixed_var2_value = 0.0;

  // 積分器設定
  IntegratorType integrator_type = IntegratorType::kSymplectic6;
  // RK45（適応ステップ）用パラメータ
  double rk45_atol = 1e-10;          // 絶対誤差許容値
  double rk45_rtol = 1e-10;          // 相対誤差許容値
  double rk45_initial_step = 0.001;  // 初期ステップ幅
  double rk45_max_step = 0.01;       // 最大ステップ幅

  // 衝突判定設定
  double collision_radius_earth = 1e-12;  // 地球への最小距離（無次元）
  double collision_radius_sun = 1e-12;    // 太陽への最小距離（無次元）

  // 出力設定
  std::string output_tag = "";
};

// ---------------------------------------------------------------------------
// ポアンカレマップの交差点を表す構造体
// ---------------------------------------------------------------------------
struct PoincareCrossing {
  double time;
  double x, y, z, vx, vy, vz;
  int crossing_index;  // 何回目の交差か
};

// ---------------------------------------------------------------------------
// コマンドライン引数をパースする
// ---------------------------------------------------------------------------
void ParseCommandLineArgs(int argc, char* argv[], bool& is_continuous, bool& skip_wait,
                          std::string& output_tag) {
  is_continuous = false;
  skip_wait = false;
  output_tag = "";

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--continuous" || arg == "-c") {
      is_continuous = true;
    } else if (arg == "--no-wait" || arg == "-n") {
      skip_wait = true;
    } else if ((arg == "--tag" || arg == "-t") && i + 1 < argc) {
      output_tag = argv[++i];
    } else if (arg == "--help" || arg == "-h") {
      std::cout
          << "Usage: " << argv[0] << " [options]\n"
          << "Options:\n"
          << "  -c, --continuous  連続シミュレーションモード（複数configファイルを順次処理）\n"
          << "  -n, --no-wait     ユーザー確認のための待機をスキップ\n"
          << "  -t, --tag <TAG>   出力フォルダに付与するタグ\n"
          << "  -h, --help        このヘルプを表示\n"
          << std::endl;
      std::exit(0);
    }
  }
}

// ---------------------------------------------------------------------------
// 変数名からインデックスを取得
// ---------------------------------------------------------------------------
int GetVarIndex(const std::string& var_name) {
  if (var_name == "x") return 0;
  if (var_name == "y") return 1;
  if (var_name == "z") return 2;
  if (var_name == "vx") return 3;
  if (var_name == "vy") return 4;
  if (var_name == "vz") return 5;
  return -1;
}

// ---------------------------------------------------------------------------
// 状態ベクトルから指定インデックスの値を取得
// ---------------------------------------------------------------------------
double GetStateValue(const State<double>& state, int index) {
  switch (index) {
    case 0:
      return state.x;
    case 1:
      return state.y;
    case 2:
      return state.z;
    case 3:
      return state.vx;
    case 4:
      return state.vy;
    case 5:
      return state.vz;
    default:
      return 0.0;
  }
}

// ---------------------------------------------------------------------------
// 状態ベクトルの指定インデックスに値を設定
// ---------------------------------------------------------------------------
void SetStateValue(State<double>* state, int index, double value) {
  switch (index) {
    case 0:
      state->x = value;
      break;
    case 1:
      state->y = value;
      break;
    case 2:
      state->z = value;
      break;
    case 3:
      state->vx = value;
      break;
    case 4:
      state->vy = value;
      break;
    case 5:
      state->vz = value;
      break;
  }
}

// ---------------------------------------------------------------------------
// コンフィグファイルを読み込む
// ---------------------------------------------------------------------------
bool LoadConfig(const std::string& filepath, PoincareMapConfig* config) {
  try {
    TomlConfigParser parser(filepath);

    // シミュレーション設定
    config->calc_time = parser.GetDouble("simulation.calc_time", 100.0);
    config->timestep = parser.GetDouble("simulation.timestep", 0.0001);

    // ポアンカレ断面設定
    config->section_var = parser.GetString("section.section_var", "y");
    config->section_value = parser.GetDouble("section.section_value", 0.0);

    std::string crossing_str = parser.GetString("section.crossing_direction", "positive");
    if (crossing_str == "positive") {
      config->crossing_direction = CrossingDirection::kPositive;
    } else if (crossing_str == "negative") {
      config->crossing_direction = CrossingDirection::kNegative;
    } else if (crossing_str == "both") {
      config->crossing_direction = CrossingDirection::kBoth;
    }

    // 初期条件設定
    config->jacobi_integral = parser.GetDouble("initial_condition.jacobi_integral", 3.0008);
    config->var1_name = parser.GetString("initial_condition.var1_name", "x");
    config->var1_min = parser.GetDouble("initial_condition.var1_min", 0.98);
    config->var1_max = parser.GetDouble("initial_condition.var1_max", 1.02);
    config->var1_division = parser.GetInt("initial_condition.var1_division", 20);
    config->var2_name = parser.GetString("initial_condition.var2_name", "z");
    config->var2_min = parser.GetDouble("initial_condition.var2_min", -0.01);
    config->var2_max = parser.GetDouble("initial_condition.var2_max", 0.01);
    config->var2_division = parser.GetInt("initial_condition.var2_division", 20);
    config->fixed_var1_name = parser.GetString("initial_condition.fixed_var1_name", "vx");
    config->fixed_var1_value = parser.GetDouble("initial_condition.fixed_var1_value", 0.0);
    config->fixed_var2_name = parser.GetString("initial_condition.fixed_var2_name", "vy");
    config->fixed_var2_value = parser.GetDouble("initial_condition.fixed_var2_value", 0.0);

    // 積分器設定
    std::string integrator_str = parser.GetString("integrator.type", "symplectic6");
    if (integrator_str == "symplectic4") {
      config->integrator_type = IntegratorType::kSymplectic4;
    } else if (integrator_str == "symplectic6") {
      config->integrator_type = IntegratorType::kSymplectic6;
    } else if (integrator_str == "rk4") {
      config->integrator_type = IntegratorType::kRK4;
    } else if (integrator_str == "rk45" || integrator_str == "dopri5") {
      config->integrator_type = IntegratorType::kRK45;
    }

    // RK45適応ステップ用パラメータ
    config->rk45_atol = parser.GetDouble("integrator.rk45_atol", 1e-10);
    config->rk45_rtol = parser.GetDouble("integrator.rk45_rtol", 1e-10);
    config->rk45_initial_step = parser.GetDouble("integrator.rk45_initial_step", 0.001);
    config->rk45_max_step = parser.GetDouble("integrator.rk45_max_step", 0.01);

    // 衝突判定設定
    config->collision_radius_earth = parser.GetDouble("collision.collision_radius_earth", 1e-12);
    config->collision_radius_sun = parser.GetDouble("collision.collision_radius_sun", 1e-12);

    // 出力設定
    config->output_tag = parser.GetString("output.tag", "");

    return true;
  } catch (const std::exception& e) {
    std::cerr << "Error loading config: " << e.what() << std::endl;
    return false;
  }
}

// ---------------------------------------------------------------------------
// 設定内容を表示
// ---------------------------------------------------------------------------
void PrintConfig(const PoincareMapConfig& config) {
  std::cout << "<>    [Config Summary]" << std::endl;
  std::cout << "<>        calc_time: " << config.calc_time << std::endl;
  std::cout << "<>        timestep: " << config.timestep << std::endl;
  std::cout << "<>        section_var: " << config.section_var << " = " << config.section_value
            << std::endl;

  std::string crossing_str;
  switch (config.crossing_direction) {
    case CrossingDirection::kPositive:
      crossing_str = "positive";
      break;
    case CrossingDirection::kNegative:
      crossing_str = "negative";
      break;
    case CrossingDirection::kBoth:
      crossing_str = "both";
      break;
  }
  std::cout << "<>        crossing_direction: " << crossing_str << std::endl;
  std::cout << "<>        jacobi_integral: " << config.jacobi_integral << std::endl;
  std::cout << "<>        var1: " << config.var1_name << " [" << config.var1_min << ", "
            << config.var1_max << "] div=" << config.var1_division << std::endl;
  std::cout << "<>        var2: " << config.var2_name << " [" << config.var2_min << ", "
            << config.var2_max << "] div=" << config.var2_division << std::endl;
  std::cout << "<>        fixed: " << config.fixed_var1_name << "=" << config.fixed_var1_value
            << ", " << config.fixed_var2_name << "=" << config.fixed_var2_value << std::endl;

  std::string integrator_str;
  switch (config.integrator_type) {
    case IntegratorType::kSymplectic4:
      integrator_str = "symplectic4";
      break;
    case IntegratorType::kSymplectic6:
      integrator_str = "symplectic6";
      break;
    case IntegratorType::kRK4:
      integrator_str = "rk4";
      break;
    case IntegratorType::kRK45:
      integrator_str = "rk45 (Dormand-Prince)";
      break;
  }
  std::cout << "<>        integrator: " << integrator_str << std::endl;
  if (config.integrator_type == IntegratorType::kRK45) {
    std::cout << "<>          rk45_atol: " << config.rk45_atol << std::endl;
    std::cout << "<>          rk45_rtol: " << config.rk45_rtol << std::endl;
    std::cout << "<>          rk45_initial_step: " << config.rk45_initial_step << std::endl;
    std::cout << "<>          rk45_max_step: " << config.rk45_max_step << std::endl;
  }
  std::cout << "<>        output_tag: " << config.output_tag << std::endl;
}

// ---------------------------------------------------------------------------
// 断面交差を検出する
// ---------------------------------------------------------------------------
bool DetectCrossing(double prev_val, double curr_val, double section_value,
                    CrossingDirection direction) {
  double prev_diff = prev_val - section_value;
  double curr_diff = curr_val - section_value;

  // 断面を通過したかチェック（符号が変化）
  if (prev_diff * curr_diff >= 0) {
    return false;  // 符号変化なし
  }

  // 方向チェック
  switch (direction) {
    case CrossingDirection::kPositive:
      return prev_diff < 0 && curr_diff > 0;
    case CrossingDirection::kNegative:
      return prev_diff > 0 && curr_diff < 0;
    case CrossingDirection::kBoth:
      return true;
  }
  return false;
}

// ---------------------------------------------------------------------------
// 線形補間で交差点を求める
// ---------------------------------------------------------------------------
State<double> InterpolateCrossing(const State<double>& prev_state, const State<double>& curr_state,
                                  int section_index, double section_value, double* crossing_time,
                                  double prev_time, double curr_time) {
  double prev_val = GetStateValue(prev_state, section_index);
  double curr_val = GetStateValue(curr_state, section_index);

  double alpha = (section_value - prev_val) / (curr_val - prev_val);

  State<double> crossing_state;
  crossing_state.x = prev_state.x + alpha * (curr_state.x - prev_state.x);
  crossing_state.y = prev_state.y + alpha * (curr_state.y - prev_state.y);
  crossing_state.z = prev_state.z + alpha * (curr_state.z - prev_state.z);
  crossing_state.vx = prev_state.vx + alpha * (curr_state.vx - prev_state.vx);
  crossing_state.vy = prev_state.vy + alpha * (curr_state.vy - prev_state.vy);
  crossing_state.vz = prev_state.vz + alpha * (curr_state.vz - prev_state.vz);

  *crossing_time = prev_time + alpha * (curr_time - prev_time);

  return crossing_state;
}

// ---------------------------------------------------------------------------
// ヤコビ積分から残り1変数を導出
// ---------------------------------------------------------------------------
bool DeriveRemainingVar(State<double>* state, int derive_index, double jacobi_target, double mu) {
  if (derive_index >= 3) {
    double x = state->x, y = state->y, z = state->z;
    // double r1 = std::sqrt((x + mu) * (x + mu) + y * y + z * z);
    // double r2 = std::sqrt((x - 1.0 + mu) * (x - 1.0 + mu) + y * y + z * z);

    // double U = x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 + mu * (1.0 - mu);

    double U = calc_potential_U(x, y, z, mu);
    double v_squared_target = 2.0 * U - jacobi_target;

    if (v_squared_target < 0) {
      return false;
    }

    double other_v_squared = 0.0;
    if (derive_index != 3) other_v_squared += state->vx * state->vx;
    if (derive_index != 4) other_v_squared += state->vy * state->vy;
    if (derive_index != 5) other_v_squared += state->vz * state->vz;

    double derive_v_squared = v_squared_target - other_v_squared;
    if (derive_v_squared < 0) {
      return false;
    }

    SetStateValue(state, derive_index, std::sqrt(derive_v_squared));
    return true;
  }

  return false;
}

// ---------------------------------------------------------------------------
// gnuplotスクリプトを生成
// ---------------------------------------------------------------------------
void GenerateGnuplotScript(const std::string& output_dir, const std::string& data_filename,
                           const std::string& var1_name, const std::string& var2_name,
                           const PoincareMapConfig& config) {
  std::string gp_path = output_dir + "/poincare_map.gp";
  std::ofstream gp(gp_path);

  double xrange_min = config.var1_min - std::sqrt(config.var1_min * config.var1_min) * 0.1;
  double xrange_max = config.var1_max + std::sqrt(config.var1_max * config.var1_max) * 0.1;
  double yrange_min = config.var2_min - std::sqrt(config.var2_min * config.var2_min) * 0.1;
  double yrange_max = config.var2_max + std::sqrt(config.var2_max * config.var2_max) * 0.1;

  gp << "# Poincare Map gnuplot script\n";
  gp << "set datafile separator \",\"\n";
  gp << "set terminal pngcairo size 800,800 enhanced font 'Arial,12'\n";
  if (xrange_min != xrange_max) gp << "set xrange[" << xrange_min << ":" << xrange_max << "]\n";
  if (yrange_min != yrange_max) gp << "set yrange[" << yrange_min << ":" << yrange_max << "]\n";
  gp << "set output '" << output_dir << "/poincare_map.png'\n";
  gp << "\n";
  gp << "set title 'Poincare Map (" << config.section_var << " = " << config.section_value
     << ")' font ',14'\n";
  gp << "set xlabel '" << var1_name << "' font ',12'\n";
  gp << "set ylabel '" << var2_name << "' font ',12'\n";
  gp << "set size ratio 1\n";
  gp << "set grid\n";
  gp << "unset key\n";
  gp << "\n";
  gp << "plot '" << data_filename << "' using 1:2 with points pt 7 ps 0.3 lc rgb '#0066CC'\n";
  gp << "\n";
  gp << "# EPS output\n";
  gp << "set terminal postscript eps enhanced color font 'Arial,12'\n";
  gp << "set output '" << output_dir << "/poincare_map.eps'\n";
  gp << "replot\n";

  gp.close();
}

// ---------------------------------------------------------------------------
// メイン関数
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP Poincare Map Generator ver1.0" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // コマンドライン引数のパース
  bool is_continuous = false;
  bool skip_wait = false;
  std::string cli_output_tag = "";
  ParseCommandLineArgs(argc, argv, is_continuous, skip_wait, cli_output_tag);

  // 天文パラメータの読み込み
  const std::string kConfigFilePath = CONFIG_DIR;
  std::string astro_param_file = kConfigFilePath + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  std::cout << "<>    mu = " << std::setprecision(10) << kMU << std::endl;

  // 選択されたモードを表示
  if (is_continuous) {
    std::cout << "<>    selected mode : continuous simulation (--continuous)\n" << std::endl;
  } else {
    std::cout << "<> >  selected mode : single simulation (default)" << std::endl;
  }

  if (skip_wait) {
    std::cout << "<>    WaitForEnter : skipped (--no-wait)" << std::endl;
  }

  // コンフィグファイルのリストを取得
  const std::string kCalcConfigPath = kConfigFilePath + "/poincare_map/";
  std::string config_prefix = "poincare_map_config";
  std::vector<std::string> config_file_list;

  const std::regex pattern(is_continuous ? "^" + config_prefix + "_\\d+\\.toml$"
                                         : "^" + config_prefix + "\\.toml$");

  try {
    for (const auto& entry : fs::directory_iterator(kCalcConfigPath)) {
      if (entry.is_regular_file()) {
        std::string filename = entry.path().filename().string();
        if (std::regex_match(filename, pattern)) {
          config_file_list.push_back(fs::absolute(entry.path()).string());
        }
      }
    }
  } catch (const fs::filesystem_error& e) {
    std::cerr << "Error accessing config directory: " << e.what() << std::endl;
    return -1;
  }

  if (config_file_list.empty()) {
    std::cerr << "<> No config files found in " << kCalcConfigPath << std::endl;
    return -1;
  }

  std::sort(config_file_list.begin(), config_file_list.end());

  std::cout << "<> Loaded config file list:" << std::endl;
  for (const auto& filepath : config_file_list) {
    std::cout << "<>    - " << filepath << std::endl;
  }

  if (!skip_wait) {
    WaitForEnter();
  }

  std::cout << "<> " << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  // OpenMP スレッド数の設定
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
    if (getcore <= Core_Max && getcore > 0) {
      OMP_Fmax = getcore;
      std::cout << "<>" << std::endl;
      std::cout << "<>     >> Number of OMP threads is " << OMP_Fmax << std::endl;
      std::cout << "<>" << std::endl;
    } else {
      OMP_Fmax = Core_Max - 1;
      std::cout << "  <>     >> Your input is INVALID. OMP threads is "
                << "automatically determined as " << OMP_Fmax << std::endl;
    }
    WaitForEnter();
  } else {
    std::cout << "  <>     >> OMP threads is automatically determined as " << OMP_Fmax << std::endl;
  }
  omp_set_num_threads(OMP_Fmax);
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  // RK4用の運動方程式ファンクタ
  EquationOfMotion<double> eom(astro_params);

  // セッション出力ディレクトリを作成
  std::string output_base_path = OUTPUT_DIR;
  std::string session_dir_name = getcurrent_date();
  if (!cli_output_tag.empty()) {
    session_dir_name += "_" + cli_output_tag;
  }
  std::string session_output_dir = output_base_path + "/poincare_map/" + session_dir_name;
  if (!fs::exists(session_output_dir)) {
    fs::create_directories(session_output_dir);
  }
  std::cout << "<>" << std::endl;
  std::cout << "<>    Session output directory: " << session_output_dir << std::endl;
  std::cout << "<>" << std::endl;

  auto start_ofall = std::chrono::system_clock::now();

  // 各コンフィグファイルを処理
  for (const auto& configfilepath : config_file_list) {
    std::cout << "\n<>----------------------------------------------------------------"
              << std::endl;
    std::cout << "<>    Processing: " << configfilepath << std::endl;

    PoincareMapConfig config;
    if (!LoadConfig(configfilepath, &config)) {
      std::cerr << "<>    Failed to load config, skipping..." << std::endl;
      continue;
    }

    PrintConfig(config);

    if (!is_continuous && !skip_wait) {
      WaitForEnter();
    }

    // 出力ディレクトリ作成
    std::string config_basename = fs::path(configfilepath).stem().string();
    std::string output_dir = session_output_dir + "/" + config_basename;
    if (!config.output_tag.empty()) {
      output_dir += "_" + config.output_tag;
    }
    fs::create_directories(output_dir);

    // 使用したconfigファイルを出力ディレクトリにコピー
    std::string config_copy_path = output_dir + "/" + fs::path(configfilepath).filename().string();
    try {
      fs::copy_file(configfilepath, config_copy_path, fs::copy_options::overwrite_existing);
      std::cout << "<>    Config file copied to: " << config_copy_path << std::endl;
    } catch (const fs::filesystem_error& e) {
      std::cerr << "<>    Warning: Failed to copy config file: " << e.what() << std::endl;
    }

    std::cout << "<>    Output directory: " << output_dir << std::endl;

    auto start = std::chrono::system_clock::now();

    // 変数インデックスの取得
    int section_index = GetVarIndex(config.section_var);
    int var1_index = GetVarIndex(config.var1_name);
    int var2_index = GetVarIndex(config.var2_name);
    int fixed1_index = GetVarIndex(config.fixed_var1_name);
    int fixed2_index = GetVarIndex(config.fixed_var2_name);

    // 導出する変数のインデックスを特定
    int derive_index = -1;
    for (int i = 0; i < 6; ++i) {
      if (i != section_index && i != var1_index && i != var2_index && i != fixed1_index &&
          i != fixed2_index) {
        derive_index = i;
        break;
      }
    }

    if (derive_index < 0) {
      std::cerr << "<>    Error: Cannot determine variable to derive from Jacobi integral"
                << std::endl;
      continue;
    }

    std::cout << "<>    Deriving variable index: " << derive_index << std::endl;

    // 初期条件グリッドの生成
    int total_points = config.var1_division * config.var2_division;
    double var1_step = (config.var1_max - config.var1_min) / std::max(1, config.var1_division - 1);
    double var2_step = (config.var2_max - config.var2_min) / std::max(1, config.var2_division - 1);

    std::cout << "<>    Generating Poincare map (" << total_points << " initial conditions)..."
              << std::endl;

    // OpenMP並列化のための共有変数
    std::vector<std::vector<PoincareCrossing>> thread_crossings(OMP_Fmax);
    int completed_count = 0;
    int display_interval = std::max(total_points / 100, 1);

// OpenMP並列化ループ
#ifndef __debug__
#pragma omp parallel shared(thread_crossings, completed_count, total_points)
#endif
    {
      int thread_id = omp_get_thread_num();

#ifndef __debug__
#pragma omp for schedule(dynamic)
#endif
      for (int idx = 0; idx < total_points; ++idx) {
        int i = idx / config.var2_division;
        int j = idx % config.var2_division;

        double var1_val = config.var1_min + i * var1_step;
        double var2_val = config.var2_min + j * var2_step;

        // 初期状態を構築
        State<double> state = {0, 0, 0, 0, 0, 0};
        SetStateValue(&state, section_index, config.section_value);
        SetStateValue(&state, var1_index, var1_val);
        SetStateValue(&state, var2_index, var2_val);
        SetStateValue(&state, fixed1_index, config.fixed_var1_value);
        SetStateValue(&state, fixed2_index, config.fixed_var2_value);

        // ヤコビ積分から残り1変数を導出
        if (!DeriveRemainingVar(&state, derive_index, config.jacobi_integral, kMU)) {
#ifndef __debug__
#pragma omp atomic
#endif
          completed_count++;
          continue;
        }
#ifdef __debug___
        std::cout << "state " << state.x << " " << state.y << " " << state.z << " " << state.vx
                  << " " << state.vy << " " << state.vz << std::endl;
#endif
        // 軌道積分とポアンカレ断面交差検出
        State<double> prev_state = state;
        double time = 0.0;
        int crossing_count = 0;

        // RK45の場合は適応的ステップ、それ以外は固定ステップ
        if (config.integrator_type == IntegratorType::kRK45) {
          // Dormand-Prince RK45 (適応的ステップサイズ制御)
          namespace odeint = boost::numeric::odeint;
          using Stepper = odeint::runge_kutta_dopri5<State<double>, double, State<double>, double,
                                                     odeint::vector_space_algebra>;
          using ControlledStepper = odeint::controlled_runge_kutta<Stepper>;

          Dopri5OrbitSystem<double> dopri5_system(astro_params);
          ControlledStepper controlled_stepper =
              odeint::make_controlled(config.rk45_atol, config.rk45_rtol, Stepper());

          double dt = config.rk45_initial_step;
          State<double> curr_state = prev_state;

          while (time < config.calc_time) {
            // ステップサイズを最大値で制限
            double step_to_take = std::min(dt, config.rk45_max_step);
            step_to_take = std::min(step_to_take, config.calc_time - time);

            prev_state = curr_state;
            double prev_time = time;

            // 適応的ステップで1ステップ進める
            auto result =
                controlled_stepper.try_step(dopri5_system, curr_state, time, step_to_take);
            if (result == odeint::controlled_step_result::success) {
              dt = step_to_take;  // 次回も同じステップ幅を試す
            } else if (result == odeint::controlled_step_result::fail) {
              dt = step_to_take / 2.0;  // ステップ幅を縮小してリトライ
              curr_state = prev_state;
              time = prev_time;
              continue;
            }

            // 衝突判定（地球・太陽への距離チェック）
            double r_earth = calc_r2(curr_state.x, curr_state.y, curr_state.z, kMU);
            double r_sun = calc_r1(curr_state.x, curr_state.y, curr_state.z, kMU);
            if (r_earth < config.collision_radius_earth || r_sun < config.collision_radius_sun) {
              // 天体に接近しすぎたため計算を打ち切り
              break;
            }

            // ポアンカレ断面交差検出
            double prev_section_val = GetStateValue(prev_state, section_index);
            double curr_section_val = GetStateValue(curr_state, section_index);

            if (DetectCrossing(prev_section_val, curr_section_val, config.section_value,
                               config.crossing_direction)) {
              double crossing_time;
              State<double> crossing_state =
                  InterpolateCrossing(prev_state, curr_state, section_index, config.section_value,
                                      &crossing_time, prev_time, time);

              PoincareCrossing crossing;
              crossing.time = crossing_time;
              crossing.x = crossing_state.x;
              crossing.y = crossing_state.y;
              crossing.z = crossing_state.z;
              crossing.vx = crossing_state.vx;
              crossing.vy = crossing_state.vy;
              crossing.vz = crossing_state.vz;
              crossing.crossing_index = ++crossing_count;

              thread_crossings[thread_id].push_back(crossing);
            }
          }
        } else {
          // 固定ステップ積分器 (Symplectic, RK4)
          int num_steps = static_cast<int>(config.calc_time / config.timestep);

          for (int step = 0; step < num_steps; ++step) {
            State<double> curr_state;

            switch (config.integrator_type) {
              case IntegratorType::kSymplectic4:
                curr_state = SymplecticStep4thOrder(kMU, prev_state, config.timestep);
                break;
              case IntegratorType::kSymplectic6:
                curr_state = SymplecticStep6thOrder(kMU, prev_state, config.timestep);
                break;
              case IntegratorType::kRK4:
                curr_state = RungeKutta4Step(eom, prev_state, time, config.timestep);
                break;
              default:
                curr_state = prev_state;
                break;
            }

            double curr_time = time + config.timestep;

            // 衝突判定（地球・太陽への距離チェック）
            double r_earth = calc_r2(curr_state.x, curr_state.y, curr_state.z, kMU);
            double r_sun = calc_r1(curr_state.x, curr_state.y, curr_state.z, kMU);
            if (r_earth < config.collision_radius_earth || r_sun < config.collision_radius_sun) {
              // 天体に接近しすぎたため計算を打ち切り
              break;
            }

            double prev_section_val = GetStateValue(prev_state, section_index);
            double curr_section_val = GetStateValue(curr_state, section_index);

            if (DetectCrossing(prev_section_val, curr_section_val, config.section_value,
                               config.crossing_direction)) {
              double crossing_time;
              State<double> crossing_state =
                  InterpolateCrossing(prev_state, curr_state, section_index, config.section_value,
                                      &crossing_time, time, curr_time);

              PoincareCrossing crossing;
              crossing.time = crossing_time;
              crossing.x = crossing_state.x;
              crossing.y = crossing_state.y;
              crossing.z = crossing_state.z;
              crossing.vx = crossing_state.vx;
              crossing.vy = crossing_state.vy;
              crossing.vz = crossing_state.vz;
              crossing.crossing_index = ++crossing_count;

              thread_crossings[thread_id].push_back(crossing);
            }

            prev_state = curr_state;
            time = curr_time;
          }
        }
#ifndef __debug__
#pragma omp atomic
#endif
        completed_count++;

        if (omp_get_thread_num() == 0) {
          if (completed_count % display_interval == 0 ||
              completed_count >= total_points - OMP_Fmax) {
            double progress = static_cast<double>(completed_count) / total_points;
            displayProgressBarThreadSafe(progress);
          }
        }
      }
    }

    displayProgressBarThreadSafe(1.0);
    std::cout << std::endl;

    // 全スレッドの結果を統合
    std::vector<PoincareCrossing> all_crossings;
    for (const auto& tc : thread_crossings) {
      all_crossings.insert(all_crossings.end(), tc.begin(), tc.end());
    }

    std::cout << "<>    Total crossings detected: " << all_crossings.size() << std::endl;

    // ポアンカレマップデータを出力
    std::string poincare_data_path = output_dir + "/poincare_map.csv";
    std::ofstream poincare_file(poincare_data_path);
    poincare_file << std::fixed << std::setprecision(15);
    poincare_file << "# Poincare Map Data\n";
    poincare_file << "# " << config.var1_name << "," << config.var2_name
                  << ",time,crossing_index\n";

    for (const auto& crossing : all_crossings) {
      double plot_x = GetStateValue(
          State<double>{crossing.x, crossing.y, crossing.z, crossing.vx, crossing.vy, crossing.vz},
          var1_index);
      double plot_y = GetStateValue(
          State<double>{crossing.x, crossing.y, crossing.z, crossing.vx, crossing.vy, crossing.vz},
          var2_index);
      poincare_file << plot_x << "," << plot_y << "," << crossing.time << ","
                    << crossing.crossing_index << "\n";
    }
    poincare_file.close();

    GenerateGnuplotScript(output_dir, poincare_data_path, config.var1_name, config.var2_name,
                          config);

    std::string gp_command = "gnuplot \"" + output_dir + "/poincare_map.gp\"";
    int gp_result = system(gp_command.c_str());
    if (gp_result == 0) {
      std::cout << "<>    Generated: poincare_map.png, poincare_map.eps" << std::endl;
    } else {
      std::cout << "<>    Warning: gnuplot execution failed (code " << gp_result << ")"
                << std::endl;
    }

    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    auto msec = duration.count() % 1000;
    auto sec = duration.count() / 1000 % 60;
    auto min = duration.count() / 1000 / 60 % 60;
    auto hour = duration.count() / 1000 / 60 / 60;
    std::cout << "<>        elapsed time : " << hour << "h " << min << "m " << sec << "s " << msec
              << "ms" << std::endl;
    std::cout << "<>        Simulation for " << configfilepath << " finished" << std::endl;
  }

  auto end_ofall = std::chrono::system_clock::now();
  auto total_duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_ofall - start_ofall);
  auto total_msec = total_duration.count() % 1000;
  auto total_sec = total_duration.count() / 1000 % 60;
  auto total_min = total_duration.count() / 1000 / 60 % 60;
  auto total_hour = total_duration.count() / 1000 / 60 / 60;

  std::cout << "\n<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    All simulations completed" << std::endl;
  std::cout << "<>    Total elapsed time: " << total_hour << "h " << total_min << "m " << total_sec
            << "s " << total_msec << "ms" << std::endl;

  return 0;
}
