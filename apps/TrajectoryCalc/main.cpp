/**
 * @file main.cpp
 * @brief 円制限三体問題の回転座標系における軌道計算アプリケーション
 *
 * 設定ファイルから初期条件を読み込み、6次シンプレクティック積分で軌道を計算し、
 * CSVとgnuplotプロットを出力する。
 *
 * 複数のconfigファイルに対応し、各ファイル内の複数COORDを処理する。
 * 出力は同じconfigファイル由来の軌道を1ファイルにまとめ、
 * gnuplotのインデックス機能用に軌道間に空行を挿入する。
 */

#include <algorithm>
#include <boost/numeric/odeint.hpp>
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

// IntegratorType と ChaosIndexType は utils.hpp で定義されているため、ローカル定義は不要
using utils::ChaosIndexType;
using utils::IntegratorType;

/**
 * @brief インパルス（速度変更）設定
 */
struct ImpulseConfig {
  bool enabled = false;                               ///< インパルス有効フラグ
  std::string trigger_type = "position";              ///< "position", "time", "distance_from_point"
  std::string trigger_axis = "x";                     ///< x, y, z, r, r_earth (position型のみ)
  double trigger_value = 0.0;                         ///< トリガー閾値
  std::string trigger_direction = "any";              ///< "any", "increasing", "decreasing"
  std::string velocity_mode = "delta";                ///< "delta", "absolute"
  std::array<double, 3> delta_v = {0, 0, 0};          ///< 速度増分 (delta モード)
  std::array<double, 3> target_velocity = {0, 0, 0};  ///< 目標速度 (absolute モード)
  std::array<double, 3> trigger_point = {0, 0, 0};    ///< トリガー位置 (distance_from_point用)
  double trigger_radius = 0.01;                       ///< トリガー距離閾値 (distance_from_point用)
  bool applied = false;                               ///< 適用済みフラグ (one_shot用)
};

/**
 * @brief 設定ファイルから初期条件を読み込む
 */
struct TrajectoryConfig {
  double calc_timestep = 0.0001;
  double time_threshold = 10.0;
  std::vector<my_type::State<double>> initial_coords;
  ChaosIndexType chaos_index_type = ChaosIndexType::NONE;  ///< カオス指標の種類
  int gali_k = 2;                                          ///< GALIの偏差ベクトル数 (2, 4, 6)
  IntegratorType integrator_type = IntegratorType::kSymplectic6th;  ///< 数値積分法
  bool enable_freq_analysis = false;                                ///< 周波数解析を有効にするか

  // LLE (最大リヤプノフ指数) 用パラメータ
  int lle_renorm_interval = 100;  ///< 再正規化の間隔（ステップ数）
  int lle_skip_steps = 1000;      ///< 初期過渡状態をスキップ（ステップ数）

  // DOPRI (Dormand-Prince) 用パラメータ
  double dopri_abs_tol = 1.0e-10;     ///< 絶対誤差許容値
  double dopri_rel_tol = 1.0e-10;     ///< 相対誤差許容値
  double dopri_initial_step = 0.001;  ///< 初期ステップサイズ
  double dopri_max_step = 0.1;        ///< 最大ステップサイズ

  // 出力制御フラグ
  bool output_trajectory = true;        ///< 軌道データのCSV出力
  bool output_chaos_index = true;       ///< カオス指標の時系列出力
  bool output_orbital_elements = true;  ///< 軌道要素の時系列出力
  bool output_freq_analysis = true;     ///< 周波数解析結果の出力
  bool output_gnuplot = true;           ///< gnuplotスクリプトと画像の生成

  // インパルス設定
  ImpulseConfig impulse;  ///< インパルス（速度変更）設定
};

// TrimString -> utils::trim() に置換

/**
 * @brief トリガー軸の値を取得
 * @param state 状態ベクトル
 * @param axis 軸名 ("x", "y", "z", "r", "r_earth")
 * @param mu 質量パラメータ
 * @return トリガー軸の値
 */
double GetTriggerAxisValue(const my_type::State<double>& state, const std::string& axis,
                           double mu) {
  if (axis == "x") return state.x;
  if (axis == "y") return state.y;
  if (axis == "z") return state.z;
  if (axis == "r") return std::sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  if (axis == "r_earth") {
    double dx = state.x - (1.0 - mu);
    return std::sqrt(dx * dx + state.y * state.y + state.z * state.z);
  }
  return 0.0;
}

/**
 * @brief インパルストリガー条件をチェック
 * @param prev_state 前ステップの状態
 * @param current_state 現在の状態
 * @param prev_time 前ステップの時刻
 * @param current_time 現在の時刻
 * @param impulse インパルス設定
 * @param mu 質量パラメータ
 * @return トリガー条件を満たした場合true
 */
bool CheckImpulseTrigger(const my_type::State<double>& prev_state,
                         const my_type::State<double>& current_state, double prev_time,
                         double current_time, const ImpulseConfig& impulse, double mu) {
  if (impulse.trigger_type == "time") {
    // 時刻ベースのトリガー
    return prev_time < impulse.trigger_value && current_time >= impulse.trigger_value;
  }
  if (impulse.trigger_type == "distance_from_point") {
    // 指定点からの距離ベースのトリガー
    double prev_dx = prev_state.x - impulse.trigger_point[0];
    double prev_dy = prev_state.y - impulse.trigger_point[1];
    double prev_dz = prev_state.z - impulse.trigger_point[2];
    double prev_dist = std::sqrt(prev_dx * prev_dx + prev_dy * prev_dy + prev_dz * prev_dz);

    double curr_dx = current_state.x - impulse.trigger_point[0];
    double curr_dy = current_state.y - impulse.trigger_point[1];
    double curr_dz = current_state.z - impulse.trigger_point[2];
    double curr_dist = std::sqrt(curr_dx * curr_dx + curr_dy * curr_dy + curr_dz * curr_dz);

    // 距離が閾値を下回った（外から内へ進入）
    return prev_dist > impulse.trigger_radius && curr_dist <= impulse.trigger_radius;
  }
  // 位置ベースのトリガー (position)
  double prev_val = GetTriggerAxisValue(prev_state, impulse.trigger_axis, mu);
  double curr_val = GetTriggerAxisValue(current_state, impulse.trigger_axis, mu);
  double threshold = impulse.trigger_value;

  bool crossed = (prev_val < threshold && curr_val >= threshold) ||
                 (prev_val > threshold && curr_val <= threshold);
  if (!crossed) return false;

  // 方向チェック
  if (impulse.trigger_direction == "increasing") {
    return prev_val < threshold && curr_val >= threshold;
  } else if (impulse.trigger_direction == "decreasing") {
    return prev_val > threshold && curr_val <= threshold;
  }
  // "any"
  return true;
}

/**
 * @brief 線形補間でトリガー位置の割合を計算
 * @param prev_state 前ステップの状態
 * @param current_state 現在の状態
 * @param prev_time 前ステップの時刻
 * @param current_time 現在の時刻
 * @param impulse インパルス設定
 * @param mu 質量パラメータ
 * @return [0, 1] の割合（0=前ステップ、1=現在ステップ）
 */
double ComputeTriggerFraction(const my_type::State<double>& prev_state,
                              const my_type::State<double>& current_state, double prev_time,
                              double current_time, const ImpulseConfig& impulse, double mu) {
  if (impulse.trigger_type == "time") {
    double dt = current_time - prev_time;
    if (dt <= 0) return 0.0;
    return (impulse.trigger_value - prev_time) / dt;
  }
  if (impulse.trigger_type == "distance_from_point") {
    // 距離ベースの線形補間
    double prev_dx = prev_state.x - impulse.trigger_point[0];
    double prev_dy = prev_state.y - impulse.trigger_point[1];
    double prev_dz = prev_state.z - impulse.trigger_point[2];
    double prev_dist = std::sqrt(prev_dx * prev_dx + prev_dy * prev_dy + prev_dz * prev_dz);

    double curr_dx = current_state.x - impulse.trigger_point[0];
    double curr_dy = current_state.y - impulse.trigger_point[1];
    double curr_dz = current_state.z - impulse.trigger_point[2];
    double curr_dist = std::sqrt(curr_dx * curr_dx + curr_dy * curr_dy + curr_dz * curr_dz);

    double ddist = curr_dist - prev_dist;
    if (std::abs(ddist) < 1e-15) return 0.0;
    return (impulse.trigger_radius - prev_dist) / ddist;
  }
  // 位置ベース (position)
  double prev_val = GetTriggerAxisValue(prev_state, impulse.trigger_axis, mu);
  double curr_val = GetTriggerAxisValue(current_state, impulse.trigger_axis, mu);
  double dval = curr_val - prev_val;
  if (std::abs(dval) < 1e-15) return 0.0;
  return (impulse.trigger_value - prev_val) / dval;
}

/**
 * @brief 状態を線形補間
 */
my_type::State<double> InterpolateState(const my_type::State<double>& s1,
                                        const my_type::State<double>& s2, double alpha) {
  return {s1.x + alpha * (s2.x - s1.x),    s1.y + alpha * (s2.y - s1.y),
          s1.z + alpha * (s2.z - s1.z),    s1.vx + alpha * (s2.vx - s1.vx),
          s1.vy + alpha * (s2.vy - s1.vy), s1.vz + alpha * (s2.vz - s1.vz)};
}

/**
 * @brief インパルスを適用
 * @param state 状態ベクトル（速度を変更）
 * @param impulse インパルス設定
 */
void ApplyImpulse(my_type::State<double>* state, const ImpulseConfig& impulse) {
  if (impulse.velocity_mode == "delta") {
    state->vx += impulse.delta_v[0];
    state->vy += impulse.delta_v[1];
    state->vz += impulse.delta_v[2];
  } else if (impulse.velocity_mode == "absolute") {
    state->vx = impulse.target_velocity[0];
    state->vy = impulse.target_velocity[1];
    state->vz = impulse.target_velocity[2];
  }
}

/**
 * @brief 設定ファイルを解析してTrajectoryConfigを返す
 */
bool LoadTrajectoryConfig(const std::string& filepath, TrajectoryConfig* config) {
  try {
    utils::TomlConfigParser parser(filepath);

    // シミュレーション設定
    config->calc_timestep = parser.GetDouble("simulation.calc_timestep", 0.0001);
    config->time_threshold = parser.GetDouble("simulation.time_threshold", 10.0);

    // カオス指標設定
    std::string chaos_type = parser.GetString("chaos.index_type", "NONE");
    if (chaos_type == "NONE" || chaos_type == "none") {
      config->chaos_index_type = ChaosIndexType::NONE;
    } else if (chaos_type == "SALI" || chaos_type == "sali") {
      config->chaos_index_type = ChaosIndexType::SALI;
      config->gali_k = 2;
    } else if (chaos_type == "GALI2" || chaos_type == "gali2") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 2;
    } else if (chaos_type == "GALI4" || chaos_type == "gali4") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 4;
    } else if (chaos_type == "GALI6" || chaos_type == "gali6") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 6;
    } else if (chaos_type == "LLE" || chaos_type == "lle") {
      config->chaos_index_type = ChaosIndexType::LLE;
    }

    // LLE用パラメータ
    config->lle_renorm_interval = parser.GetInt("chaos.lle_renorm_interval", 100);
    config->lle_skip_steps = parser.GetInt("chaos.lle_skip_steps", 1000);

    // 積分器設定
    std::string integrator_type = parser.GetString("integrator.type", "SYMPLECTIC6");
    if (integrator_type == "SYMPLECTIC4" || integrator_type == "symplectic4") {
      config->integrator_type = IntegratorType::kSymplectic4th;
    } else if (integrator_type == "SYMPLECTIC6" || integrator_type == "symplectic6") {
      config->integrator_type = IntegratorType::kSymplectic6th;
    } else if (integrator_type == "DOPRI" || integrator_type == "dopri") {
      config->integrator_type = IntegratorType::kDormandPrince;
    } else if (integrator_type == "RK4" || integrator_type == "rk4") {
      config->integrator_type = IntegratorType::kRungeKutta4th;
    }

    // DOPRI用パラメータ
    config->dopri_abs_tol = parser.GetDouble("integrator.dopri_abs_tol", 1.0e-10);
    config->dopri_rel_tol = parser.GetDouble("integrator.dopri_rel_tol", 1.0e-10);
    config->dopri_initial_step = parser.GetDouble("integrator.dopri_initial_step", 0.001);
    config->dopri_max_step = parser.GetDouble("integrator.dopri_max_step", 0.1);

    // 周波数解析設定
    config->enable_freq_analysis = parser.GetBool("analysis.enable_freq_analysis", false);

    // 出力設定
    config->output_trajectory = parser.GetBool("output.output_trajectory", true);
    config->output_chaos_index = parser.GetBool("output.output_chaos_index", true);
    config->output_orbital_elements = parser.GetBool("output.output_orbital_elements", true);
    config->output_freq_analysis = parser.GetBool("output.output_freq_analysis", true);
    config->output_gnuplot = parser.GetBool("output.output_gnuplot", true);

    // 初期座標の読み込み
    config->initial_coords = parser.GetCoordsArray("coords");

    // インパルス設定の読み込み (オプション)
    // enabledフラグがあればそれを優先、なければtrigger_typeの有無で判断
    if (parser.HasKey("impulse.enabled")) {
      config->impulse.enabled = parser.GetBool("impulse.enabled", false);
    } else {
      config->impulse.enabled = parser.HasKey("impulse.trigger_type");
    }

    if (config->impulse.enabled) {
      config->impulse.trigger_type = parser.GetString("impulse.trigger_type", "position");
      config->impulse.trigger_axis = parser.GetString("impulse.trigger_axis", "x");
      config->impulse.trigger_value = parser.GetDouble("impulse.trigger_value", 0.0);
      config->impulse.trigger_direction = parser.GetString("impulse.trigger_direction", "any");
      config->impulse.velocity_mode = parser.GetString("impulse.velocity_mode", "delta");
      auto dv = parser.GetDoubleArray("impulse.delta_v");
      if (dv.size() >= 3) {
        config->impulse.delta_v = {dv[0], dv[1], dv[2]};
      }
      auto tv = parser.GetDoubleArray("impulse.target_velocity");
      if (tv.size() >= 3) {
        config->impulse.target_velocity = {tv[0], tv[1], tv[2]};
      }
      // distance_from_point 用パラメータ
      auto tp = parser.GetDoubleArray("impulse.trigger_point");
      if (tp.size() >= 3) {
        config->impulse.trigger_point = {tp[0], tp[1], tp[2]};
      }
      config->impulse.trigger_radius = parser.GetDouble("impulse.trigger_radius", 0.01);
    }

    return true;
  } catch (const std::exception& e) {
    std::cerr << "<> !err! Error loading config: " << e.what() << std::endl;
    return false;
  }
}

/**
 * @brief gnuplotスクリプトを生成してEPS/PNGを出力する
 *
 * 複数軌道をインデックスで区別してプロットする
 */
void GenerateGnuplot(const std::string& csv_path, const std::string& output_dir,
                     const std::string& base_name, int num_trajectories) {
  std::string gnuplot_script = output_dir + "/" + base_name + ".gp";
  std::string eps_path = output_dir + "/" + base_name + ".eps";
  std::string png_path = output_dir + "/" + base_name + ".png";

  std::ofstream gp(gnuplot_script);
  if (!gp) {
    std::cerr << "<> !err! Cannot create gnuplot script: " << gnuplot_script << std::endl;
    return;
  }

  gp << "set terminal postscript eps enhanced color font 'Helvetica,14'\n";
  gp << "set output '" << eps_path << "'\n";
  gp << "set datafile separator ','\n";
  gp << "set xlabel 'x (AU)'\n";
  gp << "set ylabel 'y (AU)'\n";
  gp << "set title 'Trajectory in Rotating Frame (X-Y)'\n";
  gp << "set size ratio -1\n";
  gp << "set grid\n";

  // 複数軌道をインデックスで区別してプロット
  gp << "plot ";
  for (int i = 0; i < num_trajectories; ++i) {
    if (i > 0) {
      gp << ", \\\n     ";
    }
    gp << "'" << csv_path << "' index " << i << " using 2:3 with lines title 'Trajectory "
       << (i + 1) << "' lw 2";
  }
  gp << "\n\n";

  gp << "set terminal pngcairo enhanced font 'Helvetica,12' size 800,600\n";
  gp << "set output '" << png_path << "'\n";
  gp << "replot\n";
  gp.close();

  // gnuplotを実行
  std::string cmd = "gnuplot \"" + gnuplot_script + "\"";
  int ret = std::system(cmd.c_str());
  if (ret == 0) {
    std::cout << "<>        EPS generated: " << eps_path << std::endl;
    std::cout << "<>        PNG generated: " << png_path << std::endl;
  } else {
    std::cerr << "<> !warn! gnuplot execution failed (return code: " << ret << ")" << std::endl;
  }
}

// GetConfigFileList は utils.hpp の DiscoverConfigFilesToml に統合済み

int main(int argc, char* argv[]) {
  using namespace crtbp;
  using namespace utils;

  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv);
  bool skip_wait = args.skip_wait;
  std::string output_tag = args.output_tag;

  // ヘッダー出力 (SALI3dV2スタイル)
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP Trajectory Calculator ver2.0" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // パス設定
  std::string config_base_path = CONFIG_DIR;
  std::string output_base_path = OUTPUT_DIR;
  std::string config_dir = config_base_path + "/trajectory_calc";
  std::string output_dir = output_base_path + "/trajectory_calc";

  // 天文定数の読み込み
  std::string astro_param_file = config_base_path + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);

  const double kGMSUN = astro_params.gm_sun;
  const double kGMEARTH = astro_params.gm_earth;
  const double kMU = kGMEARTH / (kGMEARTH + kGMSUN);

  std::cout << "<>    mu parameter: " << std::setprecision(15) << kMU << std::endl;
  std::cout << "<>" << std::endl;

  // configファイル一覧を取得（_sample付きは除外）
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = args.is_continuous;  // コマンドライン引数で指定
  std::vector<std::string> config_files =
      DiscoverConfigFilesToml(config_dir, "trajectory_calc_config", discovery_opts);
  if (config_files.empty()) {
    std::cerr << "<> !err! No config files found in: " << config_dir << std::endl;
    std::cerr << "<>        Expected pattern: trajectory_calc_config"
              << (args.is_continuous ? "_*.toml" : ".toml") << std::endl;
    return -1;
  }

  std::cout << "<>    Found " << config_files.size() << " config file(s):" << std::endl;
  for (const auto& file : config_files) {
    std::cout << "<>        - " << file << std::endl;
  }
  std::cout << "<>" << std::endl;

  // 出力ディレクトリの作成（セッションディレクトリ）
  OutputDirResult output_result =
      CreateSessionOutputDir(output_base_path, "trajectory_calc", output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string run_output_dir = output_result.session_dir;

  // 全体の実行時間計測
  auto start_total = std::chrono::system_clock::now();

  // 各configファイルについて処理
  for (size_t config_idx = 0; config_idx < config_files.size(); ++config_idx) {
    const std::string& config_filepath = config_files[config_idx];
    std::string config_filename = fs::path(config_filepath).stem().string();

    std::cout << "<>================================================================" << std::endl;
    std::cout << "<>    Processing config file " << (config_idx + 1) << " / " << config_files.size()
              << std::endl;
    std::cout << "<>    File: " << config_filepath << std::endl;

    // 設定ファイル読み込み
    TrajectoryConfig config;
    if (!LoadTrajectoryConfig(config_filepath, &config)) {
      std::cerr << "<> !err! Failed to load config file, skipping..." << std::endl;
      continue;
    }

    if (config.initial_coords.empty()) {
      std::cerr << "<> !err! No COORD entries found in config file, skipping..." << std::endl;
      continue;
    }

    std::cout << "<>        CALC TIMESTEP: " << config.calc_timestep << std::endl;
    std::cout << "<>        TIME THRESHOLD: " << config.time_threshold << std::endl;
    std::cout << "<>        Number of COORDs: " << config.initial_coords.size() << std::endl;
    std::string chaos_index_str;
    switch (config.chaos_index_type) {
      case ChaosIndexType::NONE:
        chaos_index_str = "NONE";
        break;
      case ChaosIndexType::SALI:
        chaos_index_str = "SALI";
        break;
      case ChaosIndexType::GALI:
        chaos_index_str = "GALI" + std::to_string(config.gali_k);
        break;
      case ChaosIndexType::LLE:
        chaos_index_str = "LLE";
        break;
    }
    std::cout << "<>        CHAOS_INDEX: " << chaos_index_str << std::endl;

    // 積分器名を表示
    std::string integrator_str;
    switch (config.integrator_type) {
      case IntegratorType::kSymplectic4th:
        integrator_str = "Symplectic 4th Order";
        break;
      case IntegratorType::kSymplectic6th:
        integrator_str = "Symplectic 6th Order (Yoshida)";
        break;
      case IntegratorType::kDormandPrince:
        integrator_str = "Dormand-Prince 5(4) (DOPRI5)";
        break;
      case IntegratorType::kRungeKutta4th:
        integrator_str = "Runge-Kutta 4th Order (RK4)";
        break;
    }
    std::cout << "<> INTEGRATOR: " << integrator_str << std::endl;
    // DOPRI使用時は追加パラメータを表示
    if (config.integrator_type == IntegratorType::kDormandPrince) {
      std::cout << "<>        DOPRI_ABS_TOL: " << std::scientific << config.dopri_abs_tol
                << std::endl;
      std::cout << "<>        DOPRI_REL_TOL: " << std::scientific << config.dopri_rel_tol
                << std::endl;
      std::cout << "<>        DOPRI_INITIAL_STEP: " << std::fixed << config.dopri_initial_step
                << std::endl;
      std::cout << "<>        DOPRI_MAX_STEP: " << std::fixed << config.dopri_max_step << std::endl;
    }
    std::cout << "<>        FREQ_ANALYSIS: "
              << (config.enable_freq_analysis ? "ENABLED" : "DISABLED") << std::endl;

    // 出力フラグ表示
    std::cout << "<>        OUTPUT_TRAJECTORY: " << (config.output_trajectory ? "ON" : "OFF")
              << std::endl;
    std::cout << "<>        OUTPUT_CHAOS_INDEX: " << (config.output_chaos_index ? "ON" : "OFF")
              << std::endl;
    std::cout << "<>        OUTPUT_ORBITAL_ELEMENTS: "
              << (config.output_orbital_elements ? "ON" : "OFF") << std::endl;
    std::cout << "<>        OUTPUT_FREQ_ANALYSIS: " << (config.output_freq_analysis ? "ON" : "OFF")
              << std::endl;
    std::cout << "<>        OUTPUT_GNUPLOT: " << (config.output_gnuplot ? "ON" : "OFF")
              << std::endl;

    // インパルス設定表示
    if (config.impulse.enabled) {
      std::cout << "<>        IMPULSE: ENABLED" << std::endl;
      if (config.impulse.trigger_type == "distance_from_point") {
        std::cout << "<>          Trigger: distance_from_point" << std::endl;
        std::cout << "<>          Point: (" << config.impulse.trigger_point[0] << ", "
                  << config.impulse.trigger_point[1] << ", " << config.impulse.trigger_point[2]
                  << ")" << std::endl;
        std::cout << "<>          Radius: " << config.impulse.trigger_radius << std::endl;
      } else if (config.impulse.trigger_type == "time") {
        std::cout << "<>          Trigger: time = " << config.impulse.trigger_value << std::endl;
      } else {
        std::cout << "<>          Trigger: " << config.impulse.trigger_type << " ("
                  << config.impulse.trigger_axis << " = " << config.impulse.trigger_value << ", "
                  << config.impulse.trigger_direction << ")" << std::endl;
      }
      std::cout << "<>          Mode: " << config.impulse.velocity_mode << std::endl;
      if (config.impulse.velocity_mode == "delta") {
        std::cout << "<>          Delta-V: (" << config.impulse.delta_v[0] << ", "
                  << config.impulse.delta_v[1] << ", " << config.impulse.delta_v[2] << ")"
                  << std::endl;
      } else {
        std::cout << "<>          Target-V: (" << config.impulse.target_velocity[0] << ", "
                  << config.impulse.target_velocity[1] << ", " << config.impulse.target_velocity[2]
                  << ")" << std::endl;
      }
    } else {
      std::cout << "<>        IMPULSE: DISABLED" << std::endl;
    }

    // 計算のステップ数
    int num_steps = static_cast<int>(config.time_threshold / config.calc_timestep);
    std::cout << "<>    Total integration steps per trajectory: " << num_steps << std::endl;
    std::cout << "<>" << std::endl;

    auto start_config = std::chrono::system_clock::now();

    // configファイルごとに出力サブフォルダを作成（run_output_dir内）
    std::string config_output_dir = run_output_dir + "/" + config_filename;

    if (!fs::exists(config_output_dir)) {
      fs::create_directories(config_output_dir);
      std::cout << "<>    Created config output subdirectory: " << config_output_dir << std::endl;
    }

    // 出力ファイル名（configファイルごとに1つ）
    std::string base_name = config_filename;
    std::string csv_path = config_output_dir + "/" + base_name + ".csv";

    // CSVファイルを開く
    std::ofstream ofs(csv_path);
    if (!ofs) {
      std::cerr << "<> !err! Cannot create output file: " << csv_path << std::endl;
      continue;
    }

    ofs << std::setprecision(15) << std::fixed;

    // ファイルヘッダー書き込み
    ofs << "# Trajectory Calculation Output\n";
    ofs << "# Config file: " << config_filepath << "\n";
    ofs << "# CALC TIMESTEP=" << config.calc_timestep << "\n";
    ofs << "# TIME THRESHOLD=" << config.time_threshold << "\n";
    ofs << "# MU=" << kMU << "\n";
    ofs << "# Number of trajectories: " << config.initial_coords.size() << "\n";
    ofs << "# CHAOS_INDEX=" << chaos_index_str << "\n";
    if (config.chaos_index_type != ChaosIndexType::NONE) {
      ofs << "# Data format: time,x,y,z,vx,vy,vz,jacobi,chaos_index\n";
    } else {
      ofs << "# Data format: time,x,y,z,vx,vy,vz,jacobi\n";
    }
    ofs << "# Trajectories are separated by blank lines for gnuplot index\n";
    ofs << "#\n";

    // ループ前にinitial_coordsをコピー（メモリ破損回避）
    std::vector<my_type::State<double>> safe_initial_coords = config.initial_coords;
    const size_t num_coords = safe_initial_coords.size();

    // 各初期条件について計算
    for (size_t coord_idx = 0; coord_idx < num_coords; ++coord_idx) {
      std::cout << "<>----------------------------------------------------------------"
                << std::endl;
      std::cout << "<>    Calculating trajectory " << (coord_idx + 1) << " / " << num_coords
                << std::endl;

      // safe_initial_coordsから初期状態を取得
      const my_type::State<double>& initial_state = safe_initial_coords[coord_idx];
      my_type::State<double> state = initial_state;

      std::cout << "<>        Initial state: (" << state.x << ", " << state.y << ", " << state.z
                << ", " << state.vx << ", " << state.vy << ", " << state.vz << ")" << std::endl;

      auto start = std::chrono::system_clock::now();

      // 軌道ヘッダー（コメント）
      ofs << "# Trajectory " << (coord_idx + 1) << " Initial: " << state.x << "," << state.y << ","
          << state.z << "," << state.vx << "," << state.vy << "," << state.vz << "\n";

      // 初期状態を記録
      double initial_chaos_value = 2.0;  // SALI/GALI初期値
      double initial_jacobi = crtbp::calc_jacobi_integral(state, kMU);
      if (config.chaos_index_type != ChaosIndexType::NONE) {
        ofs << 0.0 << "," << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
            << state.vy << "," << state.vz << "," << initial_jacobi << "," << initial_chaos_value
            << "\n";
      } else {
        ofs << 0.0 << "," << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
            << state.vy << "," << state.vz << "," << initial_jacobi << "\n";
      }

      // 積分ループ
      double current_time = 0.0;
      int output_interval = std::max(1, num_steps / 1000);  // 最大1000点出力

      // カオス指標時系列データ
      std::vector<std::pair<double, double>> chaos_timeseries;
      if (config.chaos_index_type != ChaosIndexType::NONE) {
        chaos_timeseries.reserve(num_steps / output_interval + 1);
        chaos_timeseries.emplace_back(0.0, initial_chaos_value);
      }

      // 軌道要素時系列データ (time, a, e, i, Omega, omega, nu)
      struct OrbElemTimestep {
        double t, a, e, i, Omega, omega, nu;
      };
      std::vector<OrbElemTimestep> orb_elem_timeseries;
      orb_elem_timeseries.reserve(num_steps / output_interval + 1);
      // 初期軌道要素を計算
      auto init_orb = crtbp::ConvertToOrbitalElements(state, kMU);
      orb_elem_timeseries.push_back(
          {0.0, init_orb.a, init_orb.e, init_orb.i, init_orb.Omega, init_orb.omega, init_orb.nu});

      // 周波数解析用の軌道バッファ
      rtbp::freq_analysis::TrajectoryBuffer<double> freq_buffer;
      freq_buffer.Reserve(num_steps + 1);
      freq_buffer.Push(0.0, state);  // 初期状態を追加

      // カオス指標計算用の状態（タイプによって使い分け）
      crtbp::SaliState<double> sali_state;
      crtbp::GaliState<double, 4> gali4_state;
      crtbp::GaliState<double, 6> gali6_state;

      // 初期化
      if (config.chaos_index_type == ChaosIndexType::SALI ||
          (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 2)) {
        sali_state.state = crtbp::ConvertToCanonical(state);
        sali_state.w1 = crtbp::CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        sali_state.w2 = crtbp::CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
      } else if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 4) {
        gali4_state.state = crtbp::ConvertToCanonical(state);
        gali4_state.w[0] = crtbp::CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        gali4_state.w[1] = crtbp::CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        gali4_state.w[2] = crtbp::CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
        gali4_state.w[3] = crtbp::CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
      } else if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 6) {
        gali6_state.state = crtbp::ConvertToCanonical(state);
        gali6_state.w[0] = crtbp::CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        gali6_state.w[1] = crtbp::CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        gali6_state.w[2] = crtbp::CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
        gali6_state.w[3] = crtbp::CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        gali6_state.w[4] = crtbp::CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
        gali6_state.w[5] = crtbp::CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
      } else if (config.chaos_index_type == ChaosIndexType::LLE) {
        // LLE用に初期化（SALI状態を流用、w1のみ使用）
        sali_state.state = crtbp::ConvertToCanonical(state);
        sali_state.w1 = crtbp::CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        sali_state.w2 = crtbp::CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  // 未使用
      }

      // LLE計算用の変数
      double lle_sum = 0.0;      // log(増大率)の累積和
      int lle_renorm_count = 0;  // 再正規化回数

      // インパルス用変数
      bool impulse_applied = false;
      my_type::State<double> prev_state = state;  // 前ステップの状態保持

      for (int step = 0; step < num_steps; ++step) {
        double chaos_value = 0.0;
        double prev_time = current_time;  // インパルス検出用に時刻を保存

        switch (config.chaos_index_type) {
          case ChaosIndexType::SALI:
          case ChaosIndexType::GALI:
            if (config.gali_k == 2 || config.chaos_index_type == ChaosIndexType::SALI) {
              // SALI (K=2)
              crtbp::SymplecticStep6thOrderSALI(kMU, &sali_state, config.calc_timestep);
              sali_state.w1.Normalize();
              sali_state.w2.Normalize();
              state.x = sali_state.state.qx;
              state.y = sali_state.state.qy;
              state.z = sali_state.state.qz;
              state.vx = sali_state.state.px + sali_state.state.qy;
              state.vy = sali_state.state.py - sali_state.state.qx;
              state.vz = sali_state.state.pz;
              double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
              double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
              chaos_value = std::min(norm_plus, norm_minus);
            } else if (config.gali_k == 4) {
              // GALI4
              crtbp::SymplecticStep6thOrderGALI(kMU, &gali4_state, config.calc_timestep);
              gali4_state.NormalizeDeviationVectors();
              state.x = gali4_state.state.qx;
              state.y = gali4_state.state.qy;
              state.z = gali4_state.state.qz;
              state.vx = gali4_state.state.px + gali4_state.state.qy;
              state.vy = gali4_state.state.py - gali4_state.state.qx;
              state.vz = gali4_state.state.pz;
              chaos_value = gali4_state.ComputeGALI();
            } else if (config.gali_k == 6) {
              // GALI6
              crtbp::SymplecticStep6thOrderGALI(kMU, &gali6_state, config.calc_timestep);
              gali6_state.NormalizeDeviationVectors();
              state.x = gali6_state.state.qx;
              state.y = gali6_state.state.qy;
              state.z = gali6_state.state.qz;
              state.vx = gali6_state.state.px + gali6_state.state.qy;
              state.vy = gali6_state.state.py - gali6_state.state.qx;
              state.vz = gali6_state.state.pz;
              chaos_value = gali6_state.ComputeGALI();
            }
            break;

          case ChaosIndexType::LLE:
            // LLE計算（全積分器対応）
            switch (config.integrator_type) {
              case IntegratorType::kSymplectic4th:
                crtbp::SymplecticStep4thOrderSALI(kMU, &sali_state, config.calc_timestep);
                break;
              case IntegratorType::kSymplectic6th:
                crtbp::SymplecticStep6thOrderSALI(kMU, &sali_state, config.calc_timestep);
                break;
              case IntegratorType::kDormandPrince: {
                // DOPRI5での偏差ベクトル伝搬（固定ステップ）
                // 注意:
                // SaliStateは可変ステップ制御に必要なnorm_inf/abs操作をサポートしないため固定ステップ
                using namespace boost::numeric::odeint;
                crtbp::SaliOrbitSystem<double> sali_system(kMU);
                // vector_space_algebraを明示的に指定（SaliStateはRangeではない）
                runge_kutta_dopri5<crtbp::SaliState<double>, double, crtbp::SaliState<double>,
                                   double, vector_space_algebra>
                    stepper;
                stepper.do_step(sali_system, sali_state, 0.0, config.calc_timestep);
              } break;
              case IntegratorType::kRungeKutta4th: {
                // RK4での偏差ベクトル伝搬
                using namespace boost::numeric::odeint;
                crtbp::SaliOrbitSystem<double> sali_system(kMU);
                // vector_space_algebraを明示的に指定（SaliStateはRangeではない）
                runge_kutta4<crtbp::SaliState<double>, double, crtbp::SaliState<double>, double,
                             vector_space_algebra>
                    stepper;
                stepper.do_step(sali_system, sali_state, 0.0, config.calc_timestep);
              } break;
            }

            // 定期的に再正規化し、増大率を記録
            if ((step + 1) % config.lle_renorm_interval == 0 && step >= config.lle_skip_steps) {
              double norm_before = sali_state.w1.Norm();
              if (norm_before > 1e-12) {  // ゼロ除算回避
                lle_sum += std::log(norm_before);
                lle_renorm_count++;
              }
              sali_state.w1.Normalize();
            }

            // 主軌道の状態を更新
            state.x = sali_state.state.qx;
            state.y = sali_state.state.qy;
            state.z = sali_state.state.qz;
            state.vx = sali_state.state.px + sali_state.state.qy;
            state.vy = sali_state.state.py - sali_state.state.qx;
            state.vz = sali_state.state.pz;

            // 現在のLLE推定値を計算（時系列出力用）
            if (lle_renorm_count > 0) {
              double elapsed_time = (step + 1) * config.calc_timestep;
              chaos_value = lle_sum / elapsed_time;
            } else {
              chaos_value = 0.0;
            }
            break;

          case ChaosIndexType::NONE:
          default:
            // 積分器の選択
            switch (config.integrator_type) {
              case IntegratorType::kSymplectic4th:
                state = crtbp::SymplecticStep4thOrder(kMU, state, config.calc_timestep);
                break;
              case IntegratorType::kSymplectic6th:
                state = crtbp::SymplecticStep6thOrder(kMU, state, config.calc_timestep);
                break;
              case IntegratorType::kDormandPrince: {
                // Boost odeint を使用したDormand-Prince法（可変ステップ）
                using namespace boost::numeric::odeint;
                using StateVec = std::array<double, 6>;
                StateVec y = {state.x, state.y, state.z, state.vx, state.vy, state.vz};
                auto crtbp_eom = [mu = kMU](const StateVec& y, StateVec& dydt, double /*t*/) {
                  double r1_3 =
                      std::pow((y[0] + mu) * (y[0] + mu) + y[1] * y[1] + y[2] * y[2], 1.5);
                  double r2_3 =
                      std::pow((y[0] - 1 + mu) * (y[0] - 1 + mu) + y[1] * y[1] + y[2] * y[2], 1.5);
                  dydt[0] = y[3];
                  dydt[1] = y[4];
                  dydt[2] = y[5];
                  dydt[3] =
                      2 * y[4] + y[0] - (1 - mu) * (y[0] + mu) / r1_3 - mu * (y[0] - 1 + mu) / r2_3;
                  dydt[4] = -2 * y[3] + y[1] - (1 - mu) * y[1] / r1_3 - mu * y[1] / r2_3;
                  dydt[5] = -(1 - mu) * y[2] / r1_3 - mu * y[2] / r2_3;
                };
                // 誤差制御付きコントローラを構築
                auto controlled_stepper =
                    make_controlled(config.dopri_abs_tol, config.dopri_rel_tol,
                                    config.dopri_max_step, runge_kutta_dopri5<StateVec>());
                // 1ステップ分(calc_timestepまで)を適応的に積分
                double t_local = 0.0;
                double dt = config.dopri_initial_step;
                while (t_local < config.calc_timestep) {
                  double remaining = config.calc_timestep - t_local;
                  if (dt > remaining) {
                    dt = remaining;
                  }
                  controlled_step_result result =
                      controlled_stepper.try_step(crtbp_eom, y, t_local, dt);
                  if (result == success) {
                    // ステップ成功: t_localは更新済み
                  }
                  // failed: dtが自動調整され、再試行される
                }
                state = {y[0], y[1], y[2], y[3], y[4], y[5]};
              } break;
              case IntegratorType::kRungeKutta4th: {
                // 4次ルンゲクッタ法
                using StateVec = std::array<double, 6>;
                StateVec y = {state.x, state.y, state.z, state.vx, state.vy, state.vz};
                double h = config.calc_timestep;
                double mu = kMU;
                auto calc_deriv = [mu](const StateVec& s) -> StateVec {
                  double r1_3 =
                      std::pow((s[0] + mu) * (s[0] + mu) + s[1] * s[1] + s[2] * s[2], 1.5);
                  double r2_3 =
                      std::pow((s[0] - 1 + mu) * (s[0] - 1 + mu) + s[1] * s[1] + s[2] * s[2], 1.5);
                  return {
                      s[3],
                      s[4],
                      s[5],
                      2 * s[4] + s[0] - (1 - mu) * (s[0] + mu) / r1_3 - mu * (s[0] - 1 + mu) / r2_3,
                      -2 * s[3] + s[1] - (1 - mu) * s[1] / r1_3 - mu * s[1] / r2_3,
                      -(1 - mu) * s[2] / r1_3 - mu * s[2] / r2_3};
                };
                StateVec k1 = calc_deriv(y);
                StateVec y2, y3, y4;
                for (int i = 0; i < 6; ++i) y2[i] = y[i] + 0.5 * h * k1[i];
                StateVec k2 = calc_deriv(y2);
                for (int i = 0; i < 6; ++i) y3[i] = y[i] + 0.5 * h * k2[i];
                StateVec k3 = calc_deriv(y3);
                for (int i = 0; i < 6; ++i) y4[i] = y[i] + h * k3[i];
                StateVec k4 = calc_deriv(y4);
                for (int i = 0; i < 6; ++i) {
                  y[i] += h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
                }
                state = {y[0], y[1], y[2], y[3], y[4], y[5]};
              } break;
            }
            break;
        }
        current_time += config.calc_timestep;

        // ===== インパルス適用チェック =====
        if (config.impulse.enabled && !impulse_applied) {
          if (CheckImpulseTrigger(prev_state, state, prev_time, current_time, config.impulse,
                                  kMU)) {
            // 線形補間で正確なトリガー位置を計算
            double alpha = ComputeTriggerFraction(prev_state, state, prev_time, current_time,
                                                  config.impulse, kMU);
            double trigger_time = prev_time + alpha * config.calc_timestep;

            // トリガー位置での状態を補間（オプション: 精度を上げたい場合）
            my_type::State<double> trigger_state = InterpolateState(prev_state, state, alpha);

            // インパルス適用前のヤコビ積分
            double pre_jacobi = crtbp::calc_jacobi_integral(trigger_state, kMU);

            // 速度変更を適用
            ApplyImpulse(&trigger_state, config.impulse);

            // インパルス適用後のヤコビ積分
            double post_jacobi = crtbp::calc_jacobi_integral(trigger_state, kMU);

            impulse_applied = true;

            // 状態を更新（インパルス後の状態から継続）
            state = trigger_state;

            std::cout << "<>        [IMPULSE] Applied at t=" << std::setprecision(6) << trigger_time
                      << std::endl;
            if (config.impulse.trigger_type == "distance_from_point") {
              std::cout << "<>                  Trigger: distance from ("
                        << config.impulse.trigger_point[0] << ", "
                        << config.impulse.trigger_point[1] << ", "
                        << config.impulse.trigger_point[2]
                        << ") <= " << config.impulse.trigger_radius << std::endl;
            } else if (config.impulse.trigger_type == "time") {
              std::cout << "<>                  Trigger: time = " << config.impulse.trigger_value
                        << std::endl;
            } else {
              std::cout << "<>                  Trigger: " << config.impulse.trigger_axis << " = "
                        << config.impulse.trigger_value << std::endl;
            }
            std::cout << "<>                  Mode: " << config.impulse.velocity_mode << std::endl;
            if (config.impulse.velocity_mode == "delta") {
              std::cout << "<>                  Delta-V: (" << config.impulse.delta_v[0] << ", "
                        << config.impulse.delta_v[1] << ", " << config.impulse.delta_v[2] << ")"
                        << std::endl;
            } else {
              std::cout << "<>                  Target V: (" << config.impulse.target_velocity[0]
                        << ", " << config.impulse.target_velocity[1] << ", "
                        << config.impulse.target_velocity[2] << ")" << std::endl;
            }
            std::cout << "<>                  Jacobi: " << pre_jacobi << " -> " << post_jacobi
                      << " (delta=" << post_jacobi - pre_jacobi << ")" << std::endl;

            // カオス指標計算用の状態も更新が必要な場合（SALI/GALI使用時）
            if (config.chaos_index_type == ChaosIndexType::SALI ||
                config.chaos_index_type == ChaosIndexType::GALI ||
                config.chaos_index_type == ChaosIndexType::LLE) {
              sali_state.state = crtbp::ConvertToCanonical(state);
            }
            if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 4) {
              gali4_state.state = crtbp::ConvertToCanonical(state);
            }
            if (config.chaos_index_type == ChaosIndexType::GALI && config.gali_k == 6) {
              gali6_state.state = crtbp::ConvertToCanonical(state);
            }
          }
        }

        // 前ステップの状態を更新
        prev_state = state;

        // 周波数解析用にデータを収集
        freq_buffer.Push(current_time, state);

        // 間引いて出力
        if ((step + 1) % output_interval == 0 || step == num_steps - 1) {
          double jacobi_value = crtbp::calc_jacobi_integral(state, kMU);
          if (config.chaos_index_type != ChaosIndexType::NONE) {
            ofs << current_time << "," << state.x << "," << state.y << "," << state.z << ","
                << state.vx << "," << state.vy << "," << state.vz << "," << jacobi_value << ","
                << chaos_value << "\n";
            chaos_timeseries.emplace_back(current_time, chaos_value);
          } else {
            ofs << current_time << "," << state.x << "," << state.y << "," << state.z << ","
                << state.vx << "," << state.vy << "," << state.vz << "," << jacobi_value << "\n";
          }
          // 軌道要素を計算して保存
          auto orb = crtbp::ConvertToOrbitalElements(state, kMU);
          orb_elem_timeseries.push_back(
              {current_time, orb.a, orb.e, orb.i, orb.Omega, orb.omega, orb.nu});
        }

        // 進捗表示
        if (step % (num_steps / 10 + 1) == 0) {
          double progress = static_cast<double>(step) / num_steps * 100.0;
          std::cout << "\r<>        Progress: " << std::fixed << std::setprecision(1) << progress
                    << "%" << std::flush;
        }
      }
      std::cout << "\r<>        Progress: 100.0%" << std::endl;

      // カオス指標時系列プロット生成
      if (config.output_chaos_index && config.chaos_index_type != ChaosIndexType::NONE &&
          !chaos_timeseries.empty()) {
        std::string chaos_csv_path = config_output_dir + "/" + base_name + "_" + chaos_index_str +
                                     "_traj" + std::to_string(coord_idx + 1) + ".csv";
        std::ofstream chaos_ofs(chaos_csv_path);
        if (chaos_ofs) {
          chaos_ofs << std::setprecision(15) << std::fixed;
          chaos_ofs << "# " << chaos_index_str << " Time Series - Trajectory " << (coord_idx + 1)
                    << "\n";
          chaos_ofs << "time," << chaos_index_str << "\n";
          for (const auto& ts : chaos_timeseries) {
            chaos_ofs << ts.first << "," << ts.second << "\n";
          }
          chaos_ofs.close();
          std::cout << "<>        " << chaos_index_str << " CSV: " << chaos_csv_path << std::endl;

          // gnuplotでカオス指標時系列プロット
          if (config.output_gnuplot) {
            std::string chaos_gp_path = config_output_dir + "/" + base_name + "_" +
                                        chaos_index_str + "_traj" + std::to_string(coord_idx + 1) +
                                        ".gp";
            std::string chaos_eps_path = config_output_dir + "/" + base_name + "_" +
                                         chaos_index_str + "_traj" + std::to_string(coord_idx + 1) +
                                         ".eps";
            std::string chaos_png_path = config_output_dir + "/" + base_name + "_" +
                                         chaos_index_str + "_traj" + std::to_string(coord_idx + 1) +
                                         ".png";
            std::ofstream gp(chaos_gp_path);
            if (gp) {
              gp << "set datafile separator ','\n";
              gp << "set datafile commentschars '#'\n";
              gp << "set terminal postscript eps enhanced color font 'Helvetica,14'\n";
              gp << "set output '" << chaos_eps_path << "'\n";
              gp << "set xlabel 'Time (non-dim)'\n";
              gp << "set ylabel '" << chaos_index_str << "'\n";
              gp << "set title '" << chaos_index_str << " Time Series - Trajectory "
                 << (coord_idx + 1) << "'\n";
              gp << "set grid\n";
              gp << "set logscale xy\n";
              gp << "set format x '10^{%L}'\n";
              gp << "set format y '10^{%L}'\n";
              gp << "plot '" << chaos_csv_path
                 << "' using ($1 > 0 ? $1 : 1/0):($2 > 0 ? $2 : 1e-16) with lines lw 2 lc rgb "
                    "'blue' "
                    "title '"
                 << chaos_index_str << "'\n";
              gp << "\n";
              gp << "set terminal pngcairo enhanced font 'Helvetica,14' size 1200,600\n";
              gp << "set output '" << chaos_png_path << "'\n";
              gp << "replot\n";
              gp.close();

              std::string cmd = "gnuplot \"" + chaos_gp_path + "\"";
              int ret = std::system(cmd.c_str());
              if (ret == 0) {
                std::cout << "<>        " << chaos_index_str << " plot: " << chaos_png_path
                          << std::endl;
              }
            }
          }
        }
      }

      // 軌道要素時系列CSV出力
      if (config.output_orbital_elements && !orb_elem_timeseries.empty()) {
        std::string orb_csv_path = config_output_dir + "/" + base_name + "_orb_elem_traj" +
                                   std::to_string(coord_idx + 1) + ".csv";
        std::ofstream orb_ofs(orb_csv_path);
        if (orb_ofs) {
          orb_ofs << std::setprecision(15) << std::fixed;
          orb_ofs << "# Osculating Orbital Elements Time Series - Trajectory " << (coord_idx + 1)
                  << "\n";
          orb_ofs << "# Elements are computed relative to m2 (secondary body)\n";
          orb_ofs << "# a: semi-major axis (non-dim), e: eccentricity, i: inclination (rad)\n";
          orb_ofs << "# Omega: RAAN (rad), omega: argument of periapsis (rad), nu: true anomaly "
                     "(rad)\n";
          orb_ofs << "time,a,e,i,Omega,omega,nu\n";
          for (const auto& elem : orb_elem_timeseries) {
            orb_ofs << elem.t << "," << elem.a << "," << elem.e << "," << elem.i << ","
                    << elem.Omega << "," << elem.omega << "," << elem.nu << "\n";
          }
          orb_ofs.close();
          std::cout << "<>        Orbital elements CSV: " << orb_csv_path << std::endl;

          // gnuplotで軌道要素時系列プロット (a, e, i)
          if (config.output_gnuplot) {
            std::string orb_gp_path = config_output_dir + "/" + base_name + "_orb_elem_traj" +
                                      std::to_string(coord_idx + 1) + ".gp";
            std::string orb_png_path = config_output_dir + "/" + base_name + "_orb_elem_traj" +
                                       std::to_string(coord_idx + 1) + ".png";
            std::ofstream gp_orb(orb_gp_path);
            if (gp_orb) {
              gp_orb << "set datafile separator ','\n";
              gp_orb << "set datafile commentschars '#'\n";
              gp_orb << "set terminal pngcairo enhanced font 'Helvetica,12' size 1200,900\n";
              gp_orb << "set output '" << orb_png_path << "'\n";
              gp_orb << "set multiplot layout 3,2 title 'Orbital Elements - Trajectory "
                     << (coord_idx + 1) << "'\n";
              gp_orb << "set xlabel 'Time (non-dim)'\n";
              gp_orb << "set grid\n";
              // Semi-major axis
              gp_orb << "set ylabel 'a (semi-major axis)'\n";
              gp_orb << "plot '" << orb_csv_path
                     << "' using 1:2 with lines lw 1.5 lc rgb 'blue' "
                        "notitle\n";
              // Eccentricity
              gp_orb << "set ylabel 'e (eccentricity)'\n";
              gp_orb << "plot '" << orb_csv_path
                     << "' using 1:3 with lines lw 1.5 lc rgb 'red' "
                        "notitle\n";
              // Inclination
              gp_orb << "set ylabel 'i (inclination) [rad]'\n";
              gp_orb << "plot '" << orb_csv_path
                     << "' using 1:4 with lines lw 1.5 lc rgb 'green' "
                        "notitle\n";
              // RAAN
              gp_orb << "set ylabel 'Omega (RAAN) [rad]'\n";
              gp_orb << "plot '" << orb_csv_path
                     << "' using 1:5 with lines lw 1.5 lc rgb 'purple' "
                        "notitle\n";
              // Argument of periapsis
              gp_orb << "set ylabel 'omega (arg of peri) [rad]'\n";
              gp_orb << "plot '" << orb_csv_path
                     << "' using 1:6 with lines lw 1.5 lc rgb 'orange' "
                        "notitle\n";
              // True anomaly
              gp_orb << "set ylabel 'nu (true anomaly) [rad]'\n";
              gp_orb << "plot '" << orb_csv_path
                     << "' using 1:7 with lines lw 1.5 lc rgb 'cyan' "
                        "notitle\n";
              gp_orb << "unset multiplot\n";
              gp_orb.close();

              std::string cmd = "gnuplot \"" + orb_gp_path + "\"";
              int ret = std::system(cmd.c_str());
              if (ret == 0) {
                std::cout << "<>        Orbital elements plot: " << orb_png_path << std::endl;
              }
            }
          }
        }
      }

      // 軌道間に空行を挿入（gnuplotのインデックス機能用）
      // 最後の軌道の後には空行を入れない
      if (coord_idx < config.initial_coords.size() - 1) {
        ofs << "\n\n";  // 2つの空行でgnuplotのインデックスを区切る
      }

      // ========== 周波数解析 (Laskar's Frequency Map Analysis) ==========
      if (config.enable_freq_analysis && config.output_freq_analysis) {
        if (freq_buffer.Size() >= 100) {  // 十分なデータがある場合のみ解析
          std::cout << "<>        Performing Frequency Analysis..." << std::endl;

          // x, y, z成分の周波数解析
          auto freq_result_x = rtbp::freq_analysis::AnalyzeFromBuffer(freq_buffer, "x");
          auto freq_result_y = rtbp::freq_analysis::AnalyzeFromBuffer(freq_buffer, "y");
          auto freq_result_r = rtbp::freq_analysis::AnalyzeFromBuffer(freq_buffer, "r");

          // 結果をコンソールに出力
          std::cout << "<>        [Frequency Analysis Result]" << std::endl;
          std::cout << "<>          X-component: nu1=" << std::scientific << std::setprecision(6)
                    << freq_result_x.nu1 << ", nu2=" << freq_result_x.nu2
                    << ", log10(D)=" << std::fixed << std::setprecision(3) << freq_result_x.log10_D
                    << std::endl;
          std::cout << "<>          Y-component: nu1=" << std::scientific << std::setprecision(6)
                    << freq_result_y.nu1 << ", nu2=" << freq_result_y.nu2
                    << ", log10(D)=" << std::fixed << std::setprecision(3) << freq_result_y.log10_D
                    << std::endl;
          std::cout << "<>          R-component: nu1=" << std::scientific << std::setprecision(6)
                    << freq_result_r.nu1 << ", nu2=" << freq_result_r.nu2
                    << ", log10(D)=" << std::fixed << std::setprecision(3) << freq_result_r.log10_D
                    << std::endl;
          std::cout << "<>          Chaos level (X): "
                    << rtbp::freq_analysis::GetChaosLevelString(freq_result_x.log10_D) << std::endl;

          // 周波数解析結果をCSVに出力
          std::string freq_csv_path = config_output_dir + "/" + base_name + "_freq_analysis_traj" +
                                      std::to_string(coord_idx + 1) + ".csv";
          std::ofstream freq_ofs(freq_csv_path);
          if (freq_ofs) {
            freq_ofs << std::setprecision(15) << std::fixed;
            freq_ofs << "# Frequency Map Analysis Result - Trajectory " << (coord_idx + 1) << "\n";
            freq_ofs << "# Method: Laskar's Frequency Map Analysis\n";
            freq_ofs << "# Integration time: " << config.time_threshold << "\n";
            freq_ofs << "# Timestep: " << config.calc_timestep << "\n";
            freq_ofs << "# Number of data points: " << freq_buffer.Size() << "\n";
            freq_ofs << "component,nu1,nu2,diffusion_D,log10_D,is_chaotic,chaos_level\n";
            freq_ofs << "x," << freq_result_x.nu1 << "," << freq_result_x.nu2 << ","
                     << freq_result_x.diffusion_D << "," << freq_result_x.log10_D << ","
                     << (freq_result_x.is_chaotic ? "1" : "0") << ","
                     << rtbp::freq_analysis::GetChaosLevelString(freq_result_x.log10_D) << "\n";
            freq_ofs << "y," << freq_result_y.nu1 << "," << freq_result_y.nu2 << ","
                     << freq_result_y.diffusion_D << "," << freq_result_y.log10_D << ","
                     << (freq_result_y.is_chaotic ? "1" : "0") << ","
                     << rtbp::freq_analysis::GetChaosLevelString(freq_result_y.log10_D) << "\n";
            freq_ofs << "r," << freq_result_r.nu1 << "," << freq_result_r.nu2 << ","
                     << freq_result_r.diffusion_D << "," << freq_result_r.log10_D << ","
                     << (freq_result_r.is_chaotic ? "1" : "0") << ","
                     << rtbp::freq_analysis::GetChaosLevelString(freq_result_r.log10_D) << "\n";
            freq_ofs.close();
            std::cout << "<>        Frequency analysis CSV: " << freq_csv_path << std::endl;
          }
        } else {
          std::cout << "<>        Skipping frequency analysis (insufficient data points: "
                    << freq_buffer.Size() << ")" << std::endl;
        }
      }

      // 経過時間
      auto end = std::chrono::system_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      auto msec = duration.count() % 1000;
      auto sec = duration.count() / 1000 % 60;
      auto min = duration.count() / 1000 / 60;
      std::cout << "<>        Elapsed time: " << min << "m " << sec << "s " << msec << "ms"
                << std::endl;
    }

    ofs.close();
    std::cout << "<>" << std::endl;
    std::cout << "<>    CSV output: " << csv_path << std::endl;

    // gnuplotでプロット生成
    GenerateGnuplot(csv_path, config_output_dir, base_name, static_cast<int>(num_coords));

    // configファイルごとの経過時間
    auto end_config = std::chrono::system_clock::now();
    auto duration_config =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_config - start_config);
    auto msec = duration_config.count() % 1000;
    auto sec = duration_config.count() / 1000 % 60;
    auto min = duration_config.count() / 1000 / 60;
    std::cout << "<>    Config file elapsed time: " << min << "m " << sec << "s " << msec << "ms"
              << std::endl;
  }

  // 全体の経過時間
  auto end_total = std::chrono::system_clock::now();
  auto duration_total =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
  auto msec = duration_total.count() % 1000;
  auto sec = duration_total.count() / 1000 % 60;
  auto min = duration_total.count() / 1000 / 60;

  std::cout << "<>" << std::endl;
  std::cout << "<>================================================================" << std::endl;
  std::cout << "<>    All calculations finished" << std::endl;
  std::cout << "<>    Total elapsed time: " << min << "m " << sec << "s " << msec << "ms"
            << std::endl;
  std::cout << "<>================================================================" << std::endl;

  return 0;
}
