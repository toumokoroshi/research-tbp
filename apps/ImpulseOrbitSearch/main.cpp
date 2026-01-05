/**
 * @file main.cpp
 * @brief インパルス軌道探索アプリケーション
 *
 * 小惑星軌道上の各点からインパルスを与え、目標ヤコビ積分を達成しつつ
 * 地球近傍に長時間滞在する安定軌道を探索する。
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
#include <map>
#include <regex>
#include <rtbp.hpp>
#include <set>
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
 * @brief 設定パラメータ構造体
 */
struct ImpulseOrbitSearchConfig {
  // データ設定
  std::string orbit_data_dir;                               ///< 軌道データディレクトリ
  std::string orbit_file_pattern = "OBT_\\d+_Earth\\.txt";  ///< ファイル名正規表現パターン
  int phase_step = 1;                                       ///< 軌道データの間引きステップ

  // 目標設定
  double target_jacobi = 3.0008;  ///< 目標ヤコビ積分

  // 速度掃引設定（適応的サンプリング）
  bool enable_adaptive_sampling = true;  ///< 適応的サンプリング使用
  int theta_divisions_coarse = 12;       ///< 粗いθ分割数
  int phi_divisions_coarse = 6;          ///< 粗いφ分割数
  int theta_divisions_fine = 6;          ///< 細かいθ分割数
  int phi_divisions_fine = 6;            ///< 細かいφ分割数

  // 細分化条件
  bool refine_high_sali = true;        ///< 高SALI領域を細分化
  double sali_refine_threshold = 0.1;  ///< SALI細分化閾値
  bool refine_no_escape = true;        ///< 脱出なし条件
  bool refine_no_collision = true;     ///< 衝突なし条件

  // Δvフィルタ
  bool enable_dv_filter = false;  ///< Δvフィルタを使用
  double dv_max_threshold = 0.1;  ///< 最大Δv閾値

  // シミュレーション設定
  double calc_duration_nd = 10.0;     ///< 計算時間（無次元）
  double calc_timestep_nd = 0.00001;  ///< タイムステップ（無次元）
  double escape_dist = 0.03;          ///< 脱出判定距離
  double collision_dist = 1e-6;       ///< 衝突判定距離

  // カオス指標
  ChaosIndexType chaos_index_type = ChaosIndexType::SALI;
  int gali_k = 2;

  // 出力設定
  bool output_stable_only = false;  ///< 安定軌道のみ出力
  bool output_trajectory = false;   ///< 軌道時刻歴出力
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
 * @brief 探索結果1件分
 */
struct SearchResult {
  int orbit_idx;            ///< 軌道点インデックス
  double time_j2000;        ///< 元軌道時刻
  double x, y, z;           ///< 初期位置
  double vx, vy, vz;        ///< インパルス後速度
  double v_theta;           ///< 速度方位角(deg)
  double v_phi;             ///< 速度仰角(deg)
  double jacobi;            ///< ヤコビ積分
  double dv_x, dv_y, dv_z;  ///< Δvベクトル
  double dv_mag;            ///< Δv量
  double sali;              ///< SALI値
  int escape_flag;          ///< 脱出フラグ
  int collision_flag;       ///< 衝突フラグ
  double residence_time;    ///< 滞在時間
};

/**
 * @brief 軌道データファイルを読み込む
 */
bool LoadOrbitData(const std::string& filepath, int phase_step, std::vector<OrbitDataRow>* data) {
  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "<> !err! Cannot open orbit data file: " << filepath << std::endl;
    return false;
  }

  std::string line;
  int consecutive_empty_count = 0;
  bool second_section_started = false;
  int line_count = 0;

  while (std::getline(ifs, line)) {
    bool is_empty = line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos;

    if (is_empty) {
      consecutive_empty_count++;
      continue;
    }

    if (consecutive_empty_count >= 2 && !second_section_started) {
      second_section_started = true;
    }
    consecutive_empty_count = 0;

    if (!second_section_started) {
      continue;
    }

    if (line_count % phase_step != 0) {
      line_count++;
      continue;
    }
    line_count++;

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
 * @brief ディレクトリから軌道ファイルを発見
 * @param data_dir ディレクトリパス
 * @param pattern_str ファイル名の正規表現パターン
 * @return マッチしたファイルパスのリスト
 */
std::vector<std::string> DiscoverOrbitFiles(const std::string& data_dir,
                                            const std::string& pattern_str) {
  std::vector<std::string> files;
  if (!fs::exists(data_dir)) {
    std::cerr << "<> !err! DATA_DIR does not exist: " << data_dir << std::endl;
    return files;
  }

  std::regex pattern(pattern_str);
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
 * @brief 球座標から直交座標への変換
 * @param theta 方位角 (rad)
 * @param phi 仰角 (rad)
 * @param magnitude 大きさ
 */
State3d<double> SphericalToCartesian(double theta, double phi, double magnitude) {
  double cos_phi = std::cos(phi);
  return {magnitude * cos_phi * std::cos(theta), magnitude * cos_phi * std::sin(theta),
          magnitude * std::sin(phi)};
}

/**
 * @brief TOML設定ファイルを読み込む
 */
bool LoadConfig(const std::string& filepath, ImpulseOrbitSearchConfig* config) {
  try {
    utils::TomlConfigParser toml(filepath);

    config->orbit_data_dir = toml.GetString("data.orbit_data_dir", "");
    config->orbit_file_pattern = toml.GetString("data.orbit_file_pattern", "OBT_\\d+_Earth\\.txt");
    config->phase_step = toml.GetInt("data.phase_step", 1);

    config->target_jacobi = toml.GetDouble("target.jacobi_integral", 3.0008);

    config->enable_adaptive_sampling =
        toml.GetBool("velocity_sweep.enable_adaptive_sampling", true);
    config->theta_divisions_coarse = toml.GetInt("velocity_sweep.theta_divisions_coarse", 12);
    config->phi_divisions_coarse = toml.GetInt("velocity_sweep.phi_divisions_coarse", 6);
    config->theta_divisions_fine = toml.GetInt("velocity_sweep.theta_divisions_fine", 6);
    config->phi_divisions_fine = toml.GetInt("velocity_sweep.phi_divisions_fine", 6);

    config->refine_high_sali = toml.GetBool("velocity_sweep.refine_high_sali", true);
    config->sali_refine_threshold = toml.GetDouble("velocity_sweep.sali_refine_threshold", 0.1);
    config->refine_no_escape = toml.GetBool("velocity_sweep.refine_no_escape", true);
    config->refine_no_collision = toml.GetBool("velocity_sweep.refine_no_collision", true);

    config->enable_dv_filter = toml.GetBool("velocity_sweep.enable_dv_filter", false);
    config->dv_max_threshold = toml.GetDouble("velocity_sweep.dv_max_threshold", 0.1);

    config->calc_duration_nd = toml.GetDouble("simulation.calc_duration_nd", 10.0);
    config->calc_timestep_nd = toml.GetDouble("simulation.calc_timestep_nd", 0.00001);
    config->escape_dist = toml.GetDouble("simulation.escape_dist", 0.03);
    config->collision_dist = toml.GetDouble("simulation.collision_dist", 1e-6);

    std::string chaos_str = toml.GetString("chaos.index_type", "SALI");
    if (chaos_str == "SALI" || chaos_str == "sali") {
      config->chaos_index_type = ChaosIndexType::SALI;
      config->gali_k = 2;
    } else if (chaos_str == "GALI4" || chaos_str == "gali4") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 4;
    } else if (chaos_str == "GALI6" || chaos_str == "gali6") {
      config->chaos_index_type = ChaosIndexType::GALI;
      config->gali_k = 6;
    }

    config->output_stable_only = toml.GetBool("output.output_stable_only", false);
    config->output_trajectory = toml.GetBool("output.output_trajectory", false);

  } catch (const std::exception& e) {
    std::cerr << "<> !err! Cannot load config: " << e.what() << std::endl;
    return false;
  }
  return true;
}

/**
 * @brief 単一の初期条件でSALI計算を実行
 */
SearchResult RunSaliCalculation(const State<double>& initial_state,
                                const State<double>& original_velocity, int orbit_idx,
                                double time_j2000, double v_theta_deg, double v_phi_deg, double mu,
                                const ImpulseOrbitSearchConfig& config) {
  SearchResult result;
  result.orbit_idx = orbit_idx;
  result.time_j2000 = time_j2000;
  result.x = initial_state.x;
  result.y = initial_state.y;
  result.z = initial_state.z;
  result.vx = initial_state.vx;
  result.vy = initial_state.vy;
  result.vz = initial_state.vz;
  result.v_theta = v_theta_deg;
  result.v_phi = v_phi_deg;
  result.jacobi = calc_jacobi_integral(initial_state, mu);

  // Δv計算
  result.dv_x = initial_state.vx - original_velocity.vx;
  result.dv_y = initial_state.vy - original_velocity.vy;
  result.dv_z = initial_state.vz - original_velocity.vz;
  result.dv_mag =
      std::sqrt(result.dv_x * result.dv_x + result.dv_y * result.dv_y + result.dv_z * result.dv_z);

  result.escape_flag = 0;
  result.collision_flag = 0;
  result.residence_time = config.calc_duration_nd;
  result.sali = -1.0;

  // SALI状態初期化
  CanonicalState<double> canonical_state = ConvertToCanonical(initial_state);
  SaliState<double> sali_state;
  sali_state.state = canonical_state;
  sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  int num_steps = static_cast<int>(config.calc_duration_nd / config.calc_timestep_nd);

  for (int step = 0; step < num_steps; ++step) {
    SymplecticStep6thOrderSALI(mu, &sali_state, config.calc_timestep_nd);
    sali_state.w1.Normalize();
    sali_state.w2.Normalize();

    double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
    double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
    result.sali = std::min(norm_plus, norm_minus);

    double current_time = (step + 1) * config.calc_timestep_nd;

    // 衝突・脱出判定
    double r2 = calc_r2(sali_state.state.qx, sali_state.state.qy, sali_state.state.qz, mu);
    if (r2 < config.collision_dist) {
      result.collision_flag = 1;
      result.residence_time = current_time;
      break;
    }
    if (r2 > config.escape_dist) {
      result.escape_flag = 1;
      result.residence_time = current_time;
      break;
    }
  }

  return result;
}

int main(int argc, char* argv[]) {
  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv);
  bool skip_wait = args.skip_wait;
  std::string output_tag = args.output_tag;

  // ヘッダー出力
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            Impulse Orbit Search Application ver1.0" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // パス設定
  const std::string kConfigBasePath = CONFIG_DIR;
  const std::string kOutputBasePath = OUTPUT_DIR;
  const std::string kConfigDirPath = kConfigBasePath + "/impulse_orbit_search/";
  const std::string kAstroParamFile = kConfigBasePath + "/astro_param/astro_param.txt";

  // 天文定数読み込み
  AstroConstants<double> astro_params = loadConstants<double>(kAstroParamFile);
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);
  std::cout << "<>    mu parameter: " << std::setprecision(15) << kMU << std::endl;

  // configファイル検索
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = args.is_continuous;
  std::vector<std::string> config_file_list =
      DiscoverConfigFilesToml(kConfigDirPath, "impulse_orbit_search_config", discovery_opts);

  if (config_file_list.empty()) {
    std::cerr << "<> No config files found in " << kConfigDirPath << std::endl;
    return -1;
  }

  std::cout << "<> Loaded config file list:" << std::endl;
  for (const auto& filepath : config_file_list) {
    std::cout << "<>    - " << filepath << std::endl;
  }

  // OpenMP設定
  int core_max = omp_get_max_threads();
  int omp_threads = std::max(1, core_max - 1);
  if (!skip_wait) {
    std::cout << "<>" << std::endl;
    std::cout << "<>  [OpenMP Configuration]" << std::endl;
    std::cout << "<>    Available threads: " << core_max << std::endl;
    std::cout << "<>    Enter number of threads to use: ";
    std::cin >> omp_threads;
    if (omp_threads > core_max) omp_threads = core_max;
    if (omp_threads < 1) omp_threads = 1;
  }
  omp_set_num_threads(omp_threads);
  std::cout << "<>    Using " << omp_threads << " threads" << std::endl;

  if (!skip_wait) {
    WaitForEnter("<>    Press Enter to start...");
  }

  // 出力ディレクトリ作成
  OutputDirResult output_result =
      CreateSessionOutputDir(kOutputBasePath, "impulse_orbit_search", output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string session_output_dir = output_result.session_dir;
  std::cout << "<>    Session output directory: " << session_output_dir << std::endl;

  auto start_total = std::chrono::system_clock::now();

  // 各configファイルを処理
  for (const auto& config_filepath : config_file_list) {
    std::cout << "<>----------------------------------------------------------------" << std::endl;
    std::cout << "<>    Processing: " << config_filepath << std::endl;

    ImpulseOrbitSearchConfig config;
    if (!LoadConfig(config_filepath, &config)) {
      continue;
    }

    // 設定表示
    std::cout << "<>    Target Jacobi: " << config.target_jacobi << std::endl;
    std::cout << "<>    Orbit data dir: " << config.orbit_data_dir << std::endl;
    std::cout << "<>    File pattern: " << config.orbit_file_pattern << std::endl;
    std::cout << "<>    Adaptive sampling: " << (config.enable_adaptive_sampling ? "ON" : "OFF")
              << std::endl;
    if (config.enable_adaptive_sampling) {
      std::cout << "<>      Coarse grid: " << config.theta_divisions_coarse << " x "
                << config.phi_divisions_coarse << std::endl;
      std::cout << "<>      Refine high SALI: " << (config.refine_high_sali ? "ON" : "OFF")
                << std::endl;
      if (config.refine_high_sali) {
        std::cout << "<>      SALI threshold: " << config.sali_refine_threshold << std::endl;
      }
    }

    // 軌道ファイルを発見
    std::vector<std::string> orbit_files =
        DiscoverOrbitFiles(config.orbit_data_dir, config.orbit_file_pattern);
    if (orbit_files.empty()) {
      std::cerr << "<> !err! No orbit files found matching pattern" << std::endl;
      continue;
    }
    std::cout << "<>    Found " << orbit_files.size() << " orbit files" << std::endl;

    // configファイルごとのサブフォルダを作成
    std::string config_basename = fs::path(config_filepath).stem().string();
    std::string config_output_dir = session_output_dir + "/" + config_basename;
    fs::create_directories(config_output_dir);
    std::cout << "<>    Config output dir: " << config_output_dir << std::endl;

    // configファイルをコピー
    std::string config_copy =
        config_output_dir + "/" + fs::path(config_filepath).filename().string();
    try {
      fs::copy_file(config_filepath, config_copy, fs::copy_options::overwrite_existing);
    } catch (...) {
    }

    // 各軌道ファイルを処理
    for (const auto& orbit_file : orbit_files) {
      std::string orbit_basename = fs::path(orbit_file).stem().string();
      std::cout << "<>    Processing orbit file: " << orbit_basename << std::endl;

      // 軌道データ読み込み
      std::vector<OrbitDataRow> orbit_data;
      if (!LoadOrbitData(orbit_file, config.phase_step, &orbit_data)) {
        std::cerr << "<> !err! Failed to load orbit data" << std::endl;
        continue;
      }
      std::cout << "<>    Loaded " << orbit_data.size() << " orbit points" << std::endl;

      // 出力ファイル準備
      std::string output_csv = config_output_dir + "/" + orbit_basename + "_results.csv";
      std::ofstream ofs(output_csv);
      if (!ofs) {
        std::cerr << "<> !err! Cannot create output file: " << output_csv << std::endl;
        continue;
      }

      ofs << std::setprecision(15) << std::fixed;
      ofs << "# Impulse Orbit Search Results\n";
      ofs << "# Source: " << orbit_file << "\n";
      ofs << "# Target Jacobi: " << config.target_jacobi << "\n";
      ofs << "orbit_idx,time_j2000,x,y,z,vx,vy,vz,v_theta,v_phi,jacobi,"
          << "dv_x,dv_y,dv_z,dv_mag,sali,escape,collision,residence_time\n";

      int total_processed = 0;
      int stable_count = 0;

      // Phase 1: 粗い探索用のグリッドを準備
      std::vector<std::pair<double, double>> coarse_grid;  // (theta, phi) in radians
      double theta_step_coarse = 2.0 * M_PI / config.theta_divisions_coarse;
      double phi_step_coarse = M_PI / config.phi_divisions_coarse;
      for (int i_theta = 0; i_theta < config.theta_divisions_coarse; ++i_theta) {
        double theta = i_theta * theta_step_coarse;
        for (int i_phi = 0; i_phi <= config.phi_divisions_coarse; ++i_phi) {
          double phi = -M_PI / 2.0 + i_phi * phi_step_coarse;
          coarse_grid.emplace_back(theta, phi);
        }
      }

      std::cout << "<>    Coarse grid points: " << coarse_grid.size() << std::endl;
      std::cout << "<>    Total initial calculations: " << orbit_data.size() * coarse_grid.size()
                << std::endl;

// 各軌道点を処理
#pragma omp parallel for schedule(dynamic) reduction(+ : total_processed, stable_count)
      for (int orbit_idx = 0; orbit_idx < static_cast<int>(orbit_data.size()); ++orbit_idx) {
        const OrbitDataRow& row = orbit_data[orbit_idx];

        // J2000 → CR3BP変換
        State<double> asteroid_rot =
            ConvertInertial2RotatingV2(row.asteroid, row.earth, astro_params);
        State<double> original_velocity = {
            0, 0, 0, asteroid_rot.vx, asteroid_rot.vy, asteroid_rot.vz};

        // 目標ヤコビから速度の大きさを計算
        State3d<double> pos = {asteroid_rot.x, asteroid_rot.y, asteroid_rot.z};
        double v_abs = calc_v_abs(pos, config.target_jacobi, kMU);

#ifdef DEBUG_MODE
        // [DEBUG] 最初の数点のみデバッグ出力
        if (orbit_idx < 3) {
#pragma omp critical
          {
            std::cout << std::setprecision(12) << std::fixed;
            std::cout << "\n[DEBUG] === Orbit Point " << orbit_idx << " ===" << std::endl;
            std::cout << "[DEBUG] J2000 asteroid pos: (" << row.asteroid.x << ", " << row.asteroid.y
                      << ", " << row.asteroid.z << ")" << std::endl;
            std::cout << "[DEBUG] J2000 asteroid vel: (" << row.asteroid.vx << ", "
                      << row.asteroid.vy << ", " << row.asteroid.vz << ")" << std::endl;
            std::cout << "[DEBUG] CR3BP rotated pos: (" << asteroid_rot.x << ", " << asteroid_rot.y
                      << ", " << asteroid_rot.z << ")" << std::endl;
            std::cout << "[DEBUG] CR3BP rotated vel: (" << asteroid_rot.vx << ", "
                      << asteroid_rot.vy << ", " << asteroid_rot.vz << ")" << std::endl;
            double original_jacobi = calc_jacobi_integral(asteroid_rot, kMU);
            std::cout << "[DEBUG] Original Jacobi (before impulse): " << original_jacobi
                      << std::endl;
            std::cout << "[DEBUG] Target Jacobi: " << config.target_jacobi << std::endl;
            std::cout << "[DEBUG] Calculated v_abs from target Jacobi: " << v_abs << std::endl;
            double original_v_mag =
                std::sqrt(asteroid_rot.vx * asteroid_rot.vx + asteroid_rot.vy * asteroid_rot.vy +
                          asteroid_rot.vz * asteroid_rot.vz);
            std::cout << "[DEBUG] Original velocity magnitude: " << original_v_mag << std::endl;
          }
        }
#endif

        if (v_abs <= 0) {
          // 禁止領域（速度が虚数）
#ifdef DEBUG_MODE
          if (orbit_idx < 3) {
#pragma omp critical
            {
              std::cout << "[DEBUG] SKIPPED: v_abs <= 0 (forbidden region)" << std::endl;
            }
          }
#endif
          continue;
        }

        // ローカル結果バッファ
        std::vector<SearchResult> local_results;
        std::vector<std::pair<double, double>> refine_candidates;  // 細分化候補

        // Phase 1: 粗い探索
        for (const auto& [theta, phi] : coarse_grid) {
          State3d<double> velocity = SphericalToCartesian(theta, phi, v_abs);
          State<double> initial_state = {asteroid_rot.x, asteroid_rot.y, asteroid_rot.z,
                                         velocity.x,     velocity.y,     velocity.z};

          // Δvフィルタ
          if (config.enable_dv_filter) {
            double dv_x = velocity.x - original_velocity.vx;
            double dv_y = velocity.y - original_velocity.vy;
            double dv_z = velocity.z - original_velocity.vz;
            double dv_mag = std::sqrt(dv_x * dv_x + dv_y * dv_y + dv_z * dv_z);
            if (dv_mag > config.dv_max_threshold) {
              continue;
            }
          }

          SearchResult result =
              RunSaliCalculation(initial_state, original_velocity, orbit_idx, row.time_j2000,
                                 theta * 180.0 / M_PI, phi * 180.0 / M_PI, kMU, config);

#ifdef DEBUG_MODE
          // [DEBUG] 最初の軌道点の最初のいくつかの方向についてデバッグ出力
          if (orbit_idx == 0 && total_processed < 5) {
#pragma omp critical
            {
              std::cout << std::setprecision(12) << std::fixed;
              std::cout << "\n[DEBUG] --- Velocity direction (theta=" << theta * 180.0 / M_PI
                        << "°, phi=" << phi * 180.0 / M_PI << "°) ---" << std::endl;
              std::cout << "[DEBUG] Velocity vector: (" << velocity.x << ", " << velocity.y << ", "
                        << velocity.z << ")" << std::endl;
              double v_mag_check = std::sqrt(velocity.x * velocity.x + velocity.y * velocity.y +
                                             velocity.z * velocity.z);
              std::cout << "[DEBUG] Velocity magnitude (should be " << v_abs << "): " << v_mag_check
                        << std::endl;
              std::cout << "[DEBUG] Computed Jacobi: " << result.jacobi << std::endl;
              std::cout << "[DEBUG] Delta-V magnitude: " << result.dv_mag << std::endl;
              std::cout << "[DEBUG] SALI: " << result.sali << std::endl;
              std::cout << "[DEBUG] Escape: " << result.escape_flag
                        << ", Collision: " << result.collision_flag << std::endl;
            }
          }
#endif

          // 細分化条件のチェック
          bool should_refine = true;
          if (config.refine_no_escape && result.escape_flag != 0) {
            should_refine = false;
          }
          if (config.refine_no_collision && result.collision_flag != 0) {
            should_refine = false;
          }
          if (config.refine_high_sali && result.sali < config.sali_refine_threshold) {
            should_refine = false;
          }

          if (config.enable_adaptive_sampling && should_refine) {
            refine_candidates.emplace_back(theta, phi);
          }

          local_results.push_back(result);
          total_processed++;

          if (result.escape_flag == 0 && result.collision_flag == 0 &&
              result.sali >= config.sali_refine_threshold) {
            stable_count++;
          }
        }

        // Phase 2: 細分化探索
        if (config.enable_adaptive_sampling) {
          double theta_step_fine = theta_step_coarse / config.theta_divisions_fine;
          double phi_step_fine = phi_step_coarse / config.phi_divisions_fine;

          for (const auto& [center_theta, center_phi] : refine_candidates) {
            for (int di = -config.theta_divisions_fine / 2; di <= config.theta_divisions_fine / 2;
                 ++di) {
              for (int dj = -config.phi_divisions_fine / 2; dj <= config.phi_divisions_fine / 2;
                   ++dj) {
                if (di == 0 && dj == 0) continue;  // 既に計算済み

                double theta = center_theta + di * theta_step_fine;
                double phi = center_phi + dj * phi_step_fine;

                // 範囲チェック
                if (phi < -M_PI / 2.0 || phi > M_PI / 2.0) continue;

                State3d<double> velocity = SphericalToCartesian(theta, phi, v_abs);
                State<double> initial_state = {asteroid_rot.x, asteroid_rot.y, asteroid_rot.z,
                                               velocity.x,     velocity.y,     velocity.z};

                if (config.enable_dv_filter) {
                  double dv_x = velocity.x - original_velocity.vx;
                  double dv_y = velocity.y - original_velocity.vy;
                  double dv_z = velocity.z - original_velocity.vz;
                  double dv_mag = std::sqrt(dv_x * dv_x + dv_y * dv_y + dv_z * dv_z);
                  if (dv_mag > config.dv_max_threshold) {
                    continue;
                  }
                }

                SearchResult result =
                    RunSaliCalculation(initial_state, original_velocity, orbit_idx, row.time_j2000,
                                       theta * 180.0 / M_PI, phi * 180.0 / M_PI, kMU, config);

                local_results.push_back(result);
                total_processed++;

                if (result.escape_flag == 0 && result.collision_flag == 0 &&
                    result.sali >= config.sali_refine_threshold) {
                  stable_count++;
                }
              }
            }
          }
        }

// 結果をファイルに書き込み
#pragma omp critical
        {
          for (const auto& result : local_results) {
            if (config.output_stable_only) {
              if (result.escape_flag != 0 || result.collision_flag != 0) {
                continue;
              }
            }
            ofs << result.orbit_idx << "," << result.time_j2000 << "," << result.x << ","
                << result.y << "," << result.z << "," << result.vx << "," << result.vy << ","
                << result.vz << "," << result.v_theta << "," << result.v_phi << "," << result.jacobi
                << "," << result.dv_x << "," << result.dv_y << "," << result.dv_z << ","
                << result.dv_mag << "," << result.sali << "," << result.escape_flag << ","
                << result.collision_flag << "," << result.residence_time << "\n";
          }
        }

        // 進捗表示
        if (orbit_idx % 10 == 0) {
#pragma omp critical
          {
            std::cout << "\r<>    Progress: " << orbit_idx << " / " << orbit_data.size()
                      << " orbit points processed" << std::flush;
          }
        }
      }

      std::cout << std::endl;
      std::cout << "<>    Total calculations: " << total_processed << std::endl;
      std::cout << "<>    Stable orbits found: " << stable_count << std::endl;
      std::cout << "<>    Output: " << output_csv << std::endl;

      ofs.close();
    }  // end orbit_files loop
  }  // end config_file_list loop

  auto end_total = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_total - start_total);
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    Total elapsed time: " << elapsed.count() << " seconds" << std::endl;
  std::cout << "<>    Done." << std::endl;

  return 0;
}
