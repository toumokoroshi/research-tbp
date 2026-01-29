/**
 * @file main.cpp
 * @brief 安定多様体計算アプリケーション
 *
 * 既知のUPO（L1/L2 Lyapunov軌道）に対して安定多様体を計算し、
 * 各点の状態ベクトル（位置・速度）を出力する。
 *
 * @date 2026-01-05
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

#include "periodic_orbit.hpp"
#include "rtbp.hpp"
#include "utils.hpp"

// BOOST_NO_EXCEPTIONS が定義されている場合、boost::throw_exception を提供する必要がある
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>

#ifdef BOOST_NO_EXCEPTIONS
namespace boost {
BOOST_NORETURN void throw_exception(std::exception const& e) {
  std::cerr << "Boost exception: " << e.what() << std::endl;
  std::abort();
}

BOOST_NORETURN void throw_exception(std::exception const& e,
                                     boost::source_location const& loc) {
  std::cerr << "Boost exception at " << loc.file_name() << ":" << loc.line()
            << " - " << e.what() << std::endl;
  std::abort();
}
}  // namespace boost
#endif


namespace fs = std::filesystem;

using namespace my_type;
using namespace periodic_orbit;
using namespace utils;

// ---------------------------------------------------------------------------
// 設定構造体
// ---------------------------------------------------------------------------

struct StableManifoldConfig {
  // 入力モード
  std::string input_mode = "direct";  // "direct" or "file"
  std::string orbit_data_file;
  std::string orbit_type = "lyapunov";  // "lyapunov" or "halo"

  // 直接指定の場合の初期条件
  double initial_x = 0.9899;
  double initial_y = 0.0;
  double initial_z = 0.0;
  double initial_vx = 0.0;
  double initial_vy = 0.05;
  double initial_vz = 0.0;
  double orbit_period_guess = 3.0;

  // Newton-Raphson設定
  int newton_max_iterations = 100;
  double newton_tolerance = 1e-12;

  // 多様体計算設定
  ManifoldConfig<double> manifold_config;

  // 積分器設定
  std::string integrator_type = "rk4";
  double abs_tolerance = 1e-12;
  double rel_tolerance = 1e-12;
  double timestep = 0.0001;
  double max_integration_time = 50.0;

  // システムパラメータ
  double mu = 3.0404233870218e-06;

  // ポアンカレ断面
  std::string section_var = "y";
  double section_value = 0.0;

  // 出力設定
  std::string output_dir = "data/stable_manifold";
  std::string tag = "L1_lyapunov";
  bool save_trajectory = true;
  bool output_velocity_map = true;
};

// ---------------------------------------------------------------------------
// 設定ファイル読み込み
// ---------------------------------------------------------------------------

bool LoadConfig(const std::string& filepath, StableManifoldConfig* config) {
  try {
    auto data = toml::parse_file(filepath);

    // 入力設定
    if (auto val = data["orbit_input"]["input_mode"].value<std::string>()) {
      config->input_mode = *val;
    }
    if (auto val = data["orbit_input"]["orbit_data_file"].value<std::string>()) {
      config->orbit_data_file = *val;
    }
    if (auto val = data["orbit_input"]["orbit_type"].value<std::string>()) {
      config->orbit_type = *val;
    }
    if (auto val = data["orbit_input"]["initial_x"].value<double>()) {
      config->initial_x = *val;
    }
    if (auto val = data["orbit_input"]["initial_y"].value<double>()) {
      config->initial_y = *val;
    }
    if (auto val = data["orbit_input"]["initial_z"].value<double>()) {
      config->initial_z = *val;
    }
    if (auto val = data["orbit_input"]["initial_vx"].value<double>()) {
      config->initial_vx = *val;
    }
    if (auto val = data["orbit_input"]["initial_vy"].value<double>()) {
      config->initial_vy = *val;
    }
    if (auto val = data["orbit_input"]["initial_vz"].value<double>()) {
      config->initial_vz = *val;
    }
    if (auto val = data["orbit_input"]["orbit_period"].value<double>()) {
      config->orbit_period_guess = *val;
    }

    // Newton-Raphson設定
    if (auto val = data["newton_raphson"]["max_iterations"].value<int>()) {
      config->newton_max_iterations = *val;
    }
    if (auto val = data["newton_raphson"]["tolerance"].value<double>()) {
      config->newton_tolerance = *val;
    }

    // 多様体計算設定
    if (auto val = data["manifold"]["epsilon"].value<double>()) {
      config->manifold_config.epsilon = *val;
    }
    if (auto val = data["manifold"]["backward_time"].value<double>()) {
      config->manifold_config.backward_time = *val;
    }
    if (auto val = data["manifold"]["forward_time"].value<double>()) {
      config->manifold_config.forward_time = *val;
    }
    if (auto val = data["manifold"]["num_initial_points"].value<int>()) {
      config->manifold_config.num_initial_points = *val;
    }
    if (auto val = data["manifold"]["compute_stable"].value<bool>()) {
      config->manifold_config.compute_stable = *val;
    }
    if (auto val = data["manifold"]["compute_unstable"].value<bool>()) {
      config->manifold_config.compute_unstable = *val;
    }

    // 積分器設定
    if (auto val = data["integrator"]["type"].value<std::string>()) {
      config->integrator_type = *val;
    }
    if (auto val = data["integrator"]["abs_tolerance"].value<double>()) {
      config->abs_tolerance = *val;
    }
    if (auto val = data["integrator"]["rel_tolerance"].value<double>()) {
      config->rel_tolerance = *val;
    }
    if (auto val = data["integrator"]["timestep"].value<double>()) {
      config->timestep = *val;
    }
    if (auto val = data["integrator"]["max_integration_time"].value<double>()) {
      config->max_integration_time = *val;
    }

    config->manifold_config.integrator_type = config->integrator_type;
    config->manifold_config.abs_tolerance = config->abs_tolerance;
    config->manifold_config.rel_tolerance = config->rel_tolerance;

    // システムパラメータ
    if (auto val = data["system"]["mu"].value<double>()) {
      config->mu = *val;
    }

    // ポアンカレ断面
    if (auto val = data["poincare_section"]["section_var"].value<std::string>()) {
      config->section_var = *val;
    }
    if (auto val = data["poincare_section"]["section_value"].value<double>()) {
      config->section_value = *val;
    }

    // 出力設定
    if (auto val = data["output"]["output_dir"].value<std::string>()) {
      config->output_dir = *val;
    }
    if (auto val = data["output"]["tag"].value<std::string>()) {
      config->tag = *val;
    }
    if (auto val = data["output"]["save_trajectory"].value<bool>()) {
      config->save_trajectory = *val;
    }
    if (auto val = data["output"]["output_velocity_map"].value<bool>()) {
      config->output_velocity_map = *val;
    }

    return true;
  } catch (const toml::parse_error& err) {
    std::cerr << "Config parse error: " << err << std::endl;
    return false;
  } catch (const std::exception& e) {
    std::cerr << "Config load error: " << e.what() << std::endl;
    return false;
  }
}

// ---------------------------------------------------------------------------
// 設定表示
// ---------------------------------------------------------------------------

void PrintConfig(const StableManifoldConfig& config) {
  std::cout << "\n========================================\n";
  std::cout << "  Configuration Summary\n";
  std::cout << "========================================\n\n";

  std::cout << "[Input Mode]\n";
  std::cout << "  Mode: " << config.input_mode << "\n";
  if (config.input_mode == "direct") {
    std::cout << "  Initial state: (" << config.initial_x << ", " << config.initial_y << ", "
              << config.initial_z << ", " << config.initial_vx << ", " << config.initial_vy << ", "
              << config.initial_vz << ")\n";
  } else {
    std::cout << "  Orbit file: " << config.orbit_data_file << "\n";
  }

  std::cout << "\n[Manifold Settings]\n";
  std::cout << "  Epsilon: " << config.manifold_config.epsilon << "\n";
  std::cout << "  Backward time: " << config.manifold_config.backward_time << "\n";
  std::cout << "  Forward time: " << config.manifold_config.forward_time << "\n";
  std::cout << "  Num initial points: " << config.manifold_config.num_initial_points << "\n";
  std::cout << "  Compute stable: " << (config.manifold_config.compute_stable ? "Yes" : "No")
            << "\n";
  std::cout << "  Compute unstable: " << (config.manifold_config.compute_unstable ? "Yes" : "No")
            << "\n";

  std::cout << "\n[System Parameters]\n";
  std::cout << "  mu: " << std::scientific << config.mu << std::fixed << "\n";

  std::cout << "\n[Output Settings]\n";
  std::cout << "  Output directory: " << config.output_dir << "\n";
  std::cout << "  Tag: " << config.tag << "\n";
  std::cout << "========================================\n\n";
}

// ---------------------------------------------------------------------------
// 多様体データ保存
// ---------------------------------------------------------------------------

void SaveManifolds(const std::vector<ManifoldTrajectory<double>>& manifolds,
                   const std::string& output_dir, double mu, double target_cj,
                   const std::string& section_var, double section_value,
                   bool output_velocity_map) {
  struct CjDeviationStats {
    size_t count = 0;
    double mean = 0.0;
    double m2 = 0.0;
    double min = std::numeric_limits<double>::infinity();
    double max = 0.0;
  };

  auto UpdateStats = [](CjDeviationStats* stats, double value) {
    stats->count += 1;
    stats->min = std::min(stats->min, value);
    stats->max = std::max(stats->max, value);
    const double delta = value - stats->mean;
    stats->mean += delta / static_cast<double>(stats->count);
    const double delta2 = value - stats->mean;
    stats->m2 += delta * delta2;
  };

  CjDeviationStats cj_stats;
  // 個別の軌道ファイル
  for (size_t i = 0; i < manifolds.size(); ++i) {
    const auto& manifold = manifolds[i];
    std::string type_str =
        (manifold.type == ManifoldTrajectory<double>::Type::STABLE) ? "stable" : "unstable";

    std::ostringstream filename;
    filename << "manifold_" << type_str << "_" << std::setfill('0') << std::setw(4) << i << ".csv";

    std::ofstream file(fs::path(output_dir) / filename.str());
    file << std::setprecision(16);
    file << "# Manifold Type: " << type_str << "\n";
    file << "# time,x,y,z,vx,vy,vz,jacobi\n";

    for (size_t j = 0; j < manifold.trajectory.size(); ++j) {
      const auto& state = manifold.trajectory[j];
      double jacobi = crtbp::calc_jacobi_integral(state, mu);
      file << manifold.times[j] << "," << state.x << "," << state.y << "," << state.z << ","
           << state.vx << "," << state.vy << "," << state.vz << "," << jacobi << "\n";
    }
    file.close();
  }

  std::cout << "Saved " << manifolds.size() << " manifold trajectory files\n";

  // 速度マップ出力
  if (output_velocity_map) {
    std::ofstream map_file(fs::path(output_dir) / "velocity_map.csv");
    map_file << std::setprecision(16);
    map_file << "# Velocity Map: Position and velocity on stable manifold\n";
    map_file << "# x,y,z,vx,vy,vz,jacobi,manifold_idx,time\n";

    for (size_t i = 0; i < manifolds.size(); ++i) {
      const auto& manifold = manifolds[i];
      if (manifold.type != ManifoldTrajectory<double>::Type::STABLE) {
        continue;  // 速度マップは安定多様体のみ
      }

      for (size_t j = 0; j < manifold.trajectory.size(); ++j) {
        const auto& state = manifold.trajectory[j];
        double jacobi = crtbp::calc_jacobi_integral(state, mu);
        double cj_dev = std::abs(jacobi - target_cj);
        UpdateStats(&cj_stats, cj_dev);
        map_file << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
                 << state.vy << "," << state.vz << "," << jacobi << "," << i << ","
                 << manifold.times[j] << "\n";
      }
    }
    map_file.close();
    std::cout << "Saved velocity_map.csv\n";

    if (cj_stats.count > 0) {
      const double variance =
          cj_stats.m2 / static_cast<double>(cj_stats.count);
      const double stddev = std::sqrt(variance);

      std::ofstream stats_file(fs::path(output_dir) / "velocity_map_cj_stats.txt");
      stats_file << std::setprecision(16);
      stats_file << "Cj Deviation Statistics (|Cj - target|)\n";
      stats_file << "========================================\n";
      stats_file << "Target Cj: " << target_cj << "\n";
      stats_file << "Count: " << cj_stats.count << "\n";
      stats_file << "Min: " << cj_stats.min << "\n";
      stats_file << "Max: " << cj_stats.max << "\n";
      stats_file << "Mean: " << cj_stats.mean << "\n";
      stats_file << "StdDev: " << stddev << "\n";
      stats_file.close();
      std::cout << "Saved velocity_map_cj_stats.txt\n";
    }
  }

  // 断面交点の抽出（速度マップとは独立に出力）
  auto GetComponent = [](const State<double>& s, int index) -> double {
    switch (index) {
      case 0:
        return s.x;
      case 1:
        return s.y;
      case 2:
        return s.z;
      case 3:
        return s.vx;
      case 4:
        return s.vy;
      case 5:
        return s.vz;
      default:
        return 0.0;
    }
  };

  int section_index = -1;
  if (section_var == "x")
    section_index = 0;
  else if (section_var == "y")
    section_index = 1;
  else if (section_var == "z")
    section_index = 2;
  else if (section_var == "vx")
    section_index = 3;
  else if (section_var == "vy")
    section_index = 4;
  else if (section_var == "vz")
    section_index = 5;

  if (section_index >= 0) {
    std::ofstream section_file(fs::path(output_dir) / "manifold_section_hits.csv");
    section_file << std::setprecision(16);
    section_file
        << "# manifold_idx,type,time,x,y,z,vx,vy,vz,jacobi,section_var,section_value,sign\n";

    const double eps = 1e-14;
    for (size_t i = 0; i < manifolds.size(); ++i) {
      const auto& manifold = manifolds[i];
      const std::string type_str =
          (manifold.type == ManifoldTrajectory<double>::Type::STABLE) ? "stable" : "unstable";

      for (size_t j = 1; j < manifold.trajectory.size(); ++j) {
        const auto& prev = manifold.trajectory[j - 1];
        const auto& curr = manifold.trajectory[j];
        const double v_prev = GetComponent(prev, section_index) - section_value;
        const double v_curr = GetComponent(curr, section_index) - section_value;

        if (std::abs(v_prev) < eps && std::abs(v_curr) < eps) {
          continue;
        }

        if (v_prev == 0.0) {
          const double jacobi = crtbp::calc_jacobi_integral(prev, mu);
          const int sign = (v_curr > 0.0) ? 1 : -1;
          section_file << i << "," << type_str << "," << manifold.times[j - 1] << ","
                       << prev.x << "," << prev.y << "," << prev.z << "," << prev.vx << ","
                       << prev.vy << "," << prev.vz << "," << jacobi << "," << section_var << ","
                       << section_value << "," << sign << "\n";
          continue;
        }

        if ((v_prev > 0.0 && v_curr > 0.0) || (v_prev < 0.0 && v_curr < 0.0)) {
          continue;
        }

        const double denom = v_curr - v_prev;
        if (std::abs(denom) < eps) {
          continue;
        }
        const double alpha = -v_prev / denom;
        if (alpha < 0.0 || alpha > 1.0) {
          continue;
        }

        State<double> hit;
        hit.x = prev.x + alpha * (curr.x - prev.x);
        hit.y = prev.y + alpha * (curr.y - prev.y);
        hit.z = prev.z + alpha * (curr.z - prev.z);
        hit.vx = prev.vx + alpha * (curr.vx - prev.vx);
        hit.vy = prev.vy + alpha * (curr.vy - prev.vy);
        hit.vz = prev.vz + alpha * (curr.vz - prev.vz);
        const double hit_time = manifold.times[j - 1] +
                                alpha * (manifold.times[j] - manifold.times[j - 1]);
        const double jacobi = crtbp::calc_jacobi_integral(hit, mu);
        const int sign = (v_curr > v_prev) ? 1 : -1;

        section_file << i << "," << type_str << "," << hit_time << "," << hit.x << "," << hit.y
                     << "," << hit.z << "," << hit.vx << "," << hit.vy << "," << hit.vz << ","
                     << jacobi << "," << section_var << "," << section_value << "," << sign
                     << "\n";
      }
    }
    section_file.close();
    std::cout << "Saved manifold_section_hits.csv\n";
  }
}

// ---------------------------------------------------------------------------
// 周期軌道情報保存
// ---------------------------------------------------------------------------

void SaveOrbitInfo(const PeriodicOrbit<double>& orbit, const std::string& output_dir) {
  std::ofstream file(fs::path(output_dir) / "orbit_info.txt");
  file << std::setprecision(16);

  file << "Periodic Orbit Information\n";
  file << "==========================\n\n";

  file << "Initial State:\n";
  file << "  x  = " << orbit.initial_state.x << "\n";
  file << "  y  = " << orbit.initial_state.y << "\n";
  file << "  z  = " << orbit.initial_state.z << "\n";
  file << "  vx = " << orbit.initial_state.vx << "\n";
  file << "  vy = " << orbit.initial_state.vy << "\n";
  file << "  vz = " << orbit.initial_state.vz << "\n\n";

  file << "Period: " << orbit.period << "\n";
  file << "Jacobi Constant: " << orbit.jacobi_constant << "\n\n";

  if (orbit.stability_computed) {
    file << "Stability Analysis:\n";
    file << "  Stable: " << (orbit.is_stable ? "Yes (SPO)" : "No (UPO)") << "\n";
    file << "  Stability Index: " << orbit.stability_index << "\n\n";

    file << "Eigenvalues:\n";
    for (size_t i = 0; i < orbit.eigenvalues.size(); ++i) {
      file << "  lambda" << i + 1 << " = " << orbit.eigenvalues[i]
           << " (|lambda| = " << std::abs(orbit.eigenvalues[i]) << ")\n";
    }
  }

  file.close();
  std::cout << "Saved orbit_info.txt\n";
}

// ---------------------------------------------------------------------------
// メイン関数
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  std::cout << "========================================\n";
  std::cout << "  Stable Manifold Computation\n";
  std::cout << "  CR3BP Invariant Manifold Analysis\n";
  std::cout << "========================================\n\n";

  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv, "StableManifoldComputation");
  std::string output_tag = args.output_tag;

  // configファイル検索
  const std::string kConfigPath = std::string(CONFIG_DIR) + "/stable_manifold";
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = args.is_continuous;
  std::vector<std::string> config_files =
      DiscoverConfigFilesToml(kConfigPath, "stable_manifold_config", discovery_opts);

  if (config_files.empty()) {
    std::cerr << "No config files found in " << kConfigPath << std::endl;
    return 1;
  }

  std::cout << "Found " << config_files.size() << " config file(s):\n";
  for (const auto& f : config_files) {
    std::cout << "  - " << f << "\n";
  }

  // 最初のconfigファイルを使用
  std::string config_file = config_files[0];

  // 設定読み込み
  StableManifoldConfig config;
  if (!LoadConfig(config_file, &config)) {
    std::cerr << "Failed to load configuration file\n";
    return 1;
  }

  PrintConfig(config);

  // 出力ディレクトリ作成
  OutputDirResult output_result = CreateSessionOutputDir(
      config.output_dir, "stable_manifold", output_tag.empty() ? config.tag : output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string output_dir = output_result.session_dir;

  auto start_time = std::chrono::system_clock::now();

  // ポアンカレ断面のインデックス取得
  int section_index = -1;
  if (config.section_var == "x")
    section_index = 0;
  else if (config.section_var == "y")
    section_index = 1;
  else if (config.section_var == "z")
    section_index = 2;
  else if (config.section_var == "vx")
    section_index = 3;
  else if (config.section_var == "vy")
    section_index = 4;
  else if (config.section_var == "vz")
    section_index = 5;

  if (section_index < 0) {
    std::cerr << "Invalid section variable: " << config.section_var << "\n";
    return 1;
  }

  // 周期軌道取得
  PeriodicOrbit<double> orbit;

  if (config.input_mode == "direct") {
    std::cout << "\n[Step 1] Refining periodic orbit from initial guess...\n";

    State<double> initial_guess = {config.initial_x,  config.initial_y,  config.initial_z,
                                   config.initial_vx, config.initial_vy, config.initial_vz};

    try {
      NewtonConvergenceInfo<double> conv_info;

      if (config.orbit_type == "halo") {
        // Halo軌道は2変数Newton法（z₀固定）を使用
        // (x₀, vy₀) を調整して半周期で vx(T/2) = 0, vz(T/2) = 0 となる条件を満たす
        std::cout << "  Orbit type: Halo (z0 = " << config.initial_z << " fixed)\n";
        orbit = RefinePeriodicOrbitHalo(initial_guess, config.mu,
                                        config.newton_max_iterations, config.newton_tolerance,
                                        config.max_integration_time, config.timestep, &conv_info);
      } else {
        // Lyapunov軌道は対称性を利用した2変数Newton法を使用
        // (x₀, vy₀) を調整して半周期で vx(T/2) = 0 となる条件を満たす
        std::cout << "  Orbit type: Lyapunov (z = 0)\n";
        orbit = RefinePeriodicOrbitSymmetric(initial_guess, config.mu,
                                             config.newton_max_iterations, config.newton_tolerance,
                                             config.max_integration_time, config.timestep, &conv_info);
      }

      std::cout << "  Converged in " << conv_info.iterations << " iterations\n";
      std::cout << "  Final residual: " << std::scientific << conv_info.final_residual << std::fixed
                << "\n";
      std::cout << "  Period: " << orbit.period << "\n";
      std::cout << "  Jacobi constant: " << orbit.jacobi_constant << "\n";

    } catch (const std::exception& e) {
      std::cerr << "Failed to refine periodic orbit: " << e.what() << "\n";
      return 1;
    }
  } else {
    std::cerr << "File input mode not implemented yet\n";
    return 1;
  }

  // 安定性解析
  std::cout << "\n[Step 2] Computing monodromy matrix and analyzing stability...\n";

  try {
    orbit.monodromy_matrix = ComputeMonodromyMatrix(orbit, config.mu, config.timestep);
    AnalyzeStability(&orbit, 1.0);

    std::cout << "  Stability: " << (orbit.is_stable ? "SPO (Stable)" : "UPO (Unstable)") << "\n";
    std::cout << "  Stability index: " << orbit.stability_index << "\n";

    std::cout << "  Eigenvalues:\n";
    for (size_t i = 0; i < orbit.eigenvalues.size(); ++i) {
      std::cout << "    lambda" << i + 1 << " = " << orbit.eigenvalues[i]
                << " (|lambda| = " << std::abs(orbit.eigenvalues[i]) << ")\n";
    }

  } catch (const std::exception& e) {
    std::cerr << "Stability analysis failed: " << e.what() << "\n";
    return 1;
  }

  // 不変多様体計算
  std::cout << "\n[Step 3] Computing invariant manifolds...\n";

  std::vector<ManifoldTrajectory<double>> manifolds;
  try {
    manifolds =
        ComputeInvariantManifolds(orbit, config.mu, config.manifold_config, config.timestep);
    std::cout << "  Computed " << manifolds.size() << " manifold trajectories\n";

  } catch (const std::exception& e) {
    std::cerr << "Manifold computation failed: " << e.what() << "\n";
    return 1;
  }

  // 結果保存
  std::cout << "\n[Step 4] Saving results...\n";

  SaveOrbitInfo(orbit, output_dir);
  SaveManifolds(manifolds, output_dir, config.mu, orbit.jacobi_constant,
                config.section_var, config.section_value, config.output_velocity_map);

  auto end_time = std::chrono::system_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

  std::cout << "\n========================================\n";
  std::cout << "  Computation complete\n";
  std::cout << "  Elapsed time: " << elapsed / 1000.0 << " seconds\n";
  std::cout << "  Output: " << output_dir << "\n";
  std::cout << "========================================\n";

  return 0;
}
