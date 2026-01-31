/**
 * @file main.cpp
 * @brief 任意ヤコビ積分値における軌道・不変多様体計算アプリケーション
 *
 * 軌道族CSVから指定されたヤコビ積分値に最も近い軌道を探索し、
 * Newton法で精密化した後、不変多様体を計算する。
 *
 * @date 2026-01-29
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

struct OrbitManifoldConfig {
  // ターゲット設定
  std::string orbit_type = "halo";       // "lyapunov" or "halo"
  std::string lagrange_point = "L1";     // "L1" or "L2"
  std::string halo_branch = "north";     // "north" or "south" (haloの場合)
  double target_jacobi = 3.0005;

  // 軌道族データベース
  std::string l1_lyapunov_file = "data/l1_lyapunov_continuation/l1_lyapunov_family.csv";
  std::string l2_lyapunov_file = "data/lyapunov_continuation/l2_lyapunov_family.csv";
  std::string l1_halo_north_file = "data/halo_continuation/l1_halo_north.csv";
  std::string l1_halo_south_file = "data/halo_continuation/l1_halo_south.csv";
  std::string l2_halo_north_file = "data/halo_continuation/l2_halo_north.csv";
  std::string l2_halo_south_file = "data/halo_continuation/l2_halo_south.csv";

  // Newton-Raphson設定
  int newton_max_iterations = 100;
  double newton_tolerance = 1e-12;

  // 多様体計算設定
  ManifoldConfig<double> manifold_config;

  // 積分器設定
  double timestep = 0.0001;
  double max_integration_time = 50.0;

  // システムパラメータ
  double mu = 3.0404233870218e-06;

  // ポアンカレ断面
  std::string section_var = "y";
  double section_value = 0.0;

  // 出力設定
  std::string output_dir = "data/orbit_manifold_by_jacobi";
  bool output_velocity_map = true;
};

// ---------------------------------------------------------------------------
// 軌道データ読み込み
// ---------------------------------------------------------------------------

struct OrbitData {
  State<double> initial_state;
  double period;
  double jacobi_constant;
  bool is_stable;
  double stability_index;
};

std::vector<OrbitData> LoadOrbitFamily(const std::string& filepath) {
  std::vector<OrbitData> orbits;
  std::ifstream file(filepath);

  if (!file.is_open()) {
    std::cerr << "Warning: Could not open orbit file: " << filepath << std::endl;
    return orbits;
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    std::string token;
    std::vector<double> values;

    while (std::getline(iss, token, ',')) {
      try {
        values.push_back(std::stod(token));
      } catch (...) {
        break;
      }
    }

    // index,x0,y0,z0,vx0,vy0,vz0,period,jacobi,stable,stability_index
    if (values.size() >= 9) {
      OrbitData orbit;
      orbit.initial_state.x = values[1];
      orbit.initial_state.y = values[2];
      orbit.initial_state.z = values[3];
      orbit.initial_state.vx = values[4];
      orbit.initial_state.vy = values[5];
      orbit.initial_state.vz = values[6];
      orbit.period = values[7];
      orbit.jacobi_constant = values[8];
      orbit.is_stable = (values.size() >= 10) ? (values[9] != 0) : false;
      orbit.stability_index = (values.size() >= 11) ? values[10] : 0.0;
      orbits.push_back(orbit);
    }
  }

  return orbits;
}

// ---------------------------------------------------------------------------
// 設定ファイル読み込み
// ---------------------------------------------------------------------------

bool LoadConfig(const std::string& filepath, OrbitManifoldConfig* config) {
  try {
    auto data = toml::parse_file(filepath);

    // ターゲット設定
    if (auto val = data["target"]["orbit_type"].value<std::string>()) {
      config->orbit_type = *val;
    }
    if (auto val = data["target"]["lagrange_point"].value<std::string>()) {
      config->lagrange_point = *val;
    }
    if (auto val = data["target"]["halo_branch"].value<std::string>()) {
      config->halo_branch = *val;
    }
    if (auto val = data["target"]["target_jacobi"].value<double>()) {
      config->target_jacobi = *val;
    }

    // 軌道族データベース
    if (auto val = data["orbit_database"]["l1_lyapunov_file"].value<std::string>()) {
      config->l1_lyapunov_file = *val;
    }
    if (auto val = data["orbit_database"]["l2_lyapunov_file"].value<std::string>()) {
      config->l2_lyapunov_file = *val;
    }
    if (auto val = data["orbit_database"]["l1_halo_north_file"].value<std::string>()) {
      config->l1_halo_north_file = *val;
    }
    if (auto val = data["orbit_database"]["l1_halo_south_file"].value<std::string>()) {
      config->l1_halo_south_file = *val;
    }
    if (auto val = data["orbit_database"]["l2_halo_north_file"].value<std::string>()) {
      config->l2_halo_north_file = *val;
    }
    if (auto val = data["orbit_database"]["l2_halo_south_file"].value<std::string>()) {
      config->l2_halo_south_file = *val;
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
    if (auto val = data["integrator"]["timestep"].value<double>()) {
      config->timestep = *val;
    }
    if (auto val = data["integrator"]["max_integration_time"].value<double>()) {
      config->max_integration_time = *val;
    }

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

void PrintConfig(const OrbitManifoldConfig& config) {
  std::cout << "\n========================================\n";
  std::cout << "  Configuration Summary\n";
  std::cout << "========================================\n\n";

  std::cout << "[Target]\n";
  std::cout << "  Orbit type: " << config.orbit_type << "\n";
  std::cout << "  Lagrange point: " << config.lagrange_point << "\n";
  if (config.orbit_type == "halo") {
    std::cout << "  Halo branch: " << config.halo_branch << "\n";
  }
  std::cout << "  Target Jacobi: " << std::setprecision(10) << config.target_jacobi << "\n";

  std::cout << "\n[Manifold Settings]\n";
  std::cout << "  Epsilon: " << config.manifold_config.epsilon << "\n";
  std::cout << "  Backward time: " << config.manifold_config.backward_time << "\n";
  std::cout << "  Forward time: " << config.manifold_config.forward_time << "\n";
  std::cout << "  Num initial points: " << config.manifold_config.num_initial_points << "\n";

  std::cout << "\n[System Parameters]\n";
  std::cout << "  mu: " << std::scientific << config.mu << std::fixed << "\n";

  std::cout << "========================================\n\n";
}

// ---------------------------------------------------------------------------
// 多様体データ保存
// ---------------------------------------------------------------------------

void SaveManifolds(const std::vector<ManifoldTrajectory<double>>& manifolds,
                   const std::string& output_dir, double mu, double target_cj,
                   bool output_velocity_map) {
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
    map_file << "# Velocity Map: Position and velocity on manifold\n";
    map_file << "# x,y,z,vx,vy,vz,jacobi,manifold_idx,time,type\n";

    for (size_t i = 0; i < manifolds.size(); ++i) {
      const auto& manifold = manifolds[i];
      std::string type_str =
          (manifold.type == ManifoldTrajectory<double>::Type::STABLE) ? "stable" : "unstable";

      for (size_t j = 0; j < manifold.trajectory.size(); ++j) {
        const auto& state = manifold.trajectory[j];
        double jacobi = crtbp::calc_jacobi_integral(state, mu);
        map_file << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
                 << state.vy << "," << state.vz << "," << jacobi << "," << i << ","
                 << manifold.times[j] << "," << type_str << "\n";
      }
    }
    map_file.close();
    std::cout << "Saved velocity_map.csv\n";
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
  std::cout << "  Orbit & Manifold by Jacobi Integral\n";
  std::cout << "  CR3BP Invariant Manifold Analysis\n";
  std::cout << "========================================\n\n";

  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv, "OrbitManifoldByJacobi");

  // configファイル検索
  const std::string kConfigPath = std::string(CONFIG_DIR) + "/orbit_manifold_jacobi";
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = args.is_continuous;
  std::vector<std::string> config_files =
      DiscoverConfigFilesToml(kConfigPath, "config", discovery_opts);

  if (config_files.empty()) {
    std::cerr << "No config files found in " << kConfigPath << std::endl;
    std::cerr << "Using default configuration.\n";
  }

  // 設定読み込み
  OrbitManifoldConfig config;
  if (!config_files.empty()) {
    if (!LoadConfig(config_files[0], &config)) {
      std::cerr << "Failed to load configuration file\n";
      return 1;
    }
  }

  PrintConfig(config);

  // 出力ディレクトリ作成
  std::string tag = config.lagrange_point + "_" + config.orbit_type;
  if (config.orbit_type == "halo") {
    tag += "_" + config.halo_branch;
  }
  OutputDirResult output_result = CreateSessionOutputDir(
      config.output_dir, "orbit_manifold", args.output_tag.empty() ? tag : args.output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string output_dir = output_result.session_dir;

  auto start_time = std::chrono::system_clock::now();

  // ---------------------------------------------------------------------------
  // Step 1: 軌道族CSVから目的のC_Jに最も近い軌道を探索
  // ---------------------------------------------------------------------------
  std::cout << "\n[Step 1] Searching for orbit with target Jacobi = "
            << std::setprecision(10) << config.target_jacobi << "...\n";

  // 適切なCSVファイルを選択
  std::string orbit_file;
  if (config.orbit_type == "lyapunov") {
    if (config.lagrange_point == "L1") {
      orbit_file = config.l1_lyapunov_file;
    } else {
      orbit_file = config.l2_lyapunov_file;
    }
  } else {  // halo
    if (config.lagrange_point == "L1") {
      orbit_file = (config.halo_branch == "north") ? config.l1_halo_north_file
                                                    : config.l1_halo_south_file;
    } else {
      orbit_file = (config.halo_branch == "north") ? config.l2_halo_north_file
                                                    : config.l2_halo_south_file;
    }
  }

  std::cout << "  Loading orbit family from: " << orbit_file << "\n";
  std::vector<OrbitData> orbit_family = LoadOrbitFamily(orbit_file);

  if (orbit_family.empty()) {
    std::cerr << "  ERROR: No orbits found in " << orbit_file << "\n";
    return 1;
  }

  std::cout << "  Loaded " << orbit_family.size() << " orbits\n";

  // 最も近いC_Jの軌道を探索
  double min_diff = std::numeric_limits<double>::infinity();
  size_t best_idx = 0;
  for (size_t i = 0; i < orbit_family.size(); ++i) {
    double diff = std::abs(orbit_family[i].jacobi_constant - config.target_jacobi);
    if (diff < min_diff) {
      min_diff = diff;
      best_idx = i;
    }
  }

  OrbitData& closest_orbit = orbit_family[best_idx];
  std::cout << "  Found closest orbit at index " << best_idx << "\n";
  std::cout << "  Jacobi constant: " << closest_orbit.jacobi_constant
            << " (diff = " << min_diff << ")\n";
  std::cout << "  Initial state: x=" << closest_orbit.initial_state.x
            << ", vy=" << closest_orbit.initial_state.vy
            << ", z=" << closest_orbit.initial_state.z << "\n";

  // ---------------------------------------------------------------------------
  // Step 2: Newton法で精密化
  // ---------------------------------------------------------------------------
  std::cout << "\n[Step 2] Refining periodic orbit with Newton-Raphson...\n";

  PeriodicOrbit<double> orbit;
  try {
    NewtonConvergenceInfo<double> conv_info;

    if (config.orbit_type == "halo") {
      orbit = RefinePeriodicOrbitHalo(closest_orbit.initial_state, config.mu,
                                      config.newton_max_iterations, config.newton_tolerance,
                                      config.max_integration_time, config.timestep, &conv_info);
    } else {
      orbit = RefinePeriodicOrbitSymmetric(closest_orbit.initial_state, config.mu,
                                           config.newton_max_iterations, config.newton_tolerance,
                                           config.max_integration_time, config.timestep, &conv_info);
    }

    std::cout << "  Converged in " << conv_info.iterations << " iterations\n";
    std::cout << "  Final residual: " << std::scientific << conv_info.final_residual << std::fixed
              << "\n";
    std::cout << "  Period: " << orbit.period << "\n";
    std::cout << "  Jacobi constant: " << orbit.jacobi_constant << "\n";

  } catch (const std::exception& e) {
    std::cerr << "  Newton refinement failed: " << e.what() << "\n";
    std::cerr << "  Using original orbit data without refinement.\n";
    
    // 精密化に失敗した場合は元データをそのまま使用
    orbit.initial_state = closest_orbit.initial_state;
    orbit.period = closest_orbit.period;
    orbit.jacobi_constant = closest_orbit.jacobi_constant;
    orbit.is_stable = closest_orbit.is_stable;
    orbit.stability_index = closest_orbit.stability_index;
  }

  // ---------------------------------------------------------------------------
  // Step 3: モノドロミー行列・安定性解析
  // ---------------------------------------------------------------------------
  std::cout << "\n[Step 3] Computing monodromy matrix and analyzing stability...\n";

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
    std::cerr << "  Stability analysis failed: " << e.what() << "\n";
    return 1;
  }

  // ---------------------------------------------------------------------------
  // Step 4: 不変多様体計算
  // ---------------------------------------------------------------------------
  std::cout << "\n[Step 4] Computing invariant manifolds...\n";

  std::vector<ManifoldTrajectory<double>> manifolds;
  try {
    manifolds =
        ComputeInvariantManifolds(orbit, config.mu, config.manifold_config, config.timestep);
    std::cout << "  Computed " << manifolds.size() << " manifold trajectories\n";

  } catch (const std::exception& e) {
    std::cerr << "  Manifold computation failed: " << e.what() << "\n";
    return 1;
  }

  // ---------------------------------------------------------------------------
  // Step 5: 結果保存
  // ---------------------------------------------------------------------------
  std::cout << "\n[Step 5] Saving results...\n";

  SaveOrbitInfo(orbit, output_dir);
  SaveManifolds(manifolds, output_dir, config.mu, orbit.jacobi_constant,
                config.output_velocity_map);

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
