/**
 * @file main.cpp
 * @brief 周期軌道継続法アプリケーション
 * @details Lyapunov軌道を起点とした継続法によるUPO探索
 *          研究計画書「手順1. 周期軌道族の系統的同定」に対応
 * @date 2026-01-01
 */

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "continuation.hpp"
#include "lyapunov_initial.hpp"
#include "periodic_orbit.hpp"
#include "rtbp.hpp"
#include "utils.hpp"

// Use toml++
#include <toml++/toml.hpp>

using namespace my_type;
using namespace continuation;
using namespace lyapunov;
using namespace utils;

// ---------------------------------------------------------------------------
// 設定構造体
// ---------------------------------------------------------------------------

struct AppConfig {
  // 継続法設定
  ContinuationConfig<double> continuation;

  // 開始点設定
  std::vector<std::string> starting_points;  // "L1", "L2", "L3"

  // システムパラメータ
  double mu = 3.0404233870218e-06;  // Sun-Earth系

  // 出力設定
  std::string output_dir = "data/periodic_orbit_continuation";
  std::string output_tag = "";
  bool generate_gnuplot_script = true;
  bool save_all_orbits = true;
  bool save_bifurcation_info = true;
};

// ---------------------------------------------------------------------------
// 関数宣言
// ---------------------------------------------------------------------------

bool LoadConfig(const std::string& filepath, AppConfig* config);
void PrintConfig(const AppConfig& config);
LagrangePoint ParseLagrangePoint(const std::string& str);

// ---------------------------------------------------------------------------
// 設定ファイル読み込み
// ---------------------------------------------------------------------------

bool LoadConfig(const std::string& filepath, AppConfig* config) {
  try {
    auto data = toml::parse_file(filepath);

    // 継続法設定
    if (auto method = data["continuation"]["method"].value<std::string>()) {
      if (*method == "pseudo_arclength") {
        config->continuation.method = ContinuationMethod::PSEUDO_ARCLENGTH;
      } else {
        config->continuation.method = ContinuationMethod::NATURAL_PARAMETER;
      }
    }
    if (auto val = data["continuation"]["initial_step_size"].value<double>()) {
      config->continuation.initial_step_size = *val;
    }
    if (auto val = data["continuation"]["min_step_size"].value<double>()) {
      config->continuation.min_step_size = *val;
    }
    if (auto val = data["continuation"]["max_step_size"].value<double>()) {
      config->continuation.max_step_size = *val;
    }
    if (auto val = data["continuation"]["max_steps"].value<int>()) {
      config->continuation.max_steps = *val;
    }
    if (auto val = data["continuation"]["jacobi_min"].value<double>()) {
      config->continuation.jacobi_min = *val;
    }
    if (auto val = data["continuation"]["jacobi_max"].value<double>()) {
      config->continuation.jacobi_max = *val;
    }
    if (auto val = data["continuation"]["initial_amplitude"].value<double>()) {
      config->continuation.initial_amplitude = *val;
    }

    // 開始点
    if (auto arr = data["continuation"]["starting_points"].as_array()) {
      config->starting_points.clear();
      for (const auto& item : *arr) {
        if (auto str = item.value<std::string>()) {
          config->starting_points.push_back(*str);
        }
      }
    }

    // Newton設定
    if (auto val = data["newton"]["max_iterations"].value<int>()) {
      config->continuation.newton_max_iterations = *val;
    }
    if (auto val = data["newton"]["tolerance"].value<double>()) {
      config->continuation.newton_tolerance = *val;
    }

    // 分岐検出設定
    if (auto val = data["bifurcation"]["detect_pitchfork"].value<bool>()) {
      config->continuation.detect_bifurcations = *val;
    }
    if (auto val = data["bifurcation"]["follow_branches"].value<bool>()) {
      config->continuation.follow_branches = *val;
    }
    if (auto val = data["bifurcation"]["eigenvalue_threshold"].value<double>()) {
      config->continuation.eigenvalue_threshold = *val;
    }

    // 並列化設定
    if (auto val = data["parallel"]["num_threads"].value<int>()) {
      config->continuation.num_threads = *val;
    }

    // システムパラメータ
    if (auto val = data["system"]["mu"].value<double>()) {
      config->mu = *val;
    }

    // 出力設定
    if (auto val = data["output"]["output_dir"].value<std::string>()) {
      config->output_dir = *val;
    }
    if (auto val = data["output"]["tag"].value<std::string>()) {
      config->output_tag = *val;
    }
    if (auto val = data["output"]["generate_gnuplot_script"].value<bool>()) {
      config->generate_gnuplot_script = *val;
    }
    if (auto val = data["output"]["save_all_orbits"].value<bool>()) {
      config->save_all_orbits = *val;
    }
    if (auto val = data["output"]["save_bifurcation_info"].value<bool>()) {
      config->save_bifurcation_info = *val;
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

void PrintConfig(const AppConfig& config) {
  std::cout << "\n========================================\n";
  std::cout << "  Configuration Summary\n";
  std::cout << "========================================\n\n";

  std::cout << "[Continuation Method]\n";
  std::cout << "  Method: "
            << (config.continuation.method == ContinuationMethod::PSEUDO_ARCLENGTH
                    ? "Pseudo-arclength"
                    : "Natural Parameter")
            << "\n";
  std::cout << "  Initial step size: " << config.continuation.initial_step_size << "\n";
  std::cout << "  Max steps: " << config.continuation.max_steps << "\n";
  std::cout << "  Jacobi range: [" << config.continuation.jacobi_min << ", "
            << config.continuation.jacobi_max << "]\n";
  std::cout << "  Initial amplitude: " << config.continuation.initial_amplitude << "\n\n";

  std::cout << "[Starting Points]\n";
  for (const auto& pt : config.starting_points) {
    std::cout << "  - " << pt << "\n";
  }
  std::cout << "\n";

  std::cout << "[Newton-Raphson]\n";
  std::cout << "  Max iterations: " << config.continuation.newton_max_iterations << "\n";
  std::cout << "  Tolerance: " << config.continuation.newton_tolerance << "\n\n";

  std::cout << "[Bifurcation Detection]\n";
  std::cout << "  Enabled: " << (config.continuation.detect_bifurcations ? "yes" : "no") << "\n";
  std::cout << "  Follow branches: " << (config.continuation.follow_branches ? "yes" : "no")
            << "\n\n";

  std::cout << "[System]\n";
  std::cout << "  Mass ratio mu: " << std::scientific << config.mu << "\n\n";

  std::cout << "[Output]\n";
  std::cout << "  Directory: " << config.output_dir << "\n";
  std::cout << "========================================\n\n";
}

LagrangePoint ParseLagrangePoint(const std::string& str) {
  if (str == "L1") return LagrangePoint::L1;
  if (str == "L2") return LagrangePoint::L2;
  if (str == "L3") return LagrangePoint::L3;
  throw std::invalid_argument("Unknown Lagrange point: " + str);
}

// ---------------------------------------------------------------------------
// メイン関数
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  auto start_time = std::chrono::high_resolution_clock::now();

  std::cout << "========================================\n";
  std::cout << "  Periodic Orbit Continuation\n";
  std::cout << "  Lyapunov Orbit-Based UPO Search\n";
  std::cout << "========================================\n\n";

  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv, "PeriodicOrbitContinuation");
  std::string output_tag = args.output_tag;
  bool is_continuous = args.is_continuous;

  // 設定ファイル検索
  const std::string kConfigPath = std::string(CONFIG_DIR) + "/periodic_orbit_continuation";
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = is_continuous;
  std::vector<std::string> config_files =
      DiscoverConfigFilesToml(kConfigPath, "continuation_config", discovery_opts);

  if (config_files.empty()) {
    std::cerr << "No config files found in " << kConfigPath << std::endl;
    std::cerr << "Creating default config file..." << std::endl;

    // デフォルト設定ファイルを作成
    namespace fs = std::filesystem;
    fs::create_directories(kConfigPath);
    std::string default_config = kConfigPath + "/continuation_config_1.toml";
    std::ofstream ofs(default_config);
    ofs << R"([continuation]
method = "pseudo_arclength"
starting_points = ["L1", "L2"]
initial_amplitude = 0.001
initial_step_size = 0.001
min_step_size = 1e-6
max_step_size = 0.05
max_steps = 500
jacobi_min = 2.90
jacobi_max = 3.05

[newton]
max_iterations = 50
tolerance = 1e-12

[bifurcation]
detect_pitchfork = true
detect_period_doubling = true
follow_branches = true
eigenvalue_threshold = 1e-5

[parallel]
num_threads = 0

[system]
mu = 3.0404233870218e-06

[output]
output_dir = "data/periodic_orbit_continuation"
tag = "lyapunov_family"
generate_gnuplot_script = true
save_all_orbits = true
save_bifurcation_info = true
)";
    ofs.close();
    config_files.push_back(default_config);
  }

  std::cout << "Found " << config_files.size() << " config file(s):\n";
  for (const auto& f : config_files) {
    std::cout << "  - " << f << "\n";
  }

  // 設定ファイル読み込み
  std::string config_file = config_files[0];
  AppConfig config;
  if (!LoadConfig(config_file, &config)) {
    std::cerr << "Failed to load configuration file\n";
    return 1;
  }

  PrintConfig(config);

  // 出力ディレクトリ作成
  OutputDirResult output_result = CreateSessionOutputDir(
      config.output_dir, "continuation", output_tag.empty() ? config.output_tag : output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string output_dir = output_result.session_dir;

  // L点位置の確認
  std::cout << "=== Lagrange Point Positions ===\n";
  for (const auto& pt_str : config.starting_points) {
    try {
      LagrangePoint pt = ParseLagrangePoint(pt_str);
      double x_L = ComputeLagrangePointPosition<double>(pt, config.mu);
      std::cout << "  " << pt_str << ": x = " << std::setprecision(12) << x_L << "\n";

      auto coeff = ComputeLyapunovCoefficients<double>(pt, config.mu);
      std::cout << "    c2 = " << coeff.c2 << ", omega_xy = " << coeff.omega_xy
                << ", kappa = " << coeff.kappa << "\n";
      std::cout << "    Estimated period = " << 2.0 * M_PI / coeff.omega_xy << "\n";
    } catch (const std::exception& e) {
      std::cerr << "  " << pt_str << ": Error - " << e.what() << "\n";
    }
  }
  std::cout << "\n";

  // 継続法エンジンを作成
  OrbitContinuator<double> continuator(config.mu, config.continuation);

  // 各ラグランジュ点から継続
  std::vector<LagrangePoint> points;
  for (const auto& pt_str : config.starting_points) {
    try {
      points.push_back(ParseLagrangePoint(pt_str));
    } catch (const std::exception& e) {
      std::cerr << "Skipping invalid point: " << pt_str << " - " << e.what() << "\n";
    }
  }

  if (points.empty()) {
    std::cerr << "No valid starting points\n";
    return 1;
  }

  std::cout << "=== Starting Continuation ===\n\n";

  std::vector<PeriodicOrbitFamily<double>> families;

  if (points.size() > 1 && config.continuation.num_threads != 1) {
    // 複数L点を並列処理
    families = continuator.ContinueMultipleFamilies(points);
  } else {
    // 逐次処理
    for (const auto& pt : points) {
      auto family = continuator.ContinueLyapunovFamily(pt);
      families.push_back(family);
    }
  }

  // 結果を保存
  std::cout << "\n=== Saving Results ===\n";

  for (size_t i = 0; i < families.size(); ++i) {
    const auto& family = families[i];
    std::string family_dir = output_dir + "/" + family.family_name;
    namespace fs = std::filesystem;
    fs::create_directories(family_dir);

    std::cout << "\nFamily: " << family.family_name << "\n";
    std::cout << "  Orbits found: " << family.orbits.size() << "\n";
    std::cout << "  Bifurcations detected: " << family.bifurcations.size() << "\n";
    std::cout << "  Termination reason: " << family.termination_reason << "\n";

    if (config.save_all_orbits && !family.orbits.empty()) {
      std::string csv_path = family_dir + "/orbits.csv";
      ExportFamilyToCsv(family, csv_path);
      std::cout << "  Saved orbits to: " << csv_path << "\n";
    }

    if (config.save_bifurcation_info && !family.bifurcations.empty()) {
      std::string bif_path = family_dir + "/bifurcations.csv";
      ExportBifurcationsToCsv(family, bif_path);
      std::cout << "  Saved bifurcations to: " << bif_path << "\n";
    }

    if (config.generate_gnuplot_script && !family.orbits.empty()) {
      GenerateFamilyGnuplotScript(family, family_dir, "orbits.csv");
      std::cout << "  Generated gnuplot script: " << family_dir << "/plot_family.gp\n";
    }
  }

  // サマリー出力
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

  std::cout << "\n========================================\n";
  std::cout << "  Continuation Complete\n";
  std::cout << "========================================\n";
  std::cout << "  Total families: " << families.size() << "\n";

  int total_orbits = 0;
  int total_bifurcations = 0;
  for (const auto& family : families) {
    total_orbits += static_cast<int>(family.orbits.size());
    total_bifurcations += static_cast<int>(family.bifurcations.size());
  }
  std::cout << "  Total orbits: " << total_orbits << "\n";
  std::cout << "  Total bifurcations: " << total_bifurcations << "\n";
  std::cout << "  Elapsed time: " << duration.count() << " seconds\n";
  std::cout << "  Output directory: " << output_dir << "\n";
  std::cout << "========================================\n";

  return 0;
}
