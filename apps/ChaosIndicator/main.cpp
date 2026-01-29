/**
 * @file main.cpp
 * @brief ChaosIndicator - velocity_mapからSALI/GALIを計算するアプリケーション
 *
 * stable_manifold appで生成されたvelocity_map.csvを読み込み、
 * 各初期条件に対してSALI/GALI指標を計算する。
 *
 * @note Google C++ Style Guide準拠、Doxygenフォーマット
 */

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <rtbp.hpp>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector>

namespace fs = std::filesystem;

using namespace my_type;
using namespace crtbp;
using namespace param;
using namespace utils;

// ---------------------------------------------------------------------------
// 設定構造体
// ---------------------------------------------------------------------------

/**
 * @brief ChaosIndicator設定構造体
 */
struct ChaosIndicatorConfig {
  // 入力設定
  std::string input_mode = "velocity_map";  ///< "velocity_map" or "coords"
  std::string velocity_map_file;            ///< velocity_map.csvのパス
  int sample_step = 1;                      ///< サンプリング間隔
  int max_points = 0;                       ///< 最大点数 (0=制限なし)

  // 積分設定
  double duration = 100.0;     ///< 積分時間（無次元）
  double timestep = 0.01;      ///< タイムステップ（無次元）
  IntegratorType integrator = IntegratorType::kSymplectic6th;

  // SALI設定
  ChaosIndexType chaos_index = ChaosIndexType::SALI;
  double sali_threshold = 1e-8;    ///< カオス判定閾値
  double output_interval = 1.0;    ///< 時系列出力間隔
  int gali_k = 2;                  ///< GALI次数

  // 出力設定
  std::string output_dir = "data/chaos_indicator";
  bool save_timeseries = false;  ///< 時系列データ出力
  bool summary_only = false;     ///< サマリーのみ出力

  // 物理パラメータ
  double mu = 3.003480575402412e-06;  ///< 質量パラメータ (Sun-Earth)
};

// ---------------------------------------------------------------------------
// velocity_map読み込み用構造体
// ---------------------------------------------------------------------------

/**
 * @brief velocity_mapの1行分のデータ
 */
struct VelocityMapRow {
  State<double> state;       ///< 状態量 (x, y, z, vx, vy, vz)
  double jacobi;             ///< Jacobi定数
  int manifold_idx;          ///< 多様体インデックス
  double time;               ///< 時刻
};

// ---------------------------------------------------------------------------
// 設定読み込み
// ---------------------------------------------------------------------------

/**
 * @brief TOML設定ファイルを読み込む
 * @param[in] filepath 設定ファイルパス
 * @param[out] config 設定構造体
 * @return 成功時true
 */
bool LoadConfig(const std::string& filepath, ChaosIndicatorConfig* config) {
  try {
    TomlConfigParser toml(filepath);

    // 入力設定
    config->input_mode = toml.GetString("input.mode", "velocity_map");
    config->velocity_map_file = toml.GetString("input.velocity_map_file", "");
    config->sample_step = toml.GetInt("input.sample_step", 1);
    config->max_points = toml.GetInt("input.max_points", 0);

    // 積分設定
    config->duration = toml.GetDouble("integration.duration", 100.0);
    config->timestep = toml.GetDouble("integration.timestep", 0.01);

    std::string integrator_str = toml.GetString("integration.integrator", "SYMPLECTIC6");
    if (integrator_str == "SYMPLECTIC6" || integrator_str == "SYMPLECTIC_6TH") {
      config->integrator = IntegratorType::kSymplectic6th;
    } else if (integrator_str == "SYMPLECTIC4" || integrator_str == "SYMPLECTIC_4TH") {
      config->integrator = IntegratorType::kSymplectic4th;
    }

    // SALI設定
    std::string chaos_str = toml.GetString("sali.index_type", "SALI");
    if (chaos_str == "SALI" || chaos_str == "sali") {
      config->chaos_index = ChaosIndexType::SALI;
      config->gali_k = 2;
    } else if (chaos_str == "GALI2" || chaos_str == "gali2") {
      config->chaos_index = ChaosIndexType::GALI;
      config->gali_k = 2;
    } else if (chaos_str == "GALI4" || chaos_str == "gali4") {
      config->chaos_index = ChaosIndexType::GALI;
      config->gali_k = 4;
    } else if (chaos_str == "GALI6" || chaos_str == "gali6") {
      config->chaos_index = ChaosIndexType::GALI;
      config->gali_k = 6;
    }

    config->sali_threshold = toml.GetDouble("sali.threshold", 1e-8);
    config->output_interval = toml.GetDouble("sali.output_interval", 1.0);

    // 出力設定
    config->output_dir = toml.GetString("output.output_dir", "data/chaos_indicator");
    config->save_timeseries = toml.GetBool("output.save_timeseries", false);
    config->summary_only = toml.GetBool("output.summary_only", false);

    // 物理パラメータ
    config->mu = toml.GetDouble("physics.mu", 3.003480575402412e-06);

  } catch (const std::exception& e) {
    std::cerr << "<> !err! Cannot load config: " << e.what() << std::endl;
    return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// velocity_map読み込み
// ---------------------------------------------------------------------------

/**
 * @brief velocity_map.csvを読み込む
 * @param[in] filepath ファイルパス
 * @param[in] sample_step サンプリング間隔
 * @param[in] max_points 最大点数 (0=制限なし)
 * @param[out] data 読み込んだデータ
 * @return 成功時true
 */
bool LoadVelocityMap(const std::string& filepath, int sample_step, int max_points,
                     std::vector<VelocityMapRow>* data) {
  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "<> !err! Cannot open velocity_map file: " << filepath << std::endl;
    return false;
  }

  std::string line;
  int line_count = 0;
  int loaded_count = 0;

  while (std::getline(ifs, line)) {
    // コメント行・空行スキップ
    if (line.empty() || line[0] == '#') {
      continue;
    }

    // サンプリング
    if (line_count % sample_step != 0) {
      line_count++;
      continue;
    }
    line_count++;

    // 最大点数チェック
    if (max_points > 0 && loaded_count >= max_points) {
      break;
    }

    // CSV パース: x,y,z,vx,vy,vz,jacobi,manifold_idx,time
    std::stringstream ss(line);
    std::string token;
    std::vector<double> values;

    while (std::getline(ss, token, ',')) {
      try {
        values.push_back(std::stod(token));
      } catch (...) {
        break;
      }
    }

    if (values.size() >= 7) {
      VelocityMapRow row;
      row.state.x = values[0];
      row.state.y = values[1];
      row.state.z = values[2];
      row.state.vx = values[3];
      row.state.vy = values[4];
      row.state.vz = values[5];
      row.jacobi = values[6];
      row.manifold_idx = (values.size() > 7) ? static_cast<int>(values[7]) : 0;
      row.time = (values.size() > 8) ? values[8] : 0.0;
      data->push_back(row);
      loaded_count++;
    }
  }

  return !data->empty();
}

// ---------------------------------------------------------------------------
// 積分器名取得
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// メイン関数
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  // コマンドライン引数
  CommonArgs args = ParseCommonArgs(argc, argv);

  // ヘッダー出力
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            ChaosIndicator - SALI/GALI Calculator v1.0" << std::endl;
  std::cout << "<>            (velocity_map -> Chaos Indicator)" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n" << std::endl;

  // パス設定
  const std::string kConfigBasePath = CONFIG_DIR;
  const std::string kOutputBasePath = OUTPUT_DIR;
  const std::string kConfigDirPath = kConfigBasePath + "/chaos_indicator/";

  // 設定ファイル検索
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = args.is_continuous;
  std::vector<std::string> config_files =
      DiscoverConfigFilesToml(kConfigDirPath, "chaos_indicator_config", discovery_opts);

  if (config_files.empty()) {
    std::cerr << "<> No config files found in " << kConfigDirPath << std::endl;
    return -1;
  }

  std::cout << "<> Loaded config file: " << config_files[0] << std::endl;

  // 設定読み込み
  ChaosIndicatorConfig config;
  if (!LoadConfig(config_files[0], &config)) {
    return -1;
  }

  // 設定表示
  std::cout << "<>" << std::endl;
  std::cout << "<>****************************************************************" << std::endl;
  std::cout << "<>  [Configuration Summary]" << std::endl;
  std::cout << "<>" << std::endl;
  std::cout << "<>    Input Mode: " << config.input_mode << std::endl;
  std::cout << "<>    Velocity Map: " << config.velocity_map_file << std::endl;
  std::cout << "<>    Sample Step: " << config.sample_step << std::endl;
  std::cout << "<>    Max Points: " << (config.max_points > 0 ? std::to_string(config.max_points) : "unlimited") << std::endl;
  std::cout << "<>    Duration: " << config.duration << std::endl;
  std::cout << "<>    Timestep: " << config.timestep << std::endl;
  std::cout << "<>    Integrator: " << GetIntegratorName(config.integrator) << std::endl;
  std::cout << "<>    SALI Threshold: " << config.sali_threshold << std::endl;
  std::cout << "<>    mu: " << std::setprecision(15) << config.mu << std::endl;
  std::cout << "<>****************************************************************" << std::endl;

  // velocity_map読み込み
  std::vector<VelocityMapRow> velocity_data;
  if (!LoadVelocityMap(config.velocity_map_file, config.sample_step, config.max_points,
                       &velocity_data)) {
    std::cerr << "<> !err! Failed to load velocity_map" << std::endl;
    return -1;
  }
  std::cout << "<> Loaded " << velocity_data.size() << " points from velocity_map" << std::endl;

  // ユーザー確認
  int core_max = omp_get_max_threads();
  int omp_threads = std::max(1, core_max - 1);
  if (!args.skip_wait) {
    WaitForEnter("<> Press Enter to start calculation...");
    std::cout << "<> Enter number of threads (max " << core_max << "): ";
    std::cin >> omp_threads;
    if (omp_threads > core_max) omp_threads = core_max;
    if (omp_threads < 1) omp_threads = 1;
  }
  omp_set_num_threads(omp_threads);
  std::cout << "<> Using " << omp_threads << " threads" << std::endl;

  // 出力ディレクトリ作成
  OutputDirResult output_result =
      CreateSessionOutputDir(kOutputBasePath, "chaos_indicator", args.output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string output_dir = output_result.session_dir;

  // SALI積分器ラムダ
  const double kMU = config.mu;
  auto sali_integrator = [&](SaliState<double>* state, double h) {
    switch (config.integrator) {
      case IntegratorType::kSymplectic6th:
        SymplecticStep6thOrderSALI(kMU, state, h);
        break;
      case IntegratorType::kSymplectic4th:
        SymplecticStep4thOrderSALI(kMU, state, h);
        break;
      default:
        SymplecticStep6thOrderSALI(kMU, state, h);
        break;
    }
  };

  // 出力ファイル準備
  std::string summary_path = output_dir + "/sali_summary.csv";
  std::ofstream summary_ofs(summary_path);
  summary_ofs << std::setprecision(15) << std::fixed;
  summary_ofs << "# ChaosIndicator Summary\n";
  summary_ofs << "# Source: " << config.velocity_map_file << "\n";
  summary_ofs << "# Duration=" << config.duration << ", Timestep=" << config.timestep << "\n";
  summary_ofs << "idx,x0,y0,z0,vx0,vy0,vz0,jacobi,sali_final,threshold_time,is_chaotic\n";

  // 時系列ディレクトリ
  std::string timeseries_dir = output_dir + "/timeseries";
  if (config.save_timeseries) {
    fs::create_directories(timeseries_dir);
  }

  // 計算開始
  auto start_time = std::chrono::system_clock::now();
  int num_steps = static_cast<int>(config.duration / config.timestep);
  int completed = 0;
  int total = static_cast<int>(velocity_data.size());

  std::cout << "<>" << std::endl;
  std::cout << "<> Starting SALI calculation for " << total << " points..." << std::endl;

#pragma omp parallel for schedule(dynamic) shared(completed, summary_ofs)
  for (int idx = 0; idx < total; ++idx) {
    const VelocityMapRow& row = velocity_data[idx];

    // 正準座標に変換
    CanonicalState<double> canonical = ConvertToCanonical(row.state);

    // SALI状態初期化
    SaliState<double> sali_state;
    sali_state.state = canonical;
    sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

    double sali_final = 1.0;
    double threshold_time = config.duration;
    bool threshold_reached = false;

    // 時系列データ
    std::vector<std::pair<double, double>> timeseries;
    if (config.save_timeseries) {
      timeseries.reserve(num_steps);
    }

    // 積分ループ
    for (int step = 0; step < num_steps; ++step) {
      sali_integrator(&sali_state, config.timestep);
      sali_state.w1.Normalize();
      sali_state.w2.Normalize();

      double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
      double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
      sali_final = std::min(norm_plus, norm_minus);

      double current_time = (step + 1) * config.timestep;

      // 時系列記録
      if (config.save_timeseries) {
        timeseries.emplace_back(current_time, sali_final);
      }

      // 閾値判定
      if (!threshold_reached && sali_final < config.sali_threshold) {
        threshold_time = current_time;
        threshold_reached = true;
      }
    }

    // カオス判定
    int is_chaotic = (sali_final < config.sali_threshold) ? 1 : 0;

    // サマリー出力
#pragma omp critical
    {
      summary_ofs << idx << ","
                  << row.state.x << "," << row.state.y << "," << row.state.z << ","
                  << row.state.vx << "," << row.state.vy << "," << row.state.vz << ","
                  << row.jacobi << ","
                  << sali_final << "," << threshold_time << "," << is_chaotic << "\n";

      completed++;
      if (completed % 100 == 0 || completed == total) {
        std::cout << "<>   Progress: " << completed << " / " << total
                  << " (" << (100.0 * completed / total) << "%)" << std::endl;
      }
    }

    // 時系列ファイル出力
    if (config.save_timeseries) {
      std::ostringstream ts_filename;
      ts_filename << timeseries_dir << "/sali_ts_" << idx << ".csv";
      std::ofstream ts_ofs(ts_filename.str());
      ts_ofs << std::setprecision(15) << std::fixed;
      ts_ofs << "# SALI Time Series, idx=" << idx << "\n";
      ts_ofs << "time,sali\n";
      for (const auto& [t, s] : timeseries) {
        ts_ofs << t << "," << s << "\n";
      }
    }
  }

  summary_ofs.close();

  // 完了
  auto end_time = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

  std::cout << "<>" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<> Calculation completed!" << std::endl;
  std::cout << "<>   Total time: " << elapsed << " seconds" << std::endl;
  std::cout << "<>   Output: " << summary_path << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  if (!args.skip_wait) {
    WaitForEnter("<> Press Enter to exit...");
  }

  return 0;
}
