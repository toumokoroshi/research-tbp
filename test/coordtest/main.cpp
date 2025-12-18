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
  double sali_lower_limit = 1e-8;                          ///< SALI閾値
  bool output_sali_timeseries = false;                     ///< SALI時系列出力フラグ
  ChaosIndexType chaos_index_type = ChaosIndexType::SALI;  ///< カオス指標の種類
  int gali_k = 2;                                          ///< GALIの偏差ベクトル数 (2, 4, 6)
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
 * @brief 設定ファイルを読み込む
 * @param[in] filepath 設定ファイルパス
 * @param[out] config 設定構造体
 * @return 成功時true
 */
bool LoadTrajectorySaliConfig(const std::string& filepath, TrajectorySaliConfig* config) {
  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "<> !err! Cannot open config file: " << filepath << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;

    size_t eq_pos = line.find('=');
    if (eq_pos == std::string::npos) continue;

    std::string key = line.substr(0, eq_pos);
    std::string value = line.substr(eq_pos + 1);

    // 末尾コメント除去
    size_t comment_pos = value.find("//");
    if (comment_pos != std::string::npos) {
      value = value.substr(0, comment_pos);
    }
    // 空白除去
    while (!value.empty() && (value.back() == ' ' || value.back() == '\t')) {
      value.pop_back();
    }

    if (key == "DATA_DIR") {
      config->data_dir = value;
    } else if (key == "CALC_DURATION_ND") {
      config->calc_duration_nd = std::stod(value);
    } else if (key == "CALC_TIMESTEP_ND") {
      config->calc_timestep_nd = std::stod(value);
    } else if (key == "FORBIDDEN_RADIUS") {
      config->forbidden_radius = std::stod(value);
    } else if (key == "ESCAPE_JUDGE_HILL") {
      config->escape_judge_hill = std::stod(value);
    } else if (key == "DV_MAX") {
      config->dv_max = std::stod(value);
    } else if (key == "DV_STEP") {
      config->dv_step = std::stod(value);
    } else if (key == "PHASE_STEP") {
      config->phase_step = std::stoi(value);
    } else if (key == "INTEGRATOR") {
      if (value == "SYMPLECTIC_6TH") {
        config->integrator = IntegratorType::kSymplectic6th;
      } else if (value == "SYMPLECTIC_4TH") {
        config->integrator = IntegratorType::kSymplectic4th;
      }
    } else if (key == "SALI_LOWER_LIMIT") {
      config->sali_lower_limit = std::stod(value);
    } else if (key == "OUTPUT_SALI_TIMESERIES") {
      config->output_sali_timeseries = (value == "1" || value == "true" || value == "TRUE");
    } else if (key == "CHAOS_INDEX") {
      if (value == "SALI" || value == "sali") {
        config->chaos_index_type = ChaosIndexType::SALI;
        config->gali_k = 2;
      } else if (value == "GALI2" || value == "gali2") {
        config->chaos_index_type = ChaosIndexType::GALI;
        config->gali_k = 2;
      } else if (value == "GALI4" || value == "gali4") {
        config->chaos_index_type = ChaosIndexType::GALI;
        config->gali_k = 4;
      } else if (value == "GALI6" || value == "gali6") {
        config->chaos_index_type = ChaosIndexType::GALI;
        config->gali_k = 6;
      }
    }
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

template <typename ScalarType>
class MyObs {
 private:
  std::vector<std::array<ScalarType, 8>>& history_;
  ScalarType mu_;
  ScalarType JacobiIntegral_;

 public:
  explicit MyObs(std::vector<std::array<ScalarType, 8>>& history, ScalarType mu,
                 ScalarType JacobiIntegral)
      : history_(history), mu_(mu), JacobiIntegral_(JacobiIntegral) {}

  void operator()(const State<ScalarType>& state, ScalarType t) {
    history_.push_back(
        {t, state.x, state.y, state.z, state.vx, state.vy, state.vz, JacobiIntegral_});
  }

  const std::vector<std::array<ScalarType, 8>>& GetHistory() const { return history_; }
};

int main() {
  // ヘッダー出力
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP J2000 -> Rotating Frame Conversion test" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // パス設定
  const std::string kConfigBasePath = CONFIG_DIR;
  const std::string kOutputBasePath = OUTPUT_DIR;
  const std::string kProjectBasePath = PROJECT_BASE_DIR;
  const std::string kConfigFilePath = kConfigBasePath + "/coordtest/coordtest.txt";
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
  std::string data_dir = kProjectBasePath + "/" + config.data_dir;
  std::vector<std::string> orbit_files = DiscoverOrbitFiles(data_dir);
  if (orbit_files.empty()) {
    std::cerr << "<> !err! No orbit files found in: " << data_dir << std::endl;
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
  std::cout << "<>" << std::endl;
  std::cout << "<>    Velocity deceleration steps: "
            << static_cast<int>(config.dv_max / config.dv_step) + 1 << std::endl;
  std::cout << "<>" << std::endl;
  std::cout << "<>****************************************************************" << std::endl;
  std::cout << "<>" << std::endl;

  // ユーザー確認
  WaitForEnter("<>    Press Enter to start simulation...");

  // OpenMP設定
  int core_max = omp_get_max_threads();
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
  std::string output_base_dir = kOutputBasePath + "/coordtest";
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
    std::string csv_path = output_dir + "/" + file_basename + "_raw_converted.csv";
    std::ofstream ofs(csv_path);
    if (!ofs) {
      std::cerr << "<> !err! Cannot create output file: " << csv_path << std::endl;
      continue;
    }

    // CSVヘッダー
    ofs << std::setprecision(15) << std::fixed;
    ofs << "# coordtest Output\n";
    ofs << "# Source: " << orbit_file << "\n";
    ofs << "# CALC_DURATION_ND=" << config.calc_duration_nd << "\n";
    ofs << "# CALC_TIMESTEP_ND=" << config.calc_timestep_nd << "\n";
    ofs << "phase_idx,x0,y0,z0,vx0,vy0,vz0,jacobi\n";

    auto start_file = std::chrono::system_clock::now();
    std::stringstream buf;
    std::stringstream buf2;
    for (int phase_idx = 0; phase_idx < static_cast<int>(orbit_data.size()); ++phase_idx) {
      const OrbitDataRow& row = orbit_data[phase_idx];
      buf2 << phase_idx << "," << row.asteroid.x << "," << row.asteroid.y << "," << row.asteroid.z
           << "," << row.asteroid.vx << "," << row.asteroid.vy << "," << row.asteroid.vz << ","
           << "\n";

      // J2000→CR3BP回転座標変換
      State<double> asteroid_rot = ConvertInertial2Rotating(row.asteroid, row.earth, astro_params);

      buf << phase_idx << "," << asteroid_rot.x << "," << asteroid_rot.y << "," << asteroid_rot.z
          << "," << asteroid_rot.vx << "," << asteroid_rot.vy << "," << asteroid_rot.vz << ","
          << calc_jacobi_integral(asteroid_rot, kMU) << "\n";
    }
    ofs << buf2.str();
    ofs << std::endl;
    ofs << std::endl;
    ofs << buf.str();
    ofs << std::endl;
    ofs << std::endl;
    buf.clear();
    std::cout << "converting completed" << std::endl;

    // orbit_dataの最初の値を使って数値計算
    State<double> asteroid_init_2000 = orbit_data[0].asteroid;
    State<double> earth_init_j2000 = orbit_data[0].earth;
    State<double> asteroid_rot_init =
        ConvertInertial2Rotating(asteroid_init_2000, earth_init_j2000, astro_params);
    double jacobi_init = calc_jacobi_integral(asteroid_rot_init, kMU);

    class EquationOfMotion<double> rtbp_eom(astro_params);
    auto integrator_rk4 = [&](const my_type::State<double>& state_ptr, double time,
                              double h) -> my_type::State<double> {
      return RungeKutta4Step(rtbp_eom, state_ptr, time, h);
    };
    auto integrator_symplectic = [&](const my_type::State<double>& state_ptr, double time,
                                     double h) -> my_type::State<double> {
      return SymplecticStep4thOrder(kMU, state_ptr, h);
    };
    auto integrator_symplectic6 = [&](const my_type::State<double>& state_ptr, double time,
                                      double h) -> my_type::State<double> {
      return SymplecticStep6thOrder(kMU, state_ptr, h);
    };
    std::vector<std::array<double, 8>> history;

    MyObs<double> my_obs(history, kMU, jacobi_init);

    Integrate(asteroid_rot_init, integrator_symplectic6, my_obs, 0.0, config.calc_timestep_nd,
              config.calc_duration_nd / config.calc_timestep_nd);

    for (const auto& h : history) {
      ofs << h[0] << "," << h[1] << "," << h[2] << "," << h[3] << "," << h[4] << "," << h[5] << ","
          << h[6] << "," << h[7] << "\n";
    }
    ofs.close();

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
