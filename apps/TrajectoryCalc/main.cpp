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

/**
 * @brief カオス指標の種類
 */
enum class ChaosIndexType {
  NONE,  ///< カオス指標を計算しない
  SALI,  ///< SALI (K=2)
  GALI   ///< GALI (K可変)
};

/**
 * @brief 数値積分法の種類
 */
enum class IntegratorType {
  kSymplectic4th,  ///< 4次シンプレクティック法
  kSymplectic6th,  ///< 6次シンプレクティック法 (デフォルト)
  kDormandPrince,  ///< ドルマンプリンス法 (DOPRI5)
  kRungeKutta4th   ///< 4次ルンゲクッタ法
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
};

/**
 * @brief 文字列の前後の空白を除去する
 */
std::string TrimString(const std::string& str) {
  const char* whitespace = " \t\r\n";
  size_t start = str.find_first_not_of(whitespace);
  if (start == std::string::npos) {
    return "";
  }
  size_t end = str.find_last_not_of(whitespace);
  return str.substr(start, end - start + 1);
}

/**
 * @brief 設定ファイルを解析してTrajectoryConfigを返す
 */
bool LoadTrajectoryConfig(const std::string& filepath, TrajectoryConfig* config) {
  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "<> !err! Cannot open config file: " << filepath << std::endl;
    return false;
  }

  // デフォルト値をリセット
  config->initial_coords.clear();

  std::string line;
  while (std::getline(ifs, line)) {
    // 行の前後の空白を除去（Windows改行コード対策）
    line = TrimString(line);
    if (line.empty()) {
      continue;
    }

    // CALC TIMESTEP
    if (line.find("CALC TIMESTEP") != std::string::npos) {
      config->calc_timestep = std::stod(TrimString(line.substr(line.find("=") + 1)));
    }
    // TIME THRESHOLD
    else if (line.find("TIME THRESHOLD") != std::string::npos) {
      config->time_threshold = std::stod(TrimString(line.substr(line.find("=") + 1)));
    }
    // COORD= (カンマ区切り形式に対応)
    else if (line.find("COORD=") != std::string::npos) {
      std::string coord_str = TrimString(line.substr(line.find("=") + 1));
      std::stringstream ss(coord_str);
      std::string token;
      std::vector<double> values;
      while (std::getline(ss, token, ',')) {
        std::string trimmed = TrimString(token);
        if (!trimmed.empty()) {
          values.push_back(std::stod(trimmed));
        }
      }
      if (values.size() >= 6) {
        config->initial_coords.push_back(my_type::State<double>{values[0], values[1], values[2],
                                                                values[3], values[4], values[5]});
      }
    }
    // CHAOS_INDEX (NONE, SALI, GALI2, GALI4, GALI6)
    else if (line.find("CHAOS_INDEX") != std::string::npos) {
      std::string value = TrimString(line.substr(line.find("=") + 1));
      if (value == "NONE" || value == "none" || value == "0") {
        config->chaos_index_type = ChaosIndexType::NONE;
      } else if (value == "SALI" || value == "sali") {
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
    // OUTPUT_SALI (後方互換)
    else if (line.find("OUTPUT_SALI") != std::string::npos) {
      std::string value = TrimString(line.substr(line.find("=") + 1));
      if (value == "1" || value == "true" || value == "TRUE") {
        config->chaos_index_type = ChaosIndexType::SALI;
        config->gali_k = 2;
      }
    }
    // INTEGRATOR (SYMPLECTIC4, SYMPLECTIC6, DOPRI, RK4)
    else if (line.find("INTEGRATOR") != std::string::npos) {
      std::string value = TrimString(line.substr(line.find("=") + 1));
      if (value == "SYMPLECTIC4" || value == "symplectic4" || value == "SYMP4") {
        config->integrator_type = IntegratorType::kSymplectic4th;
      } else if (value == "SYMPLECTIC6" || value == "symplectic6" || value == "SYMP6") {
        config->integrator_type = IntegratorType::kSymplectic6th;
      } else if (value == "DOPRI" || value == "dopri" || value == "DORMANDPRINCE" ||
                 value == "dormandprince") {
        config->integrator_type = IntegratorType::kDormandPrince;
      } else if (value == "RK4" || value == "rk4" || value == "RUNGEKUTTA4") {
        config->integrator_type = IntegratorType::kRungeKutta4th;
      }
    }
    // FREQ_ANALYSIS (周波数解析)
    else if (line.find("FREQ_ANALYSIS") != std::string::npos) {
      std::string value = TrimString(line.substr(line.find("=") + 1));
      if (value == "1" || value == "true" || value == "TRUE" || value == "on" || value == "ON") {
        config->enable_freq_analysis = true;
      } else {
        config->enable_freq_analysis = false;
      }
    }
  }
  return true;
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

/**
 * @brief trajectory_calcディレクトリ内のconfigファイル一覧を取得する
 */
std::vector<std::string> GetConfigFileList(const std::string& config_dir) {
  std::vector<std::string> config_files;
  const std::regex pattern("^trajectory_calc_config.*\\.txt$");

  try {
    for (const auto& entry : fs::directory_iterator(config_dir)) {
      if (entry.is_regular_file()) {
        std::string filename = entry.path().filename().string();
        if (std::regex_match(filename, pattern)) {
          config_files.push_back(fs::absolute(entry.path()).string());
        }
      }
    }
  } catch (const fs::filesystem_error& e) {
    std::cerr << "<> !err! Error accessing directory: " << e.what() << std::endl;
  }

  // ファイル名でソート
  std::sort(config_files.begin(), config_files.end());
  return config_files;
}

int main() {
  using namespace crtbp;
  using namespace utils;

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

  // configファイル一覧を取得
  std::vector<std::string> config_files = GetConfigFileList(config_dir);
  if (config_files.empty()) {
    std::cerr << "<> !err! No config files found in: " << config_dir << std::endl;
    std::cerr << "<>        Expected pattern: trajectory_calc_config*.txt" << std::endl;
    return -1;
  }

  std::cout << "<>    Found " << config_files.size() << " config file(s):" << std::endl;
  for (const auto& file : config_files) {
    std::cout << "<>        - " << file << std::endl;
  }
  std::cout << "<>" << std::endl;

  // 出力ディレクトリの作成
  if (!fs::exists(output_dir)) {
    fs::create_directories(output_dir);
    std::cout << "<>    Created output directory: " << output_dir << std::endl;
  }

  // シミュレーション実行ごとのサブフォルダを作成（日時付き）
  std::string run_timestamp = getcurrent_date();
  std::string run_output_dir = output_dir + "/" + run_timestamp + "_run";
  if (!fs::exists(run_output_dir)) {
    fs::create_directories(run_output_dir);
    std::cout << "<>    Created run output directory: " << run_output_dir << std::endl;
  }

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
    std::cout << "<>        INTEGRATOR: " << integrator_str << std::endl;
    std::cout << "<>        FREQ_ANALYSIS: "
              << (config.enable_freq_analysis ? "ENABLED" : "DISABLED") << std::endl;

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
      }

      for (int step = 0; step < num_steps; ++step) {
        double chaos_value = 0.0;

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
                // Boost odeint を使用したDormand-Prince法
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
                runge_kutta_dopri5<StateVec> stepper;
                stepper.do_step(crtbp_eom, y, 0.0, config.calc_timestep);
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
      if (config.chaos_index_type != ChaosIndexType::NONE && !chaos_timeseries.empty()) {
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
          std::string chaos_gp_path = config_output_dir + "/" + base_name + "_" + chaos_index_str +
                                      "_traj" + std::to_string(coord_idx + 1) + ".gp";
          std::string chaos_eps_path = config_output_dir + "/" + base_name + "_" + chaos_index_str +
                                       "_traj" + std::to_string(coord_idx + 1) + ".eps";
          std::string chaos_png_path = config_output_dir + "/" + base_name + "_" + chaos_index_str +
                                       "_traj" + std::to_string(coord_idx + 1) + ".png";
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
               << "' using ($1 > 0 ? $1 : 1/0):($2 > 0 ? $2 : 1e-16) with lines lw 2 lc rgb 'blue' "
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

      // 軌道要素時系列CSV出力
      if (!orb_elem_timeseries.empty()) {
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

      // 軌道間に空行を挿入（gnuplotのインデックス機能用）
      // 最後の軌道の後には空行を入れない
      if (coord_idx < config.initial_coords.size() - 1) {
        ofs << "\n\n";  // 2つの空行でgnuplotのインデックスを区切る
      }

      // ========== 周波数解析 (Laskar's Frequency Map Analysis) ==========
      if (config.enable_freq_analysis) {
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
