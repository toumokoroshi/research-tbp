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
 * @brief 設定ファイルから初期条件を読み込む
 */
struct TrajectoryConfig {
  double calc_timestep = 0.0001;
  double time_threshold = 10.0;
  std::vector<my_type::State<double>> initial_coords;
  bool output_sali = false;  ///< SALI計算・出力フラグ
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
    // OUTPUT_SALI
    else if (line.find("OUTPUT_SALI") != std::string::npos) {
      std::string value = TrimString(line.substr(line.find("=") + 1));
      config->output_sali = (value == "1" || value == "true" || value == "TRUE");
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

  // 積分器ラムダ
  auto integrator = [&](const my_type::State<double>& state, double h) {
    return SymplecticStep6thOrder(kMU, state, h);
  };

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
    std::cout << "<>        OUTPUT_SALI: " << (config.output_sali ? "ON" : "OFF") << std::endl;

    // 計算のステップ数
    int num_steps = static_cast<int>(config.time_threshold / config.calc_timestep);
    std::cout << "<>    Total integration steps per trajectory: " << num_steps << std::endl;
    std::cout << "<>" << std::endl;

    auto start_config = std::chrono::system_clock::now();

    // configファイルごとに出力サブフォルダを作成（日付_configファイル名）
    std::string date_str = getcurrent_date();
    std::string subfolder_name = date_str + "_" + config_filename;
    std::string config_output_dir = output_dir + "/" + subfolder_name;

    if (!fs::exists(config_output_dir)) {
      fs::create_directories(config_output_dir);
      std::cout << "<>    Created output subdirectory: " << config_output_dir << std::endl;
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
    if (config.output_sali) {
      ofs << "# Data format: time,x,y,z,vx,vy,vz,sali\n";
    } else {
      ofs << "# Data format: time,x,y,z,vx,vy,vz\n";
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
      if (config.output_sali) {
        ofs << 0.0 << "," << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
            << state.vy << "," << state.vz << "," << 2.0 << "\n";  // SALI初期値=2.0
      } else {
        ofs << 0.0 << "," << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
            << state.vy << "," << state.vz << "\n";
      }

      // 積分ループ
      double current_time = 0.0;
      int output_interval = std::max(1, num_steps / 1000);  // 最大1000点出力

      // SALI計算用の状態（output_sali有効時のみ使用）
      crtbp::SaliState<double> sali_state;
      std::vector<std::pair<double, double>> sali_timeseries;  // 時系列データ

      if (config.output_sali) {
        sali_state.state = crtbp::ConvertToCanonical(state);
        sali_state.w1 = crtbp::CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        sali_state.w2 = crtbp::CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        sali_timeseries.reserve(num_steps / output_interval + 1);
        sali_timeseries.emplace_back(0.0, 2.0);  // 初期値
      }

      for (int step = 0; step < num_steps; ++step) {
        double sali_value = 0.0;

        if (config.output_sali) {
          // SALI用積分器を使用
          crtbp::SymplecticStep6thOrderSALI(kMU, &sali_state, config.calc_timestep);

          // 偏差ベクトル正規化
          sali_state.w1.Normalize();
          sali_state.w2.Normalize();

          // 状態をStateに変換して更新（CanonicalState: q=位置, p=運動量=速度）
          state.x = sali_state.state.qx;
          state.y = sali_state.state.qy;
          state.z = sali_state.state.qz;
          state.vx = sali_state.state.px;
          state.vy = sali_state.state.py;
          state.vz = sali_state.state.pz;

          // SALI計算
          double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
          double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
          sali_value = std::min(norm_plus, norm_minus);
        } else {
          // 通常の積分器を使用
          state = integrator(state, config.calc_timestep);
        }
        current_time += config.calc_timestep;

        // 間引いて出力
        if ((step + 1) % output_interval == 0 || step == num_steps - 1) {
          if (config.output_sali) {
            ofs << current_time << "," << state.x << "," << state.y << "," << state.z << ","
                << state.vx << "," << state.vy << "," << state.vz << "," << sali_value << "\n";
            sali_timeseries.emplace_back(current_time, sali_value);
          } else {
            ofs << current_time << "," << state.x << "," << state.y << "," << state.z << ","
                << state.vx << "," << state.vy << "," << state.vz << "\n";
          }
        }

        // 進捗表示
        if (step % (num_steps / 10 + 1) == 0) {
          double progress = static_cast<double>(step) / num_steps * 100.0;
          std::cout << "\r<>        Progress: " << std::fixed << std::setprecision(1) << progress
                    << "%" << std::flush;
        }
      }
      std::cout << "\r<>        Progress: 100.0%" << std::endl;

      // SALI時系列プロット生成
      if (config.output_sali && !sali_timeseries.empty()) {
        std::string sali_csv_path = config_output_dir + "/" + base_name + "_sali_traj" +
                                    std::to_string(coord_idx + 1) + ".csv";
        std::ofstream sali_ofs(sali_csv_path);
        if (sali_ofs) {
          sali_ofs << std::setprecision(15) << std::fixed;
          sali_ofs << "# SALI Time Series - Trajectory " << (coord_idx + 1) << "\n";
          sali_ofs << "time,sali\n";
          for (const auto& ts : sali_timeseries) {
            sali_ofs << ts.first << "," << ts.second << "\n";
          }
          sali_ofs.close();
          std::cout << "<>        SALI CSV: " << sali_csv_path << std::endl;

          // gnuplotでSALI時系列プロット
          std::string sali_gp_path = config_output_dir + "/" + base_name + "_sali_traj" +
                                     std::to_string(coord_idx + 1) + ".gp";
          std::string sali_eps_path = config_output_dir + "/" + base_name + "_sali_traj" +
                                      std::to_string(coord_idx + 1) + ".eps";
          std::string sali_png_path = config_output_dir + "/" + base_name + "_sali_traj" +
                                      std::to_string(coord_idx + 1) + ".png";
          std::ofstream gp(sali_gp_path);
          if (gp) {
            gp << "set datafile separator ','\n";
            gp << "set datafile commentschars '#'\n";
            gp << "set terminal postscript eps enhanced color font 'Helvetica,14'\n";
            gp << "set output '" << sali_eps_path << "'\n";
            gp << "set xlabel 'Time (non-dim)'\n";
            gp << "set ylabel 'SALI'\n";
            gp << "set title 'SALI Time Series - Trajectory " << (coord_idx + 1) << "'\n";
            gp << "set grid\n";
            gp << "set logscale y\n";
            gp << "set format y '10^{%L}'\n";
            gp << "plot '" << sali_csv_path
               << "' using 1:($2 > 0 ? $2 : 1e-16) with lines lw 2 lc rgb 'blue' title 'SALI'\n";
            gp << "\n";
            gp << "set terminal pngcairo enhanced font 'Helvetica,14' size 1200,600\n";
            gp << "set output '" << sali_png_path << "'\n";
            gp << "replot\n";
            gp.close();

            std::string cmd = "gnuplot \"" + sali_gp_path + "\"";
            int ret = std::system(cmd.c_str());
            if (ret == 0) {
              std::cout << "<>        SALI plot: " << sali_png_path << std::endl;
            }
          }
        }
      }

      // 軌道間に空行を挿入（gnuplotのインデックス機能用）
      // 最後の軌道の後には空行を入れない
      if (coord_idx < config.initial_coords.size() - 1) {
        ofs << "\n\n";  // 2つの空行でgnuplotのインデックスを区切る
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
