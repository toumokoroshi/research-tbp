/**
 * @file main.cpp
 * @brief 円制限三体問題の回転座標系における軌道計算アプリケーション
 *
 * 設定ファイルから初期条件を読み込み、6次シンプレクティック積分で軌道を計算し、
 * CSVとgnuplotプロットを出力する。
 */

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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
};

/**
 * @brief 設定ファイルを解析してTrajectoryConfigを返す
 */
bool LoadTrajectoryConfig(const std::string& filepath, TrajectoryConfig* config) {
  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "<> !err! Cannot open config file: " << filepath << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(ifs, line)) {
    // CALC TIMESTEP
    if (line.find("CALC TIMESTEP") != std::string::npos) {
      config->calc_timestep = std::stod(line.substr(line.find("=") + 1));
    }
    // TIME THRESHOLD
    else if (line.find("TIME THRESHOLD") != std::string::npos) {
      config->time_threshold = std::stod(line.substr(line.find("=") + 1));
    }
    // COORD=
    else if (line.find("COORD=") != std::string::npos) {
      std::string coord_str = line.substr(line.find("=") + 1);
      std::stringstream ss(coord_str);
      double x, y, z, vx, vy, vz;
      if (ss >> x >> y >> z >> vx >> vy >> vz) {
        config->initial_coords.push_back(my_type::State<double>{x, y, z, vx, vy, vz});
      }
    }
  }
  return true;
}

/**
 * @brief gnuplotスクリプトを生成してEPS/PNGを出力する
 */
void GenerateGnuplot(const std::string& csv_path, const std::string& output_dir,
                     const std::string& base_name) {
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
  gp << "plot '" << csv_path << "' using 2:3 with lines title 'Trajectory' lw 2\n";
  gp << "\n";
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

int main() {
  using namespace crtbp;
  using namespace utils;

  // ヘッダー出力 (SALI3dV2スタイル)
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP Trajectory Calculator ver1.0" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // パス設定
  std::string config_base_path = CONFIG_DIR;
  std::string output_base_path = OUTPUT_DIR;
  std::string config_filepath = config_base_path + "/trajectory_calc/trajectory_calc_config.txt";
  std::string output_dir = output_base_path + "/trajectory_calc";

  // 天文定数の読み込み
  std::string astro_param_file = config_base_path + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);

  const double kGMSUN = astro_params.gm_sun;
  const double kGMEARTH = astro_params.gm_earth;
  const double kMU = kGMEARTH / (kGMEARTH + kGMSUN);

  std::cout << "<>    mu parameter: " << std::setprecision(15) << kMU << std::endl;
  std::cout << "<>" << std::endl;

  // 設定ファイル読み込み
  std::cout << "<>    Loading config file: " << config_filepath << std::endl;
  TrajectoryConfig config;
  if (!LoadTrajectoryConfig(config_filepath, &config)) {
    return -1;
  }

  if (config.initial_coords.empty()) {
    std::cerr << "<> !err! No COORD entries found in config file." << std::endl;
    return -1;
  }

  std::cout << "<>        CALC TIMESTEP: " << config.calc_timestep << std::endl;
  std::cout << "<>        TIME THRESHOLD: " << config.time_threshold << std::endl;
  std::cout << "<>        Number of initial conditions: " << config.initial_coords.size()
            << std::endl;
  std::cout << "<>" << std::endl;

  // 出力ディレクトリの作成
  if (!fs::exists(output_dir)) {
    fs::create_directories(output_dir);
    std::cout << "<>    Created output directory: " << output_dir << std::endl;
  }

  // 計算のステップ数
  int num_steps = static_cast<int>(config.time_threshold / config.calc_timestep);
  std::cout << "<>    Total integration steps: " << num_steps << std::endl;
  std::cout << "<>" << std::endl;

  // 積分器ラムダ
  auto integrator = [&](const my_type::State<double>& state, double h) {
    return SymplecticStep6thOrder(kMU, state, h);
  };

  // 全体の実行時間計測
  auto start_total = std::chrono::system_clock::now();

  // 各初期条件について計算
  for (size_t coord_idx = 0; coord_idx < config.initial_coords.size(); ++coord_idx) {
    std::cout << "<>----------------------------------------------------------------" << std::endl;
    std::cout << "<>    Calculating trajectory " << (coord_idx + 1) << " / "
              << config.initial_coords.size() << std::endl;

    my_type::State<double> state = config.initial_coords[coord_idx];
    std::cout << "<>        Initial state: (" << state.x << ", " << state.y << ", " << state.z
              << ", " << state.vx << ", " << state.vy << ", " << state.vz << ")" << std::endl;

    auto start = std::chrono::system_clock::now();

    // 出力ファイル名
    std::string date_str = getcurrent_date();
    std::string base_name = "trajectory_" + std::to_string(coord_idx + 1) + "_" + date_str;
    std::string csv_path = output_dir + "/" + base_name + ".csv";

    // CSVファイルを開く
    std::ofstream ofs(csv_path);
    if (!ofs) {
      std::cerr << "<> !err! Cannot create output file: " << csv_path << std::endl;
      continue;
    }

    // ヘッダー書き込み
    ofs << std::setprecision(15) << std::fixed;
    ofs << "# Trajectory Calculation Output\n";
    ofs << "# Initial State: " << state.x << " " << state.y << " " << state.z << " " << state.vx
        << " " << state.vy << " " << state.vz << "\n";
    ofs << "# CALC TIMESTEP=" << config.calc_timestep << "\n";
    ofs << "# TIME THRESHOLD=" << config.time_threshold << "\n";
    ofs << "# MU=" << kMU << "\n";
    ofs << "time,x,y,z,vx,vy,vz\n";

    // 初期状態を記録
    ofs << 0.0 << "," << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
        << state.vy << "," << state.vz << "\n";

    // 積分ループ
    double current_time = 0.0;
    int output_interval = std::max(1, num_steps / 1000);  // 最大1000点出力

    for (int step = 0; step < num_steps; ++step) {
      state = integrator(state, config.calc_timestep);
      current_time += config.calc_timestep;

      // 間引いて出力
      if ((step + 1) % output_interval == 0 || step == num_steps - 1) {
        ofs << current_time << "," << state.x << "," << state.y << "," << state.z << "," << state.vx
            << "," << state.vy << "," << state.vz << "\n";
      }

      // 進捗表示
      if (step % (num_steps / 10 + 1) == 0) {
        double progress = static_cast<double>(step) / num_steps * 100.0;
        std::cout << "\r<>        Progress: " << std::fixed << std::setprecision(1) << progress
                  << "%" << std::flush;
      }
    }
    std::cout << "\r<>        Progress: 100.0%" << std::endl;

    ofs.close();
    std::cout << "<>        CSV output: " << csv_path << std::endl;

    // gnuplotでプロット生成
    GenerateGnuplot(csv_path, output_dir, base_name);

    // 経過時間
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    auto msec = duration.count() % 1000;
    auto sec = duration.count() / 1000 % 60;
    auto min = duration.count() / 1000 / 60;
    std::cout << "<>        Elapsed time: " << min << "m " << sec << "s " << msec << "ms"
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
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    All calculations finished" << std::endl;
  std::cout << "<>    Total elapsed time: " << min << "m " << sec << "s " << msec << "ms"
            << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  return 0;
}
