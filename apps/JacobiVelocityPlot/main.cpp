/**
 * @file main.cpp
 * @brief 速度変化量とヤコビ積分の関係を計算・可視化するアプリケーション
 *
 * 円制限三体問題において、小惑星の速度変化量の大きさとヤコビ積分の変化量の
 * 関係を計算し、4次元データ（|Δv|, ΔC, C, x）をCSV出力する。
 * gnuplotスクリプトを生成して3Dプロットを作成する。
 *
 * 軸の構成:
 *   - X軸: 速度変化量の大きさ |Δv| [-]
 *   - Y軸: ヤコビ積分の変化量 ΔC [-]
 *   - Z軸: ヤコビ積分 C [-]
 *   - 色: x軸上の距離 x [-]
 */

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
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

/**
 * @brief サンプルデータ1点分の構造体
 */
struct DataPoint {
  double delta_v;   ///< 速度変化量の大きさ |Δv| [-]
  double delta_c;   ///< ヤコビ積分の変化量 ΔC [-]
  double jacobi_c;  ///< ヤコビ積分 C [-]
  double x_pos;     ///< x軸上の位置 x [-]
};

/**
 * @brief 速度の大きさを計算する
 * @param state 状態量
 * @return 速度の大きさ |v| [-]
 */
inline double CalcVelocityMagnitude(const State<double>& state) {
  return std::sqrt(state.vx * state.vx + state.vy * state.vy + state.vz * state.vz);
}

/**
 * @brief サンプルデータを生成する
 *
 * x軸上（y=0, z=0）の異なる位置について、速度を変化させた時の
 * ヤコビ積分の変化を計算する
 *
 * @param n_samples サンプル数
 * @param mu 質量比 mu = m2/(m1+m2)
 * @param perturb_range 速度摂動の最大大きさ [-]
 * @param x_min x座標の最小値 [-]
 * @param x_max x座標の最大値 [-]
 * @param seed 乱数シード
 * @return サンプルデータのベクトル
 */
std::vector<DataPoint> GenerateSampleData(int32_t n_samples, double mu, double perturb_range,
                                          double x_min, double x_max, uint32_t seed = 42) {
  std::vector<DataPoint> data;
  data.reserve(static_cast<size_t>(n_samples));

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist_x(x_min, x_max);
  std::uniform_real_distribution<double> dist_vxy(0.0, 1);
  std::uniform_real_distribution<double> dist_vz(-0.02, 0.02);
  std::uniform_real_distribution<double> dist_perturb(-perturb_range, 0);
  std::uniform_real_distribution<double> dist_perturb_z(-perturb_range * 0.2,
                                                        0);

  for (int32_t i = 0; i < n_samples; ++i) {
    // Generate random x position on x-axis (y=0, z=0)
    const double x = dist_x(rng);
    const double y = 0.0;
    const double z = 0.0;

    // Generate random initial velocity
    const double vx = 0.0;
    const double vy = dist_vxy(rng);
    const double vz = 0.0;

    const State<double> initial_state{x, y, z, vx, vy, vz};

    // Calculate initial Jacobi integral
    double initial_jacobi;
    try {
      initial_jacobi = calc_jacobi_integral(initial_state, mu);
    } catch (...) {
      continue;  // Skip if calculation fails (too close to primary)
    }
    // Apply random velocity perturbation
    const double dvx = 0.0;
    const double dvy = dist_perturb(rng);
    const double dvz = 0.0;
    
    const State<double> perturbed_state{x, y, z, vx + dvx, vy + dvy, vz + dvz};
    
    // Calculate perturbed Jacobi integral
    double perturbed_jacobi;
    try {
      perturbed_jacobi = calc_jacobi_integral(perturbed_state, mu);
    } catch (...) {
      continue;
    }
    if(vy + dvy < 1e-4) {
    continue;
    }
    // WaitForEnter();

    // Compute derived quantities
    const double delta_v = std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
    const double delta_c = perturbed_jacobi - initial_jacobi;

    data.push_back({delta_v, delta_c, initial_jacobi, x});
  }

  return data;
}

/**
 * @brief データをCSVファイルに出力する
 * @param filepath 出力ファイルパス
 * @param data データ点のベクトル
 * @param mu 使用した質量比
 * @return 成功時true
 */
bool ExportToCSV(const std::string& filepath, const std::vector<DataPoint>& data, double mu) {
  std::ofstream ofs(filepath);
  if (!ofs) {
    std::cerr << "<> !err! Cannot create output file: " << filepath << std::endl;
    return false;
  }

  ofs << std::setprecision(10) << std::scientific;
  ofs << "# Jacobi Integral vs Velocity Change 4D Data\n";
  ofs << "# mu = " << mu << "\n";
  ofs << "# Position: y=0, z=0 (on x-axis)\n";
  ofs << "# Columns: delta_v, delta_c, jacobi_c, x_pos\n";
  ofs << "delta_v,delta_c,jacobi_c,x_pos\n";

  for (const auto& point : data) {
    ofs << point.delta_v << "," << point.delta_c << "," << point.jacobi_c << "," << point.x_pos
        << "\n";
  }

  ofs.close();
  return true;
}

/**
 * @brief gnuplotスクリプトを生成する
 * @param script_path スクリプト出力パス
 * @param csv_path CSVファイルパス
 * @param output_png PNGファイル出力パス
 * @return 成功時true
 */
bool GenerateGnuplotScript(const std::string& script_path, const std::string& csv_path,
                           const std::string& output_png) {
  std::ofstream ofs(script_path);
  if (!ofs) {
    std::cerr << "<> !err! Cannot create gnuplot script: " << script_path << std::endl;
    return false;
  }

  ofs << R"(# Gnuplot script for 4D visualization of Jacobi Integral vs Velocity Change
# Auto-generated by JacobiVelocityPlot

set datafile separator ','
set terminal pngcairo size 1200,1000 enhanced font 'Arial,12'
set output ')" << output_png
      << R"('

# Skip header lines
set datafile commentschars "#"

# Labels
set xlabel 'Velocity Change |{/Symbol D}v| [-]' offset 0,-1
set ylabel 'Jacobi Integral Change {/Symbol D}C [-]' offset 0,-1
set zlabel 'Jacobi Integral C [-]' offset 0,0 rotate by 90

set title 'CR3BP: Velocity Change vs Jacobi Integral\n(Color: x position, y=z=0)' font 'Arial,14'

# Color palette for x position
set palette defined (0 "dark-blue", 0.25 "blue", 0.5 "green", 0.75 "yellow", 1 "red")
set cbrange [*:*]
set cblabel 'x position [-]'

# 3D settings
set view 60, 45, 1, 1
set ticslevel 0
set grid

# Plot
splot ')" << csv_path
      << R"(' using 1:2:3:4 with points pt 7 ps 0.5 palette notitle

# Interactive terminal for manual exploration
# Uncomment the following lines to enable interactive mode:
# set terminal wxt size 1200,1000 enhanced
# replot
# pause -1 "Press Enter to exit"
)";

  ofs.close();
  return true;
}

/**
 * @brief gnuplotインタラクティブスクリプトを生成する（GUI操作用）
 */
bool GenerateInteractiveScript(const std::string& script_path, const std::string& csv_path) {
  std::ofstream ofs(script_path);
  if (!ofs) {
    return false;
  }

  ofs << R"(# Interactive Gnuplot script for 4D visualization
# Use mouse to rotate, zoom, and explore the 3D plot

set datafile separator ','
set terminal wxt size 1200,1000 enhanced persist title 'Jacobi Velocity 4D Plot'

# Skip header lines
set datafile commentschars "#"

# Labels
set xlabel 'Velocity Change |{/Symbol D}v| [-]' offset 0,-1
set ylabel 'Jacobi Integral Change {/Symbol D}C [-]' offset 0,-1
set zlabel 'xpos [-]' rotate by 90

set title 'CR3BP: Velocity Change vs Jacobi Integral\n(Color: x position, y=z=0)\n(Use mouse to rotate)' font 'Arial,14'

# Color palette
set palette defined (2.0 "dark-blue", 2.5 "blue", 3.0 "green", 3.5 "yellow", 4.0 "red")
set cblabel 'Jacobi Integral C [-]'
set cbrange [2.0:3.5]


# 3D settings
set view 60, 45, 1, 1
set mouse
set ticslevel 0
set grid

# Initial plot
splot ')" << csv_path
      << R"(' using 1:2:4:3 with points pt 7 ps 0.5 palette notitle

print "Interactive mode enabled!"
print "  - Left mouse drag:  Rotate view"
print "  - Right mouse drag: Zoom"
print "  - Middle mouse:     Pan"
print "  - Press 'q' to quit"

pause -1 "Press Enter or close window to exit"
)";

  ofs.close();
  return true;
}

void PrintUsage(const char* program_name) {
  std::cout << "<> Usage: " << program_name << " [options]\n";
  std::cout << "<> Options:\n";
  std::cout << "    -n <samples>      Number of samples (default: 10000)\n";
  std::cout << "    -p <range>        Velocity perturbation range (default: 1)\n";
  std::cout << "    -xmin <value>     Minimum x position (default: 0.5)\n";
  std::cout << "    -xmax <value>     Maximum x position (default: 1.5)\n";
  std::cout << "    -s <seed>         Random seed (default: 42)\n";
  std::cout << "    -o <output_dir>   Output directory (default: current)\n";
  std::cout << "    -i                Launch interactive gnuplot\n";
  std::cout << "    -h                Show this help\n";
}

int main(int argc, char* argv[]) {
  // Header
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>         Jacobi Velocity Plot Generator ver1.1" << std::endl;
  std::cout << "<>         4D visualization of velocity change and Jacobi integral"
            << std::endl;
  std::cout << "<>         (4th dim: x position on x-axis, y=z=0)"
            << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // Default parameters
  int32_t n_samples = 10000;
  double perturb_range = 0.001;
  double x_min = 0.97;
  double x_max = 1.03;
  uint32_t seed = 42;
  std::string output_dir = ".";
  bool interactive = false;

  // Parse command line arguments
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-n" && i + 1 < argc) {
      n_samples = std::stoi(argv[++i]);
    } else if (arg == "-p" && i + 1 < argc) {
      perturb_range = std::stod(argv[++i]);
    } else if (arg == "-xmin" && i + 1 < argc) {
      x_min = std::stod(argv[++i]);
    } else if (arg == "-xmax" && i + 1 < argc) {
      x_max = std::stod(argv[++i]);
    } else if (arg == "-s" && i + 1 < argc) {
      seed = static_cast<uint32_t>(std::stoul(argv[++i]));
    } else if (arg == "-o" && i + 1 < argc) {
      output_dir = argv[++i];
    } else if (arg == "-i") {
      interactive = true;
    } else if (arg == "-h") {
      PrintUsage(argv[0]);
      return 0;
    }
  }

  // Load astronomical constants
  const std::string kAstroParamFile = std::string(CONFIG_DIR) + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(kAstroParamFile);
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  // Print parameters
  std::cout << "<>    Parameters:" << std::endl;
  std::cout << "<>      mu:            " << std::setprecision(10) << kMU << std::endl;
  std::cout << "<>      n_samples:     " << n_samples << std::endl;
  std::cout << "<>      perturb_range: " << perturb_range << std::endl;
  std::cout << "<>      x_range:       [" << x_min << ", " << x_max << "]" << std::endl;
  std::cout << "<>      seed:          " << seed << std::endl;
  std::cout << "<>      output_dir:    " << output_dir << std::endl;
  std::cout << "<>      interactive:   " << (interactive ? "yes" : "no") << std::endl;
  std::cout << "<>" << std::endl;

  // Create output directory if not exists
  if (!fs::exists(output_dir)) {
    fs::create_directories(output_dir);
  }

  // Generate sample data
  std::cout << "<>    Generating " << n_samples << " samples..." << std::endl;
  std::vector<DataPoint> data = GenerateSampleData(n_samples, kMU, perturb_range, x_min, x_max, seed);
  std::cout << "<>    Generated " << data.size() << " valid data points" << std::endl;

  if (data.empty()) {
    std::cerr << "<> !err! No valid data generated" << std::endl;
    return -1;
  }

  // Compute statistics
  double min_delta_v = data[0].delta_v, max_delta_v = data[0].delta_v;
  double min_delta_c = data[0].delta_c, max_delta_c = data[0].delta_c;
  double min_jacobi = data[0].jacobi_c, max_jacobi = data[0].jacobi_c;
  double min_x_pos = data[0].x_pos, max_x_pos = data[0].x_pos;

  for (const auto& point : data) {
    min_delta_v = std::min(min_delta_v, point.delta_v);
    max_delta_v = std::max(max_delta_v, point.delta_v);
    min_delta_c = std::min(min_delta_c, point.delta_c);
    max_delta_c = std::max(max_delta_c, point.delta_c);
    min_jacobi = std::min(min_jacobi, point.jacobi_c);
    max_jacobi = std::max(max_jacobi, point.jacobi_c);
    min_x_pos = std::min(min_x_pos, point.x_pos);
    max_x_pos = std::max(max_x_pos, point.x_pos);
  }

  std::cout << "<>" << std::endl;
  std::cout << "<>    Data ranges:" << std::endl;
  std::cout << "<>      |Δv|: [" << min_delta_v << ", " << max_delta_v << "]" << std::endl;
  std::cout << "<>      ΔC:   [" << min_delta_c << ", " << max_delta_c << "]" << std::endl;
  std::cout << "<>      C:    [" << min_jacobi << ", " << max_jacobi << "]" << std::endl;
  std::cout << "<>      x:    [" << min_x_pos << ", " << max_x_pos << "]" << std::endl;

  // Output files
  const std::string sub_dir = output_dir + "/data/3dZVC";
  // Ensure output directory exists
  if (!fs::exists(sub_dir)) {
    fs::create_directories(sub_dir);
  }
  const std::string csv_path = sub_dir + "/jacobi_velocity_4d.csv";
  const std::string script_path = sub_dir + "/plot_jacobi_velocity.gp";
  const std::string interactive_script_path = sub_dir + "/plot_jacobi_velocity_interactive.gp";
  const std::string png_path = sub_dir + "/jacobi_velocity_4d_plot.png";

  // Export data
  std::cout << "<>" << std::endl;
  std::cout << "<>    Exporting data to CSV..." << std::endl;
  if (!ExportToCSV(csv_path, data, kMU)) {
    return -1;
  }
  std::cout << "<>      CSV: " << csv_path << std::endl;

  // Generate gnuplot scripts
  std::cout << "<>    Generating gnuplot scripts..." << std::endl;
  if (!GenerateGnuplotScript(script_path, csv_path, png_path)) {
    return -1;
  }
  std::cout << "<>      Script: " << script_path << std::endl;

  if (!GenerateInteractiveScript(interactive_script_path, csv_path)) {
    return -1;
  }
  std::cout << "<>      Interactive: " << interactive_script_path << std::endl;

  // Execute gnuplot to generate PNG
  std::cout << "<>" << std::endl;
  std::cout << "<>    Running gnuplot to generate PNG..." << std::endl;
  std::string gnuplot_cmd = "gnuplot \"" + script_path + "\"";
  int result = std::system(gnuplot_cmd.c_str());
  if (result == 0) {
    std::cout << "<>      PNG: " << png_path << std::endl;
  } else {
    std::cerr << "<> !warn! gnuplot execution failed (is gnuplot installed?)" << std::endl;
  }

  // Launch interactive mode if requested
  if (interactive) {
    std::cout << "<>" << std::endl;
    std::cout << "<>    Launching interactive gnuplot..." << std::endl;
    std::string interactive_cmd = "gnuplot \"" + interactive_script_path + "\"";
    std::system(interactive_cmd.c_str());
  }

  // Summary
  std::cout << "<>" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    Complete!" << std::endl;
  std::cout << "<>    Output files:" << std::endl;
  std::cout << "<>      - " << csv_path << std::endl;
  std::cout << "<>      - " << script_path << std::endl;
  std::cout << "<>      - " << interactive_script_path << std::endl;
  if (result == 0) {
    std::cout << "<>      - " << png_path << std::endl;
  }
  std::cout << "<>" << std::endl;
  std::cout << "<>    To view interactively, run:" << std::endl;
  std::cout << "<>      gnuplot \"" << interactive_script_path << "\"" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  return 0;
}
