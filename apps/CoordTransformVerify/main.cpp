/**
 * @file main.cpp
 * @brief ConvertInertial2RotatingV2 座標変換検証アプリケーション
 *
 * J2000系軌道データの座標変換関数の動作確認を行う。
 * 2つの軌道を比較：
 * 1. 変換軌道: J2000系の軌道の全時刻をそれぞれConvertInertial2RotatingV2で変換
 * 2. プロパゲート軌道: 最初の値のみ変換し、その後はCR3BPとしてプロパゲート
 *
 * @note JacobiIntegralCalc, SaliOrbPlaneを参考に作成
 */

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
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
 * @brief 軌道データ1行分（J2000系）
 */
struct OrbitDataRow {
  double time_j2000;       ///< J2000時刻
  State<double> asteroid;  ///< 小惑星状態量（AU, AU/day）
  State<double> earth;     ///< 地球状態量（AU, AU/day）
};

/**
 * @brief 軌道データファイルを読み込む
 *
 * TrajectorySaliと同様のフォーマット：
 * 時刻 小惑星(x,y,z,vx,vy,vz) 地球(x,y,z,vx,vy,vz) = 13カラム
 *
 * ファイルは2つのセクションに分かれており、2行の連続した空行で区切られている。
 * 2番目のセクション（2行連続空行の後）のデータを読み込む。
 *
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

    // 間引き
    if (line_count % phase_step != 0) {
      line_count++;
      continue;
    }
    line_count++;

    // データパース
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
 * @param[in] data_dir データディレクトリ
 * @return ファイルパスのリスト
 */
std::vector<std::string> DiscoverOrbitFiles(const std::string& data_dir) {
  std::vector<std::string> files;
  if (!fs::exists(data_dir)) {
    std::cerr << "<> !err! DATA_DIR does not exist: " << data_dir << std::endl;
    return files;
  }

  const std::regex pattern(R"(OBT_\d+_?Earth\.txt)");
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
 * @brief gnuplotスクリプトを生成
 */
void GenerateGnuplotScript(const std::string& output_dir, const std::string& base_name,
                           const std::string& converted_csv, const std::string& propagated_csv) {
  std::string gp_path = output_dir + "/" + base_name + "_3d_compare.gp";
  std::string png_path = output_dir + "/" + base_name + "_3d_compare.png";

  std::ofstream gp(gp_path);
  if (!gp) {
    std::cerr << "<> !warn! Cannot create gnuplot script: " << gp_path << std::endl;
    return;
  }

  gp << "# 座標変換検証 3Dプロット\n";
  gp << "# Converted trajectory vs CR3BP propagated trajectory\n\n";
  gp << "set datafile separator ','\n";
  gp << "set datafile commentschars '#'\n\n";

  gp << "set xlabel 'x (non-dim)'\n";
  gp << "set ylabel 'y (non-dim)'\n";
  gp << "set zlabel 'z (non-dim)'\n";
  gp << "set title 'Coordinate Transform Verification'\n\n";

  gp << "set grid\n";
  gp << "set ticslevel 0\n\n";

  // PNG出力
  gp << "set terminal pngcairo enhanced font 'Helvetica,12' size 1400,1000\n";
  gp << "set output '" << png_path << "'\n";
  gp << "set view 60, 30\n\n";

  gp << "splot '" << converted_csv << "' using 2:3:4 with lines lw 1.5 lc rgb 'blue' "
     << "title 'Converted (J2000->Rotating)', \\\n";
  gp << "      '" << propagated_csv << "' using 2:3:4 with lines lw 1.5 lc rgb 'red' "
     << "title 'CR3BP Propagated'\n\n";

  // インタラクティブ表示用
  gp << "# --- Interactive view ---\n";
  gp << "# set terminal wxt size 1200,1000\n";
  gp << "# replot\n";
  gp << "# pause -1 'Press Enter to close...'\n";

  gp.close();
  std::cout << "<>    Gnuplot script: " << gp_path << std::endl;

  // gnuplotを実行
  std::string cmd = "gnuplot \"" + gp_path + "\"";
  int ret = std::system(cmd.c_str());
  if (ret == 0) {
    std::cout << "<>    PNG generated: " << png_path << std::endl;
  } else {
    std::cerr << "<> !warn! gnuplot execution failed. Run manually: gnuplot \"" << gp_path << "\""
              << std::endl;
  }
}

int main(int argc, char* argv[]) {
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>        ConvertInertial2RotatingV2 Verification Tool v1.0" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n";

  // コマンドライン引数チェック
  if (argc < 2) {
    std::cerr << "<> Usage: " << argv[0] << " <orbit_data_file_or_directory>" << std::endl;
    std::cerr << "<>        Processes OBT_*Earth.txt files." << std::endl;
    return -1;
  }

  std::string input_path = argv[1];
  std::cout << "<>    Input: " << input_path << std::endl;

  // 天文定数読み込み
  const std::string kAstroParamFile = std::string(CONFIG_DIR) + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(kAstroParamFile);
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  std::cout << "<>    mu parameter: " << std::setprecision(15) << kMU << std::endl;

  // ファイルリスト取得
  std::vector<std::string> orbit_files;
  if (fs::is_directory(input_path)) {
    orbit_files = DiscoverOrbitFiles(input_path);
  } else if (fs::is_regular_file(input_path)) {
    orbit_files.push_back(input_path);
  } else {
    std::cerr << "<> !err! Invalid input path: " << input_path << std::endl;
    return -1;
  }

  if (orbit_files.empty()) {
    std::cerr << "<> !err! No orbit files found." << std::endl;
    return -1;
  }

  std::cout << "<>    Found " << orbit_files.size() << " orbit file(s)" << std::endl;
  std::cout << "<>" << std::endl;

  // 最初のファイルを処理
  const std::string& orbit_file = orbit_files[0];
  std::string filename = fs::path(orbit_file).filename().string();
  std::string base_name = fs::path(orbit_file).stem().string();
  std::string output_dir = fs::path(orbit_file).parent_path().string();

  std::cout << "<>    Processing: " << filename << std::endl;

  // 軌道データ読み込み
  std::vector<OrbitDataRow> orbit_data;
  if (!LoadOrbitData(orbit_file, 1, &orbit_data)) {
    std::cerr << "<> !err! Failed to load orbit data" << std::endl;
    return -1;
  }

  std::cout << "<>    Loaded " << orbit_data.size() << " data points" << std::endl;

  // ========================================
  // 1. 変換軌道: 全時刻でConvertInertial2RotatingV2を適用
  // ========================================
  std::vector<State<double>> converted_trajectory;
  converted_trajectory.reserve(orbit_data.size());

  double time_offset = orbit_data[0].time_j2000;

  for (const auto& row : orbit_data) {
    State<double> converted = ConvertInertial2RotatingV2(row.asteroid, row.earth, astro_params);
    converted_trajectory.push_back(converted);
  }

  std::cout << "<>    Converted trajectory: " << converted_trajectory.size() << " points"
            << std::endl;

  // ========================================
  // 2. プロパゲート軌道: 最初の値からCR3BPプロパゲート
  // ========================================
  std::vector<State<double>> propagated_trajectory;

  // 時間スケール変換: J2000の日数 -> 無次元時間
  // 1年 = 2π (無次元時間)
  constexpr double kDaysPerYear = 365.25;
  constexpr double kTwoPi = 2.0 * 3.14159265358979323846;

  // 時間間隔の計算（最初の2点から推定）
  double dt_j2000 = 0.0;
  if (orbit_data.size() >= 2) {
    dt_j2000 = orbit_data[1].time_j2000 - orbit_data[0].time_j2000;  // [days]
  } else {
    dt_j2000 = 1.0;  // デフォルト1日
  }

  // 無次元時間へ変換
  double dt_nd = dt_j2000 / kDaysPerYear * kTwoPi;  // [non-dimensional time]

  std::cout << "<>    Time step: " << dt_j2000 << " days = " << dt_nd << " (non-dim)" << std::endl;

  // シンプレクティック積分のための細かいステップ
  constexpr double kIntegrationSubstep = 0.0001;  // 無次元時間での細かいステップ
  int substeps_per_output = static_cast<int>(std::ceil(dt_nd / kIntegrationSubstep));
  double actual_substep = dt_nd / substeps_per_output;

  std::cout << "<>    Substeps per output: " << substeps_per_output << std::endl;

  // 初期状態（変換軌道の最初の点）
  State<double> state = converted_trajectory[0];
  propagated_trajectory.push_back(state);

  // data点数-1回分プロパゲート
  for (size_t i = 1; i < orbit_data.size(); ++i) {
    // 各出力点間をsubstepsで積分
    for (int j = 0; j < substeps_per_output; ++j) {
      state = SymplecticStep6thOrder(kMU, state, actual_substep);
    }
    propagated_trajectory.push_back(state);
  }

  std::cout << "<>    Propagated trajectory: " << propagated_trajectory.size() << " points"
            << std::endl;

  // ========================================
  // 3. 結果出力
  // ========================================
  std::string converted_csv = output_dir + "/" + base_name + "_converted.csv";
  std::string propagated_csv = output_dir + "/" + base_name + "_propagated.csv";

  // 変換軌道出力
  std::ofstream ofs_conv(converted_csv);
  ofs_conv << std::setprecision(15) << std::fixed;
  ofs_conv << "# Converted trajectory (J2000 -> Rotating)\n";
  ofs_conv << "time_nd,x,y,z,vx,vy,vz,jacobi\n";
  for (size_t i = 0; i < converted_trajectory.size(); ++i) {
    double t_nd = (orbit_data[i].time_j2000 - time_offset) / kDaysPerYear * kTwoPi;
    const auto& s = converted_trajectory[i];
    double jacobi = calc_jacobi_integral(s, kMU);
    ofs_conv << t_nd << "," << s.x << "," << s.y << "," << s.z << "," << s.vx << "," << s.vy << ","
             << s.vz << "," << jacobi << "\n";
  }
  ofs_conv.close();
  std::cout << "<>    Output: " << converted_csv << std::endl;

  // プロパゲート軌道出力
  std::ofstream ofs_prop(propagated_csv);
  ofs_prop << std::setprecision(15) << std::fixed;
  ofs_prop << "# CR3BP Propagated trajectory\n";
  ofs_prop << "time_nd,x,y,z,vx,vy,vz,jacobi\n";
  for (size_t i = 0; i < propagated_trajectory.size(); ++i) {
    double t_nd = (orbit_data[i].time_j2000 - time_offset) / kDaysPerYear * kTwoPi;
    const auto& s = propagated_trajectory[i];
    double jacobi = calc_jacobi_integral(s, kMU);
    ofs_prop << t_nd << "," << s.x << "," << s.y << "," << s.z << "," << s.vx << "," << s.vy << ","
             << s.vz << "," << jacobi << "\n";
  }
  ofs_prop.close();
  std::cout << "<>    Output: " << propagated_csv << std::endl;

  // ========================================
  // 4. gnuplotスクリプト生成・実行
  // ========================================
  GenerateGnuplotScript(output_dir, base_name, converted_csv, propagated_csv);

  // ========================================
  // 5. 差異の統計
  // ========================================
  std::cout << "\n<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    Difference Analysis" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  double max_pos_diff = 0.0;
  double max_vel_diff = 0.0;
  double sum_pos_diff = 0.0;
  double sum_vel_diff = 0.0;

  for (size_t i = 0; i < converted_trajectory.size(); ++i) {
    const auto& c = converted_trajectory[i];
    const auto& p = propagated_trajectory[i];

    double dx = c.x - p.x;
    double dy = c.y - p.y;
    double dz = c.z - p.z;
    double pos_diff = std::sqrt(dx * dx + dy * dy + dz * dz);

    double dvx = c.vx - p.vx;
    double dvy = c.vy - p.vy;
    double dvz = c.vz - p.vz;
    double vel_diff = std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

    max_pos_diff = std::max(max_pos_diff, pos_diff);
    max_vel_diff = std::max(max_vel_diff, vel_diff);
    sum_pos_diff += pos_diff;
    sum_vel_diff += vel_diff;
  }

  double avg_pos_diff = sum_pos_diff / converted_trajectory.size();
  double avg_vel_diff = sum_vel_diff / converted_trajectory.size();

  std::cout << std::scientific << std::setprecision(6);
  std::cout << "<>    Position difference (non-dim):" << std::endl;
  std::cout << "<>        Max:     " << max_pos_diff << std::endl;
  std::cout << "<>        Average: " << avg_pos_diff << std::endl;
  std::cout << "<>    Velocity difference (non-dim):" << std::endl;
  std::cout << "<>        Max:     " << max_vel_diff << std::endl;
  std::cout << "<>        Average: " << avg_vel_diff << std::endl;

  std::cout << "\n<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    Verification complete!" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  return 0;
}
