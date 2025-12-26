/**
 * @file main.cpp
 * @brief 円制限三体問題のSALI計算アプリケーション（軌道面メッシュ版）
 *
 * 任意の平面上に同心円状のメッシュを生成し、ヤコビ積分から速度を計算して
 * SALI/GALI指標を計算する。
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
 * @brief シミュレーション設定構造体
 */
struct SimulationConfig {
  State3d<double> mesh_center;     ///< メッシュ中心座標
  State3d<double> plane_norm_vec;  ///< 平面法線ベクトル
  double r_min;                    ///< 最小半径
  double r_max;                    ///< 最大半径
  int radial_div;                  ///< 半径方向分割数
  int base_angular_div;            ///< 基準角度分割数

  double velo_theta;             ///< 接線方向からの速度回転角(deg)
  double escape_dist = 3.0;      ///< 脱出判定距離
  double collision_dist = 1e-5;  ///< 衝突判定距離
  double jacobi;                 ///< ヤコビ積分

  double timestep = 0.00001;                                        ///< 計算タイムステップ
  double t_end = 10.0;                                              ///< 計算終了時刻
  IntegratorType integrator_type = IntegratorType::kSymplectic6th;  ///< 数値積分法

  ChaosIndexType chaos_index_type = ChaosIndexType::SALI;  ///< カオス指標の種類
  int gali_k;                                              ///< GALIの次数

  // 軌道データモード関連
  bool use_trajectory_plane = false;  ///< true: 軌道データから平面算出, false: norm_vector使用
  std::string orbit_data_path;        ///< 軌道データファイルパス
  int phase_step = 1;                 ///< 軌道データの間引きステップ

  // メッシュ設定
  std::string mesh_type = "concentric";  ///< メッシュタイプ: concentric, uniform, rectangular

  // 出力設定
  bool output_gnuplot = true;  ///< gnuplotスクリプト出力
  bool output_png = true;      ///< PNG画像出力
  bool output_eps = true;      ///< EPS画像出力
};

/**
 * @brief 軌道データ1行分（J2000系）
 */
struct OrbitDataRow {
  double time_j2000;       ///< J2000時刻
  State<double> asteroid;  ///< 小惑星状態量（AU, AU/day）
  State<double> earth;     ///< 地球状態量（AU, AU/day）
};

/**
 * @brief SALI出力データ行構造体
 */
struct SaliOutputRow {
  int mesh_num;
  double sali;
  double x, y, z;
  double vx, vy, vz;
  int collision;
  int escape;
  double calc_time;
};

// 前方宣言
bool ParseSaliDataLine(const std::string& line, SaliOutputRow* output_row);
bool WriteSaliOutputsSortedByValue(const std::string& input_filename,
                                   std::string* sorted_output_filename);

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
  int consecutive_empty_count = 0;      // 連続した空行のカウント
  bool second_section_started = false;  // 2番目のセクションに到達したか
  int line_count = 0;

  while (std::getline(ifs, line)) {
    bool is_empty = line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos;

    if (is_empty) {
      // 空行カウントを増やす
      consecutive_empty_count++;
      continue;
    }

    // 非空行に遭遇
    // 2行連続の空行があった場合、2番目のセクションに移行
    if (consecutive_empty_count >= 2 && !second_section_started) {
      second_section_started = true;
    }
    consecutive_empty_count = 0;  // 非空行なのでカウントリセット

    // 2番目のセクションに到達していない場合はスキップ
    if (!second_section_started) {
      continue;
    }

    // 2番目のセクション：データを読み込む
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
 * @brief 回転座標系の軌道データから軌道面法線ベクトルを算出
 *
 * 角運動量ベクトル h = r × v を複数点で計算し、平均化して法線を求める
 *
 * @param[in] rotated_states 回転座標系での状態量ベクター
 * @param[in] center メッシュ中心座標（このまわりの相対位置で計算）
 * @return 軌道面法線ベクトル（正規化済み）
 */
State3d<double> CalcOrbitalPlaneNormal(const std::vector<State<double>>& rotated_states,
                                       const State3d<double>& center) {
  State3d<double> h_sum{0.0, 0.0, 0.0};

  for (const auto& state : rotated_states) {
    // 中心からの相対位置
    double rx = state.x - center.x;
    double ry = state.y - center.y;
    double rz = state.z - center.z;

    // 角運動量 h = r × v
    double hx = ry * state.vz - rz * state.vy;
    double hy = rz * state.vx - rx * state.vz;
    double hz = rx * state.vy - ry * state.vx;

    h_sum.x += hx;
    h_sum.y += hy;
    h_sum.z += hz;
  }

  // 正規化
  double norm = std::sqrt(h_sum.x * h_sum.x + h_sum.y * h_sum.y + h_sum.z * h_sum.z);
  if (norm < 1e-12) {
    // デフォルトでz軸方向を返す
    return {0.0, 0.0, 1.0};
  }

  return {h_sum.x / norm, h_sum.y / norm, h_sum.z / norm};
}

/**
 * @brief 均一密度メッシュを生成（対数スケール半径 + 固定角度分割）
 *
 * r_max/r_minの比が大きい場合でも均一な点密度を維持する。
 * 半径は対数スケールで分割し、角度分割数は全リングで一定。
 *
 * @param center メッシュ中心座標
 * @param plane_normal 平面法線ベクトル
 * @param r_min 最小半径
 * @param r_max 最大半径
 * @param radial_divisions 半径方向分割数
 * @param angular_divisions 角度分割数（全リング共通）
 * @return メッシュ点のベクトル
 */
std::vector<State3d<double>> CreateUniformDensityMesh(const State3d<double>& center,
                                                      const State3d<double>& plane_normal,
                                                      double r_min, double r_max,
                                                      int radial_divisions, int angular_divisions) {
  std::vector<State3d<double>> meshPoints;

  // 法線ベクトルを正規化
  double n_norm = std::sqrt(plane_normal.x * plane_normal.x + plane_normal.y * plane_normal.y +
                            plane_normal.z * plane_normal.z);
  if (n_norm < 1e-12) {
    std::cerr << "<> !err! Plane normal vector is zero" << std::endl;
    return meshPoints;
  }
  State3d<double> n{plane_normal.x / n_norm, plane_normal.y / n_norm, plane_normal.z / n_norm};

  // 平面上の直交基底ベクトルを構築
  State3d<double> u, v;
  if (std::abs(n.z) < 0.9) {
    // z軸とクロス積でuを作成
    u = {-n.y, n.x, 0.0};
  } else {
    // x軸とクロス積でuを作成
    u = {0.0, -n.z, n.y};
  }
  double u_norm = std::sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
  u.x /= u_norm;
  u.y /= u_norm;
  u.z /= u_norm;

  // v = n × u
  v = {n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z, n.x * u.y - n.y * u.x};

  // メモリ予約
  meshPoints.reserve(static_cast<size_t>(radial_divisions * angular_divisions));

  // 対数スケールで半径を分割
  double log_r_min = std::log(r_min);
  double log_r_max = std::log(r_max);
  double log_r_step = (log_r_max - log_r_min) / static_cast<double>(radial_divisions - 1);

  double theta_step = (2.0 * M_PI) / static_cast<double>(angular_divisions);

  for (int ir = 0; ir < radial_divisions; ++ir) {
    // 対数スケールで半径を計算
    double r = std::exp(log_r_min + static_cast<double>(ir) * log_r_step);

    for (int itheta = 0; itheta < angular_divisions; ++itheta) {
      double theta = static_cast<double>(itheta) * theta_step;

      // 平面上の極座標から3次元座標へ変換
      double local_x = r * std::cos(theta);
      double local_y = r * std::sin(theta);

      double x = center.x + local_x * u.x + local_y * v.x;
      double y = center.y + local_x * u.y + local_y * v.y;
      double z = center.z + local_x * u.z + local_y * v.z;

      meshPoints.push_back({x, y, z});
    }
  }

  return meshPoints;
}

/**
 * @brief 矩形グリッドメッシュを生成（均一密度）
 *
 * 軌道面上に均一な矩形グリッドを生成する。
 * mesh_centerを中心に、±r_maxの範囲で格子点を配置。
 *
 * @param center メッシュ中心座標
 * @param plane_normal 平面法線ベクトル
 * @param r_max 最大範囲（中心からの距離）
 * @param grid_size グリッドの片側分割数（総分割数は2*grid_size+1）
 * @return メッシュ点のベクトル
 */
std::vector<State3d<double>> CreateRectangularMesh(const State3d<double>& center,
                                                   const State3d<double>& plane_normal,
                                                   double r_max, int grid_size) {
  std::vector<State3d<double>> meshPoints;

  // 法線ベクトルを正規化
  double n_norm = std::sqrt(plane_normal.x * plane_normal.x + plane_normal.y * plane_normal.y +
                            plane_normal.z * plane_normal.z);
  if (n_norm < 1e-12) {
    std::cerr << "<> !err! Plane normal vector is zero" << std::endl;
    return meshPoints;
  }
  State3d<double> n{plane_normal.x / n_norm, plane_normal.y / n_norm, plane_normal.z / n_norm};

  // 平面上の直交基底ベクトルを構築
  State3d<double> u, v;
  if (std::abs(n.z) < 0.9) {
    u = {-n.y, n.x, 0.0};
  } else {
    u = {0.0, -n.z, n.y};
  }
  double u_norm = std::sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
  u.x /= u_norm;
  u.y /= u_norm;
  u.z /= u_norm;

  // v = n × u
  v = {n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z, n.x * u.y - n.y * u.x};

  // グリッドのステップサイズ
  double step = r_max / static_cast<double>(grid_size);
  int total_points = (2 * grid_size + 1) * (2 * grid_size + 1);
  meshPoints.reserve(static_cast<size_t>(total_points));

  // 矩形グリッド生成
  for (int i = -grid_size; i <= grid_size; ++i) {
    for (int j = -grid_size; j <= grid_size; ++j) {
      double local_x = static_cast<double>(i) * step;
      double local_y = static_cast<double>(j) * step;

      double x = center.x + local_x * u.x + local_y * v.x;
      double y = center.y + local_x * u.y + local_y * v.y;
      double z = center.z + local_x * u.z + local_y * v.z;

      meshPoints.push_back({x, y, z});
    }
  }

  return meshPoints;
}

/**
 * @brief 軌道面上での速度ベクトルを計算
 *
 * mesh_centerを中心とし、plane_norm_vecを法線とする平面上において、
 * mesh_centerとメッシュ点を結ぶ方向に直交する接線方向を基準として、
 * velo_theta_degだけ平面上で回転させた速度ベクトルを返す。
 *
 * @param point メッシュ点の位置
 * @param center メッシュ中心
 * @param plane_normal 平面法線ベクトル
 * @param v_abs 速度の大きさ
 * @param velo_theta_deg 接線方向からの回転角(度)
 * @return 速度ベクトル
 */
State3d<double> CalcVelocityOnPlane(const State3d<double>& point, const State3d<double>& center,
                                    const State3d<double>& plane_normal, double v_abs,
                                    double velo_theta_deg) {
  // 法線ベクトルを正規化
  const double norm_n =
      std::sqrt(plane_normal.x * plane_normal.x + plane_normal.y * plane_normal.y +
                plane_normal.z * plane_normal.z);
  if (norm_n < 1e-12) {
    return {0.0, 0.0, 0.0};
  }
  const State3d<double> n{plane_normal.x / norm_n, plane_normal.y / norm_n,
                          plane_normal.z / norm_n};

  // 中心からメッシュ点への方向ベクトル (r)
  State3d<double> r{point.x - center.x, point.y - center.y, point.z - center.z};
  const double r_norm = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
  if (r_norm < 1e-12) {
    // 中心点の場合は任意方向を返す（または速度0）
    return {0.0, 0.0, 0.0};
  }
  r.x /= r_norm;
  r.y /= r_norm;
  r.z /= r_norm;

  // 接線方向 t = n × r（右手系、反時計回り方向）
  State3d<double> t{n.y * r.z - n.z * r.y, n.z * r.x - n.x * r.z, n.x * r.y - n.y * r.x};
  const double t_norm = std::sqrt(t.x * t.x + t.y * t.y + t.z * t.z);
  if (t_norm < 1e-12) {
    return {0.0, 0.0, 0.0};
  }
  t.x /= t_norm;
  t.y /= t_norm;
  t.z /= t_norm;

  // velo_theta_degで平面上を回転
  // 回転後の速度方向 = cos(theta) * t + sin(theta) * r
  const double theta_rad = velo_theta_deg * M_PI / 180.0;
  const double cos_theta = std::cos(theta_rad);
  const double sin_theta = std::sin(theta_rad);

  State3d<double> v_dir{cos_theta * t.x + sin_theta * r.x, cos_theta * t.y + sin_theta * r.y,
                        cos_theta * t.z + sin_theta * r.z};

  return {v_abs * v_dir.x, v_abs * v_dir.y, v_abs * v_dir.z};
}

/**
 * @brief SALIカラーコンタ図のgnuplotスクリプトを生成
 *
 * @param data_path 入力データファイルパス（CSV形式）
 * @param output_dir 出力ディレクトリ
 * @param base_name ベースファイル名
 * @param output_png PNG出力するか
 * @param output_eps EPS出力するか
 * @param chaos_index_str カオス指標名（SALI/GALI等）
 * @param orbit_data_path 座標変換後の軌道データパス（オプション、空の場合は使用しない）
 */
void GenerateSaliGnuplot(const std::string& data_path, const std::string& output_dir,
                         const std::string& base_name, bool output_png, bool output_eps,
                         const std::string& chaos_index_str,
                         const std::string& orbit_data_path = "") {
  std::string gp_path = output_dir + "/" + base_name + "_contour.gp";
  std::string png_path = output_dir + "/" + base_name + "_contour.png";
  std::string eps_path = output_dir + "/" + base_name + "_contour.eps";
  std::string png_3d_path = output_dir + "/" + base_name + "_3d.png";

  std::ofstream gp(gp_path);
  if (!gp) {
    std::cerr << "<> !warn! Cannot create gnuplot script: " << gp_path << std::endl;
    return;
  }

  gp << "# SALI/GALI 3D Color Contour Plot\n";
  gp << "# Generated by SaliOrbPlane\n\n";
  gp << "set datafile separator ','\n";
  gp << "set datafile commentschars '#'\n\n";

  // カラーバー設定（対数スケール）
  gp << "set logscale cb\n";
  gp << "set cbrange [1e-16:1]\n";
  gp << "set format cb '10^{%L}'\n";
  gp << "set cblabel '" << chaos_index_str << "'\n\n";

  // カラーパレット（青→緑→黄→赤）
  gp << "set palette defined (0 'blue', 0.25 'cyan', 0.5 'green', 0.75 'yellow', 1 'red')\n\n";

  // 軸ラベル
  gp << "set xlabel 'x (non-dim)'\n";
  gp << "set ylabel 'y (non-dim)'\n";
  gp << "set title '" << chaos_index_str << " Distribution on Orbital Plane'\n\n";

  // アスペクト比
  gp << "set size ratio -1\n";
  gp << "set grid\n\n";

  // EPS出力
  if (output_eps) {
    gp << "set terminal postscript eps enhanced color font 'Helvetica,14'\n";
    gp << "set output '" << eps_path << "'\n";
    gp << "plot '" << data_path
       << "' using 3:4:($2 > 0 ? $2 : 1e-16) with points pt 7 ps 0.5 palette notitle\n\n";
  }

  // PNG出力
  if (output_png) {
    gp << "set terminal pngcairo enhanced font 'Helvetica,12' size 1200,1000\n";
    gp << "set output '" << png_path << "'\n";
    gp << "replot\n\n";
  }

  // ================================================
  // 3次元プロット用スクリプト（座標変換後の軌道 + SALIカラーコンタ）
  // escapeフラグが0の点のみSALIカラーコンタを表示
  // ================================================
  gp << "# ================================================\n";
  gp << "# 3D splot: Transformed orbit + SALI color contour\n";
  gp << "# Only non-escaped orbits (escape=0) are displayed with SALI colors\n";
  gp << "# ================================================\n\n";

  gp << "# --- 3D PNG output ---\n";
  gp << "# set terminal pngcairo enhanced font 'Helvetica,12' size 1400,1000\n";
  gp << "# set output '" << png_3d_path << "'\n";
  gp << "# set view 60, 30\n";
  gp << "# set xlabel 'x (non-dim)'\n";
  gp << "# set ylabel 'y (non-dim)'\n";
  gp << "# set zlabel 'z (non-dim)'\n";
  gp << "# set title '" << chaos_index_str << " 3D Distribution (non-escaped orbits only)'\n";
  gp << "# set logscale cb\n";
  gp << "# set cbrange [1e-16:1]\n";
  gp << "# set palette defined (0 'blue', 0.25 'cyan', 0.5 'green', 0.75 'yellow', 1 'red')\n";
  gp << "# set grid\n";
  gp << "# set ticslevel 0\n\n";

  // 軌道データがある場合のプロットコマンド
  if (!orbit_data_path.empty()) {
    gp << "# --- With transformed orbit trajectory ---\n";
    gp << "# splot '" << orbit_data_path << "' using 2:3:4 with lines lw 1.5 lc rgb 'gray50' "
       << "title 'Transformed Orbit', \\\n";
    gp << "#       '" << data_path << "' using 3:4:5:(($10 == 0) ? (($2 > 0) ? $2 : 1e-16) : 1/0) "
       << "with points pt 7 ps 0.5 palette title '" << chaos_index_str << " (non-escaped)'\n\n";
  } else {
    gp << "# --- SALI color contour only (no orbit trajectory) ---\n";
    gp << "# splot '" << data_path << "' using 3:4:5:(($10 == 0) ? (($2 > 0) ? $2 : 1e-16) : 1/0) "
       << "with points pt 7 ps 0.5 palette notitle\n\n";
  }

  gp << "# --- GUI interactive 3D visualization ---\n";
  gp << "# set terminal wxt size 1200,1000    # or qt, x11\n";
  if (!orbit_data_path.empty()) {
    gp << "# splot '" << orbit_data_path << "' using 2:3:4 with lines lw 1.5 lc rgb 'gray50' "
       << "title 'Transformed Orbit', \\\n";
    gp << "#       '" << data_path << "' using 3:4:5:(($10 == 0) ? (($2 > 0) ? $2 : 1e-16) : 1/0) "
       << "with points pt 7 ps 0.5 palette title '" << chaos_index_str << " (non-escaped)'\n";
  } else {
    gp << "# splot '" << data_path << "' using 3:4:5:(($10 == 0) ? (($2 > 0) ? $2 : 1e-16) : 1/0) "
       << "with points pt 7 ps 0.5 palette notitle\n";
  }
  gp << "# pause -1 'Press Enter to close...'\n\n";

  gp.close();
  std::cout << "<>    Gnuplot script: " << gp_path << std::endl;

  // gnuplotを実行
  std::string cmd = "gnuplot \"" + gp_path + "\"";
  int ret = std::system(cmd.c_str());
  if (ret == 0) {
    if (output_eps) {
      std::cout << "<>    EPS generated: " << eps_path << std::endl;
    }
    if (output_png) {
      std::cout << "<>    PNG generated: " << png_path << std::endl;
    }
  } else {
    std::cerr << "<> !warn! gnuplot execution failed (return code: " << ret << ")" << std::endl;
    std::cerr << "<>        You can run the script manually: gnuplot \"" << gp_path << "\""
              << std::endl;
  }
}

int main(int argc, char* argv[]) {
  using namespace param;
  using namespace crtbp;
  using namespace utils;
  namespace fs = std::filesystem;

  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv);
  bool is_continuous = args.is_continuous;
  bool skip_wait = args.skip_wait;
  std::string output_tag = args.output_tag;

  // CMakeから渡されたCONFIG_DIRマクロを使用
  const std::string kConfigFilePath = CONFIG_DIR;
  const std::string kCalcConfigPath = kConfigFilePath + "/sali_orb_plane/";
  std::string calc_config_prefix = "sali_orb_plane";
  std::string astro_param_file = kConfigFilePath + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);

  const double kAU = astro_params.au;                 // astronomical unit in meters
  const double kGMSUN = astro_params.gm_sun;          // heliocentric gravitational constant m3 s-2
  const double kGMEARTH = astro_params.gm_earth;      // geocentric gravitational constant m3 s-2
  const double kMU = kGMEARTH / (kGMEARTH + kGMSUN);  // mu parameter of Earth-Sun
  std::cout << "-" << std::endl;

  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            CRTBP SALI Calculation on Orbital Planes v1.0" << std::endl;
  std::cout << "<>----------------------------------------------------------------"
               "---\n\n"
            << std::endl;

  // 選択されたモードを表示
  if (is_continuous) {
    std::cout << "<>    selected mode : continuous simulation (--continuous)\n" << std::endl;
  } else {
    std::cout << "<> >  selected mode : single simulation (default)" << std::endl;
  }

  if (skip_wait) {
    std::cout << "<>    WaitForEnter : skipped (--no-wait)" << std::endl;
  }

  // 設定ファイルを検索
  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = is_continuous;
  std::vector<std::string> config_file_list =
      DiscoverConfigFilesToml(kCalcConfigPath, calc_config_prefix, discovery_opts);

  if (config_file_list.empty()) {
    std::cerr << "<> No config files found in " << kCalcConfigPath << std::endl;
    return -1;
  }

  std::cout << "<> Loaded config file list:" << std::endl;
  for (const auto& filename : config_file_list) {
    std::cout << "<>    - " << filename << std::endl;
  }

  if (!skip_wait) {
    WaitForEnter();
  }
  std::cout << "<> " << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  int Core_Max = omp_get_max_threads();
  int OMP_Fmax = Core_Max - 1;
  if (!skip_wait) {
    std::cout << "<>  [OpenMP preparation]" << std::endl;
    std::cout << "<> " << std::endl;
    std::cout << "<> On your PC, " << Core_Max
              << " threads can be used for parallel computing employing OMP." << std::endl;
    std::cout << "<> >  " << std::endl;
    std::cout << "<>   * How many threads do you want use for simulation? "
              << "(input an integer)" << std::endl;
    std::cout << "<>   * （※最大コア数を指定すると計算が終わるまでPCが" << std::endl;
    std::cout << "<>   *    激重になるので，最大値-1くらいが良いかも？）> " << std::endl;
    std::cout << "<> >>> ";
    int getcore = 0;
    std::cin >> getcore;
    if (getcore <= Core_Max && getcore > 0) {
      OMP_Fmax = getcore;
      std::cout << "<>" << std::endl;
      std::cout << "<>     >> Number of OMP threads is " << OMP_Fmax << std::endl;
      std::cout << "<>" << std::endl;
    } else {
      OMP_Fmax = Core_Max;
      std::cout << "  <>     >> Your input is INVALID. OMP threads is "
                << "automatically determined as " << OMP_Fmax << std::endl;
    }
    WaitForEnter();
  } else {
    std::cout << "  <>     >> OMP threads is automatically determined as " << OMP_Fmax << std::endl;
  }
  omp_set_num_threads(OMP_Fmax);
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  std::ifstream ifs;
  // SALI積分器
  auto sali_integrator = [&](SaliState<double>* state_ptr, double h) {
    SymplecticStep6thOrderSALI(kMU, state_ptr, h);
  };
  // GALI4積分器
  auto gali4_integrator = [&](GaliState<double, 4>* state_ptr, double h) {
    SymplecticStep6thOrderGALI(kMU, state_ptr, h);
  };
  // GALI6積分器
  auto gali6_integrator = [&](GaliState<double, 6>* state_ptr, double h) {
    SymplecticStep6thOrderGALI(kMU, state_ptr, h);
  };
  //  実行時間の計測
  auto start_ofall = std::chrono::system_clock::now();

  // 実行ごとのセッション出力ディレクトリを作成
  OutputDirResult output_result = CreateSessionOutputDir(OUTPUT_DIR, "sali_orb_plane", output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string session_output_dir = output_result.session_dir;
  std::cout << "<>" << std::endl;
  std::cout << "<>    Session output directory: " << session_output_dir << std::endl;
  std::cout << "<>" << std::endl;

  // --------  configファイルの数だけSALI計算全体を繰り返す-----------------------
  for (const auto& configfilepath : config_file_list) {
    std::cout << "<>        Next config file : " << configfilepath << std::endl;
    ifs.open(configfilepath);
    if (!ifs) {
      std::cerr << "<> !err! config : " << configfilepath << " does NOT EXIST" << std::endl;
      continue;
    }
    double progress = 0;
    auto start = std::chrono::system_clock::now();
    std::string str;
    std::cout << std::setprecision(10);
    std::cout << "<>    loaded config >>>" << std::endl;

    //--------- 設定ファイル読み込み部分---------
    // TomlConfigParserを使用した設定読み込み
    std::unique_ptr<utils::TomlConfigParser> config_ptr;
    try {
      config_ptr = std::make_unique<utils::TomlConfigParser>(configfilepath);
    } catch (const std::exception& e) {
      std::cerr << "\n<> !ERROR! Failed to parse config file: " << configfilepath << std::endl;
      std::cerr << "<> !ERROR! Reason: " << e.what() << std::endl;
      std::cerr << "<> !ERROR! Please check if all required values are properly set." << std::endl;
      std::cerr << "<> !ERROR! Common issues:" << std::endl;
      std::cerr << "<>         - Empty values (e.g., 'jacobi_integral = ')" << std::endl;
      std::cerr << "<>         - Missing quotes around string values" << std::endl;
      std::cerr << "<>         - Invalid number format" << std::endl;
      std::cerr << "<> Skipping this config file...\n" << std::endl;
      ifs.close();
      continue;
    }
    auto& config = *config_ptr;

    // 必須パラメータのバリデーション
    std::vector<std::string> validation_errors;

    // jacobi_integral は必須（特別な値-999.0をセンチネルとして使用）
    const double kSentinelValue = -999.0;
    double jacobi_check = config.GetDouble("boundary.jacobi_integral", kSentinelValue);
    if (jacobi_check == kSentinelValue || !config.HasKey("boundary.jacobi_integral")) {
      validation_errors.push_back("boundary.jacobi_integral is missing or invalid");
    }

    // 軌道データモードかどうかチェック
    bool use_trajectory_plane = config.GetBool("trajectory.use_trajectory_plane", false);

    auto mesh_center_check = config.GetDoubleArray("mesh.center");
    if (mesh_center_check.size() < 3) {
      validation_errors.push_back("mesh.center must have 3 values");
    }

    // use_trajectory_plane == false の場合のみ norm_vector は必須
    if (!use_trajectory_plane) {
      auto plane_norm_vec_check = config.GetDoubleArray("mesh.norm_vector");
      if (plane_norm_vec_check.size() < 3) {
        validation_errors.push_back(
            "mesh.norm_vector must have 3 values (or set trajectory.use_trajectory_plane = true)");
      }
    }

    auto radial_division_check = config.GetDouble("mesh.radial_division");
    if (radial_division_check == 0.0) {
      validation_errors.push_back("mesh.radial_division must be positive value");
    }
    auto base_angular_division_check = config.GetDouble("mesh.base_angular_division");
    if (base_angular_division_check == 0.0) {
      validation_errors.push_back("mesh.base_angular_division must be positive value");
    }

    // バリデーションエラーがある場合
    if (!validation_errors.empty()) {
      std::cerr << "\n<> !ERROR! Config validation failed for: " << configfilepath << std::endl;
      for (const auto& err : validation_errors) {
        std::cerr << "<> !ERROR!   - " << err << std::endl;
      }
      std::cerr << "<> Skipping this config file...\n" << std::endl;
      ifs.close();
      continue;
    }

    SimulationConfig sim_config;

    // 軌道データモードの読み込み
    sim_config.use_trajectory_plane = use_trajectory_plane;
    sim_config.orbit_data_path = config.GetString("trajectory.orbit_data_path", "");
    sim_config.phase_step = config.GetInt("trajectory.phase_step", 1);

    // mesh.center の読み込み
    auto mesh_center_arr = config.GetDoubleArray("mesh.center");
    if (mesh_center_arr.size() >= 3) {
      sim_config.mesh_center = {mesh_center_arr[0], mesh_center_arr[1], mesh_center_arr[2]};
    }

    // plane_norm_vec の読み込み（軌道データモードでない場合のみ）
    if (!sim_config.use_trajectory_plane) {
      auto plane_norm_arr = config.GetDoubleArray("mesh.norm_vector");
      if (plane_norm_arr.size() >= 3) {
        sim_config.plane_norm_vec = {plane_norm_arr[0], plane_norm_arr[1], plane_norm_arr[2]};
      }
    }

    sim_config.radial_div = static_cast<int>(config.GetDouble("mesh.radial_division"));
    sim_config.base_angular_div = static_cast<int>(config.GetDouble("mesh.base_angular_division"));
    sim_config.r_min = config.GetDouble("mesh.r_min");
    sim_config.r_max = config.GetDouble("mesh.r_max");

    sim_config.timestep = config.GetDouble("simulation.calc_timestep", 0.000001);
    sim_config.t_end = config.GetDouble("simulation.time_threshold", 10.0);
    sim_config.velo_theta = config.GetDouble("simulation.velo_theta_deg", 0.0);

    sim_config.escape_dist = config.GetDouble("boundary.escape_dist", 0.03);
    sim_config.collision_dist = config.GetDouble("boundary.collision_dist", 1e-5);
    sim_config.jacobi = config.GetDouble("boundary.jacobi_integral", 3.0008);

    std::string chaos_str = config.GetString("chaos.index_type", "SALI");
    sim_config.chaos_index_type = ChaosIndexType::SALI;
    sim_config.gali_k = 2;
    if (chaos_str == "SALI" || chaos_str == "sali") {
      sim_config.chaos_index_type = ChaosIndexType::SALI;
      sim_config.gali_k = 2;
    } else if (chaos_str == "GALI2" || chaos_str == "gali2") {
      sim_config.chaos_index_type = ChaosIndexType::GALI;
      sim_config.gali_k = 2;
    } else if (chaos_str == "GALI4" || chaos_str == "gali4") {
      sim_config.chaos_index_type = ChaosIndexType::GALI;
      sim_config.gali_k = 4;
    } else if (chaos_str == "GALI6" || chaos_str == "gali6") {
      sim_config.chaos_index_type = ChaosIndexType::GALI;
      sim_config.gali_k = 6;
    }

    // 出力設定の読み込み
    sim_config.output_gnuplot = config.GetBool("output.gnuplot", true);
    sim_config.output_png = config.GetBool("output.png", true);
    sim_config.output_eps = config.GetBool("output.eps", true);

    // メッシュ設定の読み込み
    sim_config.mesh_type = config.GetString("mesh.type", "concentric");

    ifs.close();

    // 座標変換後の軌道データパス（gnuplot用）
    std::string orbit_rotating_path_for_plot;

    // 軌道データモードの場合：軌道データを読み込み、座標変換して法線ベクトルを算出
    if (sim_config.use_trajectory_plane) {
      std::cout << "<>    [Trajectory Plane Mode]" << std::endl;
      std::cout << "<>        Loading orbit data from: " << sim_config.orbit_data_path << std::endl;

      std::vector<OrbitDataRow> orbit_data;
      if (!LoadOrbitData(sim_config.orbit_data_path, sim_config.phase_step, &orbit_data)) {
        std::cerr << "<> !err! Failed to load orbit data, skipping config..." << std::endl;
        continue;
      }
      std::cout << "<>        Loaded " << orbit_data.size() << " orbit data points" << std::endl;

      // J2000 -> 回転座標系に変換
      std::vector<State<double>> rotated_states;
      rotated_states.reserve(orbit_data.size());
      for (const auto& row : orbit_data) {
        State<double> rotated = ConvertInertial2RotatingV2(row.asteroid, row.earth, astro_params);
        rotated_states.push_back(rotated);
      }
      std::cout << "<>        Coordinate conversion completed" << std::endl;

      // ================================================
      // 軌道データの重心を計算してmesh_centerとして使用
      // ================================================
      State3d<double> orbit_centroid{0.0, 0.0, 0.0};
      for (const auto& state : rotated_states) {
        orbit_centroid.x += state.x;
        orbit_centroid.y += state.y;
        orbit_centroid.z += state.z;
      }
      double n_points = static_cast<double>(rotated_states.size());
      orbit_centroid.x /= n_points;
      orbit_centroid.y /= n_points;
      orbit_centroid.z /= n_points;

      // configのmesh_centerを軌道データ重心で上書き
      State3d<double> config_mesh_center = sim_config.mesh_center;  // 元の値を保存（デバッグ用）
      sim_config.mesh_center = orbit_centroid;
      std::cout << "<>        Orbit centroid calculated: (" << orbit_centroid.x << ", "
                << orbit_centroid.y << ", " << orbit_centroid.z << ")" << std::endl;
      std::cout << "<>        (Config mesh_center was: (" << config_mesh_center.x << ", "
                << config_mesh_center.y << ", " << config_mesh_center.z << "), now overridden)"
                << std::endl;

      // 軌道面法線ベクトルを算出（重心を中心として計算）
      sim_config.plane_norm_vec = CalcOrbitalPlaneNormal(rotated_states, sim_config.mesh_center);
      std::cout << "<>        Calculated orbital plane normal: (" << sim_config.plane_norm_vec.x
                << ", " << sim_config.plane_norm_vec.y << ", " << sim_config.plane_norm_vec.z << ")"
                << std::endl;

      // configファイル名からベース名を取得（軌道データ出力用）
      std::string orbit_config_basename = fs::path(configfilepath).stem().string();

      // 座標変換前後の軌道をファイル出力（configごとに別ファイル）
      std::string orbit_j2000_path =
          session_output_dir + "/orbit_j2000_" + orbit_config_basename + ".dat";
      std::string orbit_rotating_path =
          session_output_dir + "/orbit_rotating_" + orbit_config_basename + ".dat";

      std::ofstream ofs_j2000(orbit_j2000_path);
      std::ofstream ofs_rotating(orbit_rotating_path);

      if (ofs_j2000 && ofs_rotating) {
        ofs_j2000 << "# J2000 Inertial Frame Trajectory (Asteroid)\n";
        ofs_j2000 << "# time_j2000, x, y, z, vx, vy, vz\n";
        ofs_rotating << "# Rotating Frame Trajectory (Converted)\n";
        ofs_rotating << "# idx, x, y, z, vx, vy, vz, jacobi\n";

        for (size_t i = 0; i < orbit_data.size(); ++i) {
          const auto& row = orbit_data[i];
          const auto& rotated = rotated_states[i];

          ofs_j2000 << std::setprecision(12) << row.time_j2000 << ", " << row.asteroid.x << ", "
                    << row.asteroid.y << ", " << row.asteroid.z << ", " << row.asteroid.vx << ", "
                    << row.asteroid.vy << ", " << row.asteroid.vz << "\n";

          // 回転座標系でのヤコビ積分を計算
          double jacobi_val = crtbp::calc_jacobi_integral(rotated, kMU);

          ofs_rotating << i << ", " << std::setprecision(12) << rotated.x << ", " << rotated.y
                       << ", " << rotated.z << ", " << rotated.vx << ", " << rotated.vy << ", "
                       << rotated.vz << ", " << jacobi_val << "\n";
        }

        ofs_j2000.close();
        ofs_rotating.close();
        std::cout << "<>        Orbit J2000 output: " << orbit_j2000_path << std::endl;
        std::cout << "<>        Orbit Rotating output: " << orbit_rotating_path << std::endl;

        // gnuplot用に軌道データパスを保存
        orbit_rotating_path_for_plot = orbit_rotating_path;
      }

      // mesh_centerをデータから更新するオプション（最初の変換済み点を使用、またはconfig指定を維持）
      // ここではconfigのmesh_centerを維持（ユーザーが指定した中心を優先）
    }

    std::cout << "<>        MESH CENTER : " << sim_config.mesh_center.x << " ,"
              << sim_config.mesh_center.y << " ," << sim_config.mesh_center.z << std::endl;
    std::cout << "<>        PLANE NORMAL : " << sim_config.plane_norm_vec.x << " ,"
              << sim_config.plane_norm_vec.y << " ," << sim_config.plane_norm_vec.z << std::endl;
    std::cout << "<>        MESH DIVISION RAD : " << sim_config.radial_div << std::endl;
    std::cout << "<>        MESH DIVISION ANG : " << sim_config.base_angular_div << std::endl;
    std::cout << "<>        MESH R MIN : " << sim_config.r_min << std::endl;
    std::cout << "<>        MESH R MAX : " << sim_config.r_max << std::endl;
    std::cout << "<>        CALC TIMESTEP : " << sim_config.timestep << std::endl;
    std::cout << "<>        CALC TIME THRESHOLD : " << sim_config.t_end << std::endl;
    std::cout << "<>        ESCAPE RADIUS : " << sim_config.escape_dist << std::endl;
    std::cout << "<>        COLLISION DIST : " << sim_config.collision_dist << std::endl;
    std::cout << "<>        JACOBI INTEGRAL : " << sim_config.jacobi << std::endl;
    std::cout << "<>        VELO THETA (deg) : " << sim_config.velo_theta << std::endl;
    if (sim_config.use_trajectory_plane) {
      std::cout << "<>        TRAJECTORY MODE : ON (plane normal from orbit data)" << std::endl;
    }

    std::string chaos_index_str;
    switch (sim_config.chaos_index_type) {
      case ChaosIndexType::GALI:
        chaos_index_str = "GALI" + std::to_string(sim_config.gali_k);
        break;
      default:
        chaos_index_str = "SALI";
        break;
    }
    std::cout << "<>        CHAOS INDEX : " << chaos_index_str << std::endl;
    std::cout << "<>    config file read successfully\n" << std::endl;
    if (!is_continuous && !skip_wait) {
      WaitForEnter();
    }

    std::cout << "<>    -- Start SALI calculation for " << configfilepath << std::endl;
    std::cout << "<>        Generating mesh ";
    std::vector<State3d<double>> meshPoints;
    if (sim_config.mesh_type == "uniform") {
      std::cout << "(uniform density mode)" << std::endl;
      meshPoints = CreateUniformDensityMesh(sim_config.mesh_center, sim_config.plane_norm_vec,
                                            sim_config.r_min, sim_config.r_max,
                                            sim_config.radial_div, sim_config.base_angular_div);
    } else if (sim_config.mesh_type == "rectangular") {
      std::cout << "(rectangular grid mode)" << std::endl;
      meshPoints = CreateRectangularMesh(sim_config.mesh_center, sim_config.plane_norm_vec,
                                         sim_config.r_max, sim_config.radial_div);
    } else {
      std::cout << "(concentric mode)" << std::endl;
      meshPoints = CreateConcentricDimensionlessMesh(
          sim_config.mesh_center, sim_config.plane_norm_vec, sim_config.r_min, sim_config.r_max,
          sim_config.radial_div, sim_config.base_angular_div);
    }

    int countt = static_cast<int>(meshPoints.size());
    std::cout << "<>        " << countt << " mesh generated successfully" << std::endl;
    std::cout << "<>        Start calculation" << std::endl;

    // ---------出力ファイル設定---------
    // configファイル名からベース名を取得
    std::string config_basename = fs::path(configfilepath).stem().string();
    std::string filename = session_output_dir + "/" + config_basename + ".dat";
    std::ofstream ofs1(filename);
    if (!ofs1) {
      std::filesystem::path filepath(filename);
      std::filesystem::path dir = filepath.parent_path();
      if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
      }
      ofs1.open(filename);
      if (!ofs1) {
        std::cerr << "Failed to open file even after creating directory: " << filename << std::endl;
        return -1;
      }
    }

    // configファイルを出力ディレクトリにコピー
    std::string config_copy_path =
        session_output_dir + "/" + fs::path(configfilepath).filename().string();
    try {
      fs::copy_file(configfilepath, config_copy_path, fs::copy_options::overwrite_existing);
      std::cout << "<>    Config file copied to: " << config_copy_path << std::endl;
    } catch (const std::exception& e) {
      std::cerr << "<> !warn! Failed to copy config file: " << e.what() << std::endl;
    }

    // ヘッダーを書き込む
    ofs1 << "# MESH CENTER=" << sim_config.mesh_center.x << " " << sim_config.mesh_center.y << " "
         << sim_config.mesh_center.z << std::endl;
    ofs1 << "# PLANE NORMAL=" << sim_config.plane_norm_vec.x << " " << sim_config.plane_norm_vec.y
         << " " << sim_config.plane_norm_vec.z << std::endl;
    ofs1 << "# MESH DIVISION RAD=" << sim_config.radial_div << std::endl;
    ofs1 << "# MESH DIVISION ANG=" << sim_config.base_angular_div << std::endl;
    ofs1 << "# MESH R MIN=" << sim_config.r_min << std::endl;
    ofs1 << "# MESH R MAX=" << sim_config.r_max << std::endl;
    ofs1 << "# CALCULATION TIMESTEP=" << sim_config.timestep << std::endl;
    ofs1 << "# SIMULATION TIME=" << sim_config.t_end << std::endl;
    ofs1 << "# ESCAPE DIST=" << sim_config.escape_dist << std::endl;
    ofs1 << "# COLLISION DIST=" << sim_config.collision_dist << std::endl;
    ofs1 << "# INITIAL JACOBI INTEGRAL=" << sim_config.jacobi << std::endl;
    ofs1 << "# VELO THETA (deg)=" << sim_config.velo_theta << std::endl;
    ofs1 << "# CHAOS INDEX=" << chaos_index_str << std::endl;
    ofs1 << "Time,SALI,x,y,z,vx,vy,vz,collision,escape,calc_time\n";

    // 計算のステップ数
    int num_step = static_cast<int>(sim_config.t_end / sim_config.timestep);
    int totalIterations = static_cast<int>(meshPoints.size());
    // 進捗カウンタ
    int completed_count = 0;
    // 進捗表示間隔（ループ外で計算してオーバーヘッド削減）
    int display_interval = std::max(totalIterations / 100, 1);
    // ファイル書き込み間隔
    constexpr int kWriteInterval = 1000;

// OpenMP並列化ループ
#pragma omp parallel shared(meshPoints, completed_count, totalIterations, progress, ofs1)
    {
      // 各スレッドがSALIの結果を一時的に保存する文字列
      std::stringstream local_output_buffer;
      local_output_buffer << std::fixed << std::setprecision(15);
#pragma omp for schedule(dynamic)
      for (int idx = 0; idx < static_cast<int>(meshPoints.size()); ++idx) {
        const auto& point = meshPoints[idx];
        int mesh_num = idx + 1;
        bool velo_err = false;
        double sali = -1;
        double vx = 0.0, vy = 0.0, vz = 0.0;
        int collision_flag = 0;
        int escape_flag = 0;
        double calc_time = sim_config.t_end;

        double v_abs = calc_v_abs(point, sim_config.jacobi, kMU);

        if (v_abs > 0) {
          State3d<double> velocity =
              CalcVelocityOnPlane(point, sim_config.mesh_center, sim_config.plane_norm_vec, v_abs,
                                  sim_config.velo_theta);
          vx = velocity.x;
          vy = velocity.y;
          vz = velocity.z;
          velo_err = false;
        } else {
          velo_err = true;
          calc_time = 0.0;
        }

        if (!velo_err) {
          State<double> initial_state = {point.x, point.y, point.z, vx, vy, vz};
          CanonicalState<double> canonical_state = ConvertToCanonical(initial_state);

          // 使用するstateポインタ（衝突/脱出判定用）
          CanonicalState<double>* current_state_ptr = nullptr;

          // SALI/GALI状態の初期化
          SaliState<double> sali_state;
          GaliState<double, 4> gali4_state;
          GaliState<double, 6> gali6_state;

          if (sim_config.chaos_index_type == ChaosIndexType::SALI ||
              (sim_config.chaos_index_type == ChaosIndexType::GALI && sim_config.gali_k == 2)) {
            sali_state.state = canonical_state;
            sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
            current_state_ptr = &sali_state.state;
          } else if (sim_config.chaos_index_type == ChaosIndexType::GALI &&
                     sim_config.gali_k == 4) {
            gali4_state.state = canonical_state;
            gali4_state.w[0] = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            gali4_state.w[1] = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
            gali4_state.w[2] = CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
            gali4_state.w[3] = CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
            current_state_ptr = &gali4_state.state;
          } else if (sim_config.chaos_index_type == ChaosIndexType::GALI &&
                     sim_config.gali_k == 6) {
            gali6_state.state = canonical_state;
            gali6_state.w[0] = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            gali6_state.w[1] = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
            gali6_state.w[2] = CanonicalState<double>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
            gali6_state.w[3] = CanonicalState<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
            gali6_state.w[4] = CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
            gali6_state.w[5] = CanonicalState<double>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
            current_state_ptr = &gali6_state.state;
          }

          // 積分ループ
          for (int step = 0; step < num_step; ++step) {
            if (sim_config.chaos_index_type == ChaosIndexType::SALI ||
                (sim_config.chaos_index_type == ChaosIndexType::GALI && sim_config.gali_k == 2)) {
              sali_integrator(&sali_state, sim_config.timestep);
              sali_state.w1.Normalize();
              sali_state.w2.Normalize();
              double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
              double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
              sali = std::min(norm_plus, norm_minus);
              current_state_ptr = &sali_state.state;
            } else if (sim_config.chaos_index_type == ChaosIndexType::GALI &&
                       sim_config.gali_k == 4) {
              gali4_integrator(&gali4_state, sim_config.timestep);
              gali4_state.NormalizeDeviationVectors();
              sali = gali4_state.ComputeGALI();
              current_state_ptr = &gali4_state.state;
            } else if (sim_config.chaos_index_type == ChaosIndexType::GALI &&
                       sim_config.gali_k == 6) {
              gali6_integrator(&gali6_state, sim_config.timestep);
              gali6_state.NormalizeDeviationVectors();
              sali = gali6_state.ComputeGALI();
              current_state_ptr = &gali6_state.state;
            }

            // Check for collision or escape
            if (current_state_ptr != nullptr) {
              double r2 =
                  calc_r2(current_state_ptr->qx, current_state_ptr->qy, current_state_ptr->qz, kMU);
              if (r2 < sim_config.collision_dist) {
                collision_flag = 1;
                calc_time = step * sim_config.timestep;
                break;
              }
              if (r2 > sim_config.escape_dist) {
                escape_flag = 1;
                calc_time = step * sim_config.timestep;
                break;
              }
            }
          }
        }
        local_output_buffer << mesh_num << "," << sali << "," << point.x << "," << point.y << ","
                            << point.z << "," << vx << "," << vy << "," << vz << ","
                            << collision_flag << "," << escape_flag << "," << calc_time << "\n";
        // 一定件数ごとにバッファをファイルに書き込む (排他制御)
        if (idx % kWriteInterval == 0) {
#pragma omp critical(file_write)
          {
            ofs1 << local_output_buffer.str();
            local_output_buffer.str("");  // バッファをクリア
            local_output_buffer.clear();  // ストリーム状態もクリア
          }
        }
#pragma omp atomic
        completed_count++;

        // スレッド0のみが進捗表示（criticalを回避）
        if (omp_get_thread_num() == 0) {
          if (completed_count % display_interval == 0 ||
              completed_count >= totalIterations - OMP_Fmax) {
            double current_progress = static_cast<double>(completed_count) / totalIterations;
            displayProgressBarThreadSafe(current_progress);
          }
        }
      }
      // ループ終了後に残りのバッファを必ず書き込む
#pragma omp critical
      {
        ofs1 << local_output_buffer.str();
      }
    }  // end of parallel region

    // 最終的な進捗バー表示
    displayProgressBarThreadSafe(1.0);
    std::cout << std::endl;

    ofs1.close();
    std::string sorted_output_filename;
    if (WriteSaliOutputsSortedByValue(filename, &sorted_output_filename)) {
      std::cout << "<>    SALI sorted (desc) file: " << sorted_output_filename << std::endl;
    } else {
      std::cout << "<>    SALI sorting was skipped due to parse error." << std::endl;
    }
    std::cout << "<>    Output File:" << filename << std::endl;

    // gnuplotによるカラーコンタ図の生成
    if (sim_config.output_gnuplot) {
      GenerateSaliGnuplot(filename, session_output_dir, config_basename, sim_config.output_png,
                          sim_config.output_eps, chaos_index_str, orbit_rotating_path_for_plot);
    }

    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    auto msec = duration.count() % 1000;
    auto sec = duration.count() / 1000 % 60;
    auto min = duration.count() / 1000 / 60 % 60;
    auto hour = duration.count() / 1000 / 60 / 60;

    std::cout << "<>        elapsed time : " << hour << "h " << min << "m " << sec << "s " << msec
              << "ms" << std::endl;
    std::cout << "<>        Simulation for " << configfilepath << " finished" << std::endl;
  }

  // 実行時間の計測
  auto end = std::chrono::system_clock::now();

  // 実行時間の計算
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_ofall);

  // 時、分、秒、ミリ秒に分解
  auto msec = duration.count() % 1000;
  auto sec = duration.count() / 1000 % 60;
  auto min = duration.count() / 1000 / 60 % 60;
  auto hour = duration.count() / 1000 / 60 / 60;

  std::cout << "<>" << std::endl;
  std::cout << "<>        Calculation finished" << std::endl;
  std::cout << "<>        Total elapsed time : " << hour << "h " << min << "m " << sec << "s "
            << msec << "ms" << std::endl;

  return 0;
}

bool ParseSaliDataLine(const std::string& line, SaliOutputRow* output_row) {
  if (output_row == nullptr) {
    return false;
  }

  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string token;
  while (std::getline(ss, token, ',')) {
    fields.push_back(token);
  }

  if (fields.size() != 11) {
    return false;
  }

  try {
    output_row->mesh_num = static_cast<int>(std::stod(fields[0]));
    output_row->sali = std::stod(fields[1]);
    output_row->x = std::stod(fields[2]);
    output_row->y = std::stod(fields[3]);
    output_row->z = std::stod(fields[4]);
    output_row->vx = std::stod(fields[5]);
    output_row->vy = std::stod(fields[6]);
    output_row->vz = std::stod(fields[7]);
    output_row->collision = std::stoi(fields[8]);
    output_row->escape = std::stoi(fields[9]);
    output_row->calc_time = std::stod(fields[10]);
  } catch (const std::exception&) {
    return false;
  }

  return true;
}

bool WriteSaliOutputsSortedByValue(const std::string& input_filename,
                                   std::string* sorted_output_filename) {
  if (sorted_output_filename == nullptr) {
    return false;
  }

  std::ifstream input(input_filename);
  if (!input.is_open()) {
    std::cerr << "Failed to open file for SALI sorting: " << input_filename << std::endl;
    return false;
  }

  std::vector<std::string> headers;
  std::vector<SaliOutputRow> rows;
  std::string line;
  bool header_finished = false;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    if (!header_finished && (line[0] == '#' || line.rfind("Time", 0) == 0)) {
      headers.push_back(line);
      if (line[0] != '#') {
        header_finished = true;
      }
      continue;
    }
    SaliOutputRow parsed;
    if (ParseSaliDataLine(line, &parsed)) {
      rows.push_back(parsed);
    }
  }

  if (rows.empty()) {
    std::cerr << "No SALI data lines to sort in " << input_filename << std::endl;
    return false;
  }

  std::sort(rows.begin(), rows.end(), [](const SaliOutputRow& lhs, const SaliOutputRow& rhs) {
    // SALIで降順ソート（高い方が上位＝より規則的な軌道）
    return lhs.sali > rhs.sali;
  });

  std::filesystem::path input_path(input_filename);
  std::filesystem::path sorted_path =
      input_path.parent_path() /
      (input_path.stem().string() + "_sorted" + input_path.extension().string());

  std::ofstream output(sorted_path);
  if (!output.is_open()) {
    std::cerr << "Failed to create SALI sorted file: " << sorted_path << std::endl;
    return false;
  }
  output << std::fixed << std::setprecision(15);
  for (const auto& header : headers) {
    output << header << '\n';
  }
  for (const auto& row : rows) {
    output << row.mesh_num << "," << row.sali << "," << row.x << "," << row.y << "," << row.z << ","
           << row.vx << "," << row.vy << "," << row.vz << "," << row.collision << "," << row.escape
           << "," << row.calc_time << "\n";
  }

  *sorted_output_filename = sorted_path.string();
  return true;
}
