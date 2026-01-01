/**
 * @file main.cpp
 * @brief 小惑星軌道と交差する周期/準周期軌道探索アプリケーション
 * @details
 *   Phase 1: 小惑星軌道からポアンカレ断面交差点を抽出
 *   Phase 2: 交差点近傍で周期軌道を探索（Newton-Raphson法）
 *   Phase 3: 元の軌道との交点検証
 * @date 2025-12-29
 */

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector>

#include "periodic_orbit.hpp"
#include "rtbp.hpp"

namespace fs = std::filesystem;
using namespace crtbp;
using namespace utils;
using namespace my_type;
using namespace periodic_orbit;

// ---------------------------------------------------------------------------
// 断面交差方向の列挙型（PoincareMapから移植）
// ---------------------------------------------------------------------------
enum class CrossingDirection {
  kPositive,  // 正方向のみ（変数が負→正）
  kNegative,  // 負方向のみ（変数が正→負）
  kBoth       // 両方向
};

// ---------------------------------------------------------------------------
// コンフィグ構造体
// ---------------------------------------------------------------------------
struct IntersectionConfig {
  // 軌道データ設定
  std::string orbit_data_folder;             ///< 小惑星軌道データフォルダパス
  std::string orbit_file_pattern = "*.txt";  ///< ファイルパターン（例: "*.txt", "OBT_*Earth.txt"）
  int phase_step = 1;                        ///< 軌道データの間引きステップ

  // 初期条件生成方式
  std::string ic_source =
      "trajectory";          ///< 初期条件ソース: "trajectory"（軌道点）, "mesh"（メッシュ）
  int trajectory_step = 10;  ///< 軌道点を使う場合の間引きステップ

  // メッシュ生成設定（ic_source="mesh"の場合のみ使用）
  std::string mesh_type = "concentric";  ///< メッシュタイプ: "concentric", "rectangular"
  double r_min = 0.001;                  ///< 最小半径（無次元）
  double r_max = 0.02;                   ///< 最大半径（無次元）
  int n_radial = 10;                     ///< 半径方向分割数
  int n_angular = 36;                    ///< 角度方向分割数
  int grid_size = 10;                    ///< 矩形メッシュの片側分割数

  // ターゲットヤコビ積分設定
  double target_jacobi = 3.0;         ///< ターゲットヤコビ積分値
  double velocity_angle_min = 0.0;    ///< 速度方向角度の最小値（度）
  double velocity_angle_max = 0.0;    ///< 速度方向角度の最大値（度）
  double velocity_angle_step = 10.0;  ///< 速度方向角度のステップ（度）

  // ポアンカレ断面設定（Phase 2用）
  std::string section_var = "y";  ///< 断面変数 (x, y, z)
  double section_value = 0.0;     ///< 断面の値
  CrossingDirection crossing_direction = CrossingDirection::kBoth;

  // Newton-Raphson設定
  int newton_max_iterations = 50;
  double newton_tolerance = 1e-12;
  double max_integration_time = 100.0;
  double dt = 0.0001;

  // 交点検証設定
  double intersection_threshold = 0.001;  ///< 交差判定閾値（無次元距離）
  double verification_time = 500.0;       ///< 検証積分時間

  // 積分器設定
  IntegratorType integrator_type = IntegratorType::kSymplectic6th;

  // 出力設定
  std::string output_tag = "";
  bool enable_phase2 = true;  ///< Phase2（周期軌道精密化）を実行
  bool enable_phase3 = true;  ///< Phase3（交点検証）を実行
};

// ---------------------------------------------------------------------------
// 軌道データ1行分（J2000系）
// ---------------------------------------------------------------------------
struct OrbitDataRow {
  double time_j2000;       ///< J2000時刻
  State<double> asteroid;  ///< 小惑星状態量（AU, AU/day）
  State<double> earth;     ///< 地球状態量（AU, AU/day）
};

// ---------------------------------------------------------------------------
// ポアンカレ断面交差点
// ---------------------------------------------------------------------------
struct PoincareCrossing {
  double time;                           ///< 交差時刻（無次元）
  State<double> state;                   ///< 交差点での状態
  double jacobi;                         ///< ヤコビ積分
  int crossing_index;                    ///< 交差インデックス
  CrossingDirection detected_direction;  ///< 検出された交差方向
};

// ---------------------------------------------------------------------------
// 交点情報
// ---------------------------------------------------------------------------
struct IntersectionEvent {
  int candidate_idx;              ///< 候補軌道インデックス
  double time;                    ///< 交点検出時刻
  State<double> candidate_state;  ///< 候補軌道側の状態
  State<double> reference_state;  ///< 参照軌道側の状態
  double distance;                ///< 2軌道間距離
};

// ---------------------------------------------------------------------------
// 変数名からインデックスを取得
// ---------------------------------------------------------------------------
int GetVarIndex(const std::string& var_name) {
  if (var_name == "x") return 0;
  if (var_name == "y") return 1;
  if (var_name == "z") return 2;
  if (var_name == "vx") return 3;
  if (var_name == "vy") return 4;
  if (var_name == "vz") return 5;
  return -1;
}

// ---------------------------------------------------------------------------
// 状態ベクトルから指定インデックスの値を取得
// ---------------------------------------------------------------------------
double GetStateValue(const State<double>& state, int index) {
  switch (index) {
    case 0:
      return state.x;
    case 1:
      return state.y;
    case 2:
      return state.z;
    case 3:
      return state.vx;
    case 4:
      return state.vy;
    case 5:
      return state.vz;
    default:
      return 0.0;
  }
}

// ---------------------------------------------------------------------------
// 軌道面法線ベクトルを計算（SaliOrbPlaneから移植）
// ---------------------------------------------------------------------------
State3d<double> CalcOrbitalPlaneNormal(const std::vector<State<double>>& rotated_states,
                                       const State3d<double>& center) {
  State3d<double> h_sum{0.0, 0.0, 0.0};
  for (const auto& state : rotated_states) {
    double rx = state.x - center.x;
    double ry = state.y - center.y;
    double rz = state.z - center.z;
    double hx = ry * state.vz - rz * state.vy;
    double hy = rz * state.vx - rx * state.vz;
    double hz = rx * state.vy - ry * state.vx;
    h_sum.x += hx;
    h_sum.y += hy;
    h_sum.z += hz;
  }
  double norm = std::sqrt(h_sum.x * h_sum.x + h_sum.y * h_sum.y + h_sum.z * h_sum.z);
  if (norm < 1e-12) {
    return {0.0, 0.0, 1.0};
  }
  return {h_sum.x / norm, h_sum.y / norm, h_sum.z / norm};
}

// ---------------------------------------------------------------------------
// 同心円メッシュを生成（SaliOrbPlaneから移植）
// ---------------------------------------------------------------------------
std::vector<State3d<double>> CreateConcentricMesh(const State3d<double>& center,
                                                  const State3d<double>& plane_normal, double r_min,
                                                  double r_max, int radial_divisions,
                                                  int angular_divisions) {
  std::vector<State3d<double>> meshPoints;
  double n_norm = std::sqrt(plane_normal.x * plane_normal.x + plane_normal.y * plane_normal.y +
                            plane_normal.z * plane_normal.z);
  if (n_norm < 1e-12) return meshPoints;
  State3d<double> n{plane_normal.x / n_norm, plane_normal.y / n_norm, plane_normal.z / n_norm};

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
  v = {n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z, n.x * u.y - n.y * u.x};

  meshPoints.reserve(radial_divisions * angular_divisions);
  double log_r_min = std::log(r_min);
  double log_r_max = std::log(r_max);
  double log_r_step = (log_r_max - log_r_min) / (radial_divisions - 1);
  double theta_step = (2.0 * M_PI) / angular_divisions;

  for (int ir = 0; ir < radial_divisions; ++ir) {
    double r = std::exp(log_r_min + ir * log_r_step);
    for (int itheta = 0; itheta < angular_divisions; ++itheta) {
      double theta = itheta * theta_step;
      double local_x = r * std::cos(theta);
      double local_y = r * std::sin(theta);
      meshPoints.push_back({center.x + local_x * u.x + local_y * v.x,
                            center.y + local_x * u.y + local_y * v.y,
                            center.z + local_x * u.z + local_y * v.z});
    }
  }
  return meshPoints;
}

// ---------------------------------------------------------------------------
// 矩形メッシュを生成（SaliOrbPlaneから移植）
// ---------------------------------------------------------------------------
std::vector<State3d<double>> CreateRectangularMesh(const State3d<double>& center,
                                                   const State3d<double>& plane_normal,
                                                   double r_max, int grid_size) {
  std::vector<State3d<double>> meshPoints;
  double n_norm = std::sqrt(plane_normal.x * plane_normal.x + plane_normal.y * plane_normal.y +
                            plane_normal.z * plane_normal.z);
  if (n_norm < 1e-12) return meshPoints;
  State3d<double> n{plane_normal.x / n_norm, plane_normal.y / n_norm, plane_normal.z / n_norm};

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
  v = {n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z, n.x * u.y - n.y * u.x};

  double step = r_max / grid_size;
  for (int i = -grid_size; i <= grid_size; ++i) {
    for (int j = -grid_size; j <= grid_size; ++j) {
      double local_x = i * step;
      double local_y = j * step;
      meshPoints.push_back({center.x + local_x * u.x + local_y * v.x,
                            center.y + local_x * u.y + local_y * v.y,
                            center.z + local_x * u.z + local_y * v.z});
    }
  }
  return meshPoints;
}

// ---------------------------------------------------------------------------
// 軌道面上での速度ベクトルを計算
// ---------------------------------------------------------------------------
State3d<double> CalcVelocityOnPlane(const State3d<double>& point, const State3d<double>& center,
                                    const State3d<double>& plane_normal, double v_abs,
                                    double velo_theta_deg) {
  const double norm_n =
      std::sqrt(plane_normal.x * plane_normal.x + plane_normal.y * plane_normal.y +
                plane_normal.z * plane_normal.z);
  if (norm_n < 1e-12) return {0.0, 0.0, 0.0};
  const State3d<double> n{plane_normal.x / norm_n, plane_normal.y / norm_n,
                          plane_normal.z / norm_n};

  State3d<double> r{point.x - center.x, point.y - center.y, point.z - center.z};
  const double r_norm = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
  if (r_norm < 1e-12) return {0.0, 0.0, 0.0};
  r.x /= r_norm;
  r.y /= r_norm;
  r.z /= r_norm;

  State3d<double> t{n.y * r.z - n.z * r.y, n.z * r.x - n.x * r.z, n.x * r.y - n.y * r.x};
  const double t_norm = std::sqrt(t.x * t.x + t.y * t.y + t.z * t.z);
  if (t_norm < 1e-12) return {0.0, 0.0, 0.0};
  t.x /= t_norm;
  t.y /= t_norm;
  t.z /= t_norm;

  const double theta_rad = velo_theta_deg * M_PI / 180.0;
  const double cos_theta = std::cos(theta_rad);
  const double sin_theta = std::sin(theta_rad);
  State3d<double> v_dir{cos_theta * t.x + sin_theta * r.x, cos_theta * t.y + sin_theta * r.y,
                        cos_theta * t.z + sin_theta * r.z};
  return {v_abs * v_dir.x, v_abs * v_dir.y, v_abs * v_dir.z};
}

// ---------------------------------------------------------------------------
// ヤコビ積分から速度の大きさを計算
// ---------------------------------------------------------------------------
double CalcVelocityFromJacobi(double target_jacobi, double x, double y, double z, double mu) {
  double U = calc_potential_U(x, y, z, mu);
  double v_squared = 2.0 * U - target_jacobi;
  if (v_squared < 0) return -1.0;  // 禁止領域
  return std::sqrt(v_squared);
}

// ---------------------------------------------------------------------------
// 軌道データファイルを読み込む（SaliOrbPlaneから移植）
// ---------------------------------------------------------------------------
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

    if (line_count % phase_step != 0) {
      line_count++;
      continue;
    }
    line_count++;

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

// ---------------------------------------------------------------------------
// 指定フォルダ内の軌道データファイルを探索
// ---------------------------------------------------------------------------
std::vector<std::string> DiscoverOrbitDataFiles(const std::string& folder_path,
                                                const std::string& pattern) {
  std::vector<std::string> files;

  if (!fs::exists(folder_path) || !fs::is_directory(folder_path)) {
    std::cerr << "<> !err! Folder does not exist: " << folder_path << std::endl;
    return files;
  }

  // パターンをglob形式からregex形式に変換
  // 例: "*.txt" -> ".*\\.txt$", "OBT_*Earth.txt" -> "OBT_.*Earth\\.txt$"
  std::string regex_pattern = pattern;
  // エスケープが必要な文字を処理
  std::string result;
  for (char c : regex_pattern) {
    switch (c) {
      case '.':
        result += "\\.";
        break;
      case '*':
        result += ".*";
        break;
      case '?':
        result += ".";
        break;
      default:
        result += c;
        break;
    }
  }
  regex_pattern = "^" + result + "$";

  std::regex file_regex(regex_pattern, std::regex::icase);

  // フォルダ内を再帰的に探索
  for (const auto& entry : fs::recursive_directory_iterator(folder_path)) {
    if (entry.is_regular_file()) {
      std::string filename = entry.path().filename().string();
      if (std::regex_match(filename, file_regex)) {
        files.push_back(entry.path().string());
      }
    }
  }

  // ファイル名でソート
  std::sort(files.begin(), files.end());

  return files;
}

// ---------------------------------------------------------------------------
// コンフィグファイルを読み込む
// ---------------------------------------------------------------------------
bool LoadConfig(const std::string& filepath, IntersectionConfig* config) {
  try {
    TomlConfigParser parser(filepath);

    // 軌道データ設定
    config->orbit_data_folder = parser.GetString("trajectory.orbit_data_folder", "");
    config->orbit_file_pattern = parser.GetString("trajectory.orbit_file_pattern", "*.txt");
    config->phase_step = parser.GetInt("trajectory.phase_step", 1);

    // 初期条件生成方式
    config->ic_source = parser.GetString("initial_conditions.source", "trajectory");
    config->trajectory_step = parser.GetInt("initial_conditions.trajectory_step", 10);

    // メッシュ生成設定
    config->mesh_type = parser.GetString("mesh.type", "concentric");
    config->r_min = parser.GetDouble("mesh.r_min", 0.001);
    config->r_max = parser.GetDouble("mesh.r_max", 0.02);
    config->n_radial = parser.GetInt("mesh.n_radial", 10);
    config->n_angular = parser.GetInt("mesh.n_angular", 36);
    config->grid_size = parser.GetInt("mesh.grid_size", 10);

    // ターゲットヤコビ積分設定
    config->target_jacobi = parser.GetDouble("jacobi.target_jacobi", 3.0);
    config->velocity_angle_min = parser.GetDouble("jacobi.velocity_angle_min", 0.0);
    config->velocity_angle_max = parser.GetDouble("jacobi.velocity_angle_max", 0.0);
    config->velocity_angle_step = parser.GetDouble("jacobi.velocity_angle_step", 10.0);

    // ポアンカレ断面設定
    config->section_var = parser.GetString("poincare_section.section_plane", "y");
    config->section_value = parser.GetDouble("poincare_section.section_value", 0.0);

    std::string crossing_str = parser.GetString("poincare_section.crossing_direction", "both");
    if (crossing_str == "positive") {
      config->crossing_direction = CrossingDirection::kPositive;
    } else if (crossing_str == "negative") {
      config->crossing_direction = CrossingDirection::kNegative;
    } else {
      config->crossing_direction = CrossingDirection::kBoth;
    }

    // Newton-Raphson設定
    config->newton_max_iterations = parser.GetInt("newton_raphson.max_iterations", 50);
    config->newton_tolerance = parser.GetDouble("newton_raphson.tolerance", 1e-12);
    config->max_integration_time = parser.GetDouble("newton_raphson.max_integration_time", 100.0);
    config->dt = parser.GetDouble("newton_raphson.dt", 0.0001);

    // 交点検証設定
    config->intersection_threshold = parser.GetDouble("intersection.threshold_distance", 0.001);
    config->verification_time = parser.GetDouble("intersection.verification_time", 500.0);

    // 積分器設定
    std::string integrator_str = parser.GetString("integrator.type", "symplectic6");
    if (integrator_str == "symplectic4") {
      config->integrator_type = IntegratorType::kSymplectic4th;
    } else if (integrator_str == "symplectic6") {
      config->integrator_type = IntegratorType::kSymplectic6th;
    } else if (integrator_str == "rk4") {
      config->integrator_type = IntegratorType::kRungeKutta4th;
    }

    // 出力設定
    config->output_tag = parser.GetString("output.tag", "");
    config->enable_phase2 = parser.GetBool("output.enable_phase2", true);
    config->enable_phase3 = parser.GetBool("output.enable_phase3", true);

    return true;
  } catch (const std::exception& e) {
    std::cerr << "Error loading config: " << e.what() << std::endl;
    return false;
  }
}

// ---------------------------------------------------------------------------
// 設定内容を表示
// ---------------------------------------------------------------------------
void PrintConfig(const IntersectionConfig& config) {
  std::cout << "<>    [Config Summary]" << std::endl;
  std::cout << "<>        orbit_data_folder: " << config.orbit_data_folder << std::endl;
  std::cout << "<>        orbit_file_pattern: " << config.orbit_file_pattern << std::endl;
  std::cout << "<>        phase_step: " << config.phase_step << std::endl;

  std::cout << "<>    [Initial Conditions]" << std::endl;
  std::cout << "<>        ic_source: " << config.ic_source << std::endl;
  if (config.ic_source == "trajectory") {
    std::cout << "<>        trajectory_step: " << config.trajectory_step << std::endl;
  } else {
    std::cout << "<>        mesh_type: " << config.mesh_type << std::endl;
    std::cout << "<>        r_min/r_max: " << config.r_min << " / " << config.r_max << std::endl;
    if (config.mesh_type == "concentric") {
      std::cout << "<>        n_radial x n_angular: " << config.n_radial << " x "
                << config.n_angular << std::endl;
    } else {
      std::cout << "<>        grid_size: " << config.grid_size << std::endl;
    }
  }

  std::cout << "<>    [Target Jacobi]" << std::endl;
  std::cout << "<>        target_jacobi: " << config.target_jacobi << std::endl;
  std::cout << "<>        velocity_angle: " << config.velocity_angle_min << " to "
            << config.velocity_angle_max << " (step=" << config.velocity_angle_step << " deg)"
            << std::endl;

  std::cout << "<>    [Phase Options]" << std::endl;
  std::cout << "<>        enable_phase2 (N-R refine): " << (config.enable_phase2 ? "true" : "false")
            << std::endl;
  std::cout << "<>        enable_phase3 (intersection): "
            << (config.enable_phase3 ? "true" : "false") << std::endl;
}

// ---------------------------------------------------------------------------
// 断面交差を検出する
// ---------------------------------------------------------------------------
bool DetectCrossing(double prev_val, double curr_val, double section_value,
                    CrossingDirection direction, CrossingDirection* detected_dir) {
  double prev_diff = prev_val - section_value;
  double curr_diff = curr_val - section_value;

  if (prev_diff * curr_diff >= 0) {
    return false;
  }

  // 検出された方向を記録
  if (prev_diff < 0 && curr_diff > 0) {
    *detected_dir = CrossingDirection::kPositive;
  } else {
    *detected_dir = CrossingDirection::kNegative;
  }

  switch (direction) {
    case CrossingDirection::kPositive:
      return prev_diff < 0 && curr_diff > 0;
    case CrossingDirection::kNegative:
      return prev_diff > 0 && curr_diff < 0;
    case CrossingDirection::kBoth:
      return true;
  }
  return false;
}

// ---------------------------------------------------------------------------
// 線形補間で交差点を求める
// ---------------------------------------------------------------------------
State<double> InterpolateCrossing(const State<double>& prev_state, const State<double>& curr_state,
                                  int section_index, double section_value, double* crossing_time,
                                  double prev_time, double curr_time) {
  double prev_val = GetStateValue(prev_state, section_index);
  double curr_val = GetStateValue(curr_state, section_index);

  double alpha = (section_value - prev_val) / (curr_val - prev_val);

  State<double> crossing_state;
  crossing_state.x = prev_state.x + alpha * (curr_state.x - prev_state.x);
  crossing_state.y = prev_state.y + alpha * (curr_state.y - prev_state.y);
  crossing_state.z = prev_state.z + alpha * (curr_state.z - prev_state.z);
  crossing_state.vx = prev_state.vx + alpha * (curr_state.vx - prev_state.vx);
  crossing_state.vy = prev_state.vy + alpha * (curr_state.vy - prev_state.vy);
  crossing_state.vz = prev_state.vz + alpha * (curr_state.vz - prev_state.vz);

  *crossing_time = prev_time + alpha * (curr_time - prev_time);

  return crossing_state;
}

// ---------------------------------------------------------------------------
// 2点間の距離を計算
// ---------------------------------------------------------------------------
double CalcDistance(const State<double>& s1, const State<double>& s2) {
  double dx = s1.x - s2.x;
  double dy = s1.y - s2.y;
  double dz = s1.z - s2.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// ---------------------------------------------------------------------------
// Phase 1: ポアンカレ断面交差点を抽出
// ---------------------------------------------------------------------------
std::vector<PoincareCrossing> ExtractPoincareCrossings(
    const std::vector<State<double>>& rotated_states, const IntersectionConfig& config, double mu) {
  std::vector<PoincareCrossing> crossings;

  int section_index = GetVarIndex(config.section_var);
  if (section_index < 0 || section_index > 2) {
    std::cerr << "<> !err! Invalid section_var: " << config.section_var << std::endl;
    return crossings;
  }

  int crossing_count = 0;
  for (size_t i = 1; i < rotated_states.size(); ++i) {
    const State<double>& prev_state = rotated_states[i - 1];
    const State<double>& curr_state = rotated_states[i];

    double prev_val = GetStateValue(prev_state, section_index);
    double curr_val = GetStateValue(curr_state, section_index);

    CrossingDirection detected_dir;
    if (DetectCrossing(prev_val, curr_val, config.section_value, config.crossing_direction,
                       &detected_dir)) {
      double crossing_time;
      State<double> crossing_state =
          InterpolateCrossing(prev_state, curr_state, section_index, config.section_value,
                              &crossing_time, static_cast<double>(i - 1), static_cast<double>(i));

      PoincareCrossing crossing;
      crossing.time = crossing_time;
      crossing.state = crossing_state;
      crossing.jacobi = calc_jacobi_integral(crossing_state, mu);
      crossing.crossing_index = ++crossing_count;
      crossing.detected_direction = detected_dir;

      crossings.push_back(crossing);
    }
  }

  return crossings;
}

// ---------------------------------------------------------------------------
// Phase 3: 参照軌道との交点検証
// ---------------------------------------------------------------------------
std::vector<IntersectionEvent> VerifyIntersection(
    const State<double>& initial_state, int candidate_idx,
    const std::vector<State<double>>& reference_trajectory, const IntersectionConfig& config,
    double mu) {
  std::vector<IntersectionEvent> events;

  State<double> state = initial_state;
  double time = 0.0;
  int num_steps = static_cast<int>(config.verification_time / config.dt);

  for (int step = 0; step < num_steps; ++step) {
    // 1ステップ積分
    State<double> next_state;
    switch (config.integrator_type) {
      case IntegratorType::kSymplectic4th:
        next_state = SymplecticStep4thOrder(mu, state, config.dt);
        break;
      case IntegratorType::kSymplectic6th:
        next_state = SymplecticStep6thOrder(mu, state, config.dt);
        break;
      default:
        next_state = SymplecticStep6thOrder(mu, state, config.dt);
        break;
    }

    // 参照軌道との最小距離を計算
    for (const auto& ref_point : reference_trajectory) {
      double dist = CalcDistance(next_state, ref_point);
      if (dist < config.intersection_threshold) {
        IntersectionEvent event;
        event.candidate_idx = candidate_idx;
        event.time = time + config.dt;
        event.candidate_state = next_state;
        event.reference_state = ref_point;
        event.distance = dist;
        events.push_back(event);
      }
    }

    state = next_state;
    time += config.dt;

    // 脱出判定（L2点の外側に出たら終了）
    double r = std::sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
    if (r > 3.0) break;
  }

  return events;
}

// ---------------------------------------------------------------------------
// メイン関数
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>  Asteroid Orbit Intersection Finder v1.0" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // コマンドライン引数のパース
  CommonArgs args = ParseCommonArgs(argc, argv);
  bool is_continuous = args.is_continuous;
  bool skip_wait = args.skip_wait;
  std::string cli_output_tag = args.output_tag;

  // 天文パラメータの読み込み
  const std::string kConfigFilePath = CONFIG_DIR;
  std::string astro_param_file = kConfigFilePath + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(astro_param_file);
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  std::cout << "<>    mu = " << std::setprecision(10) << kMU << std::endl;

  // コンフィグファイルのリストを取得
  const std::string kCalcConfigPath = kConfigFilePath + "/asteroid_orbit_intersection/";
  std::string config_prefix = "intersection_config";

  ConfigDiscoveryOptions discovery_opts;
  discovery_opts.exclude_sample = true;
  discovery_opts.continuous_mode = is_continuous;
  std::vector<std::string> config_file_list =
      DiscoverConfigFilesToml(kCalcConfigPath, config_prefix, discovery_opts);

  if (config_file_list.empty()) {
    std::cerr << "<> No config files found in " << kCalcConfigPath << std::endl;
    return -1;
  }

  std::cout << "<> Loaded config file list:" << std::endl;
  for (const auto& filepath : config_file_list) {
    std::cout << "<>    - " << filepath << std::endl;
  }

  if (!skip_wait) {
    WaitForEnter();
  }

  // OpenMP スレッド数の設定
  int Core_Max = omp_get_max_threads();
  int OMP_Fmax = skip_wait ? (Core_Max - 1) : Core_Max;
  if (!skip_wait) {
    std::cout << "<>  [OpenMP preparation]" << std::endl;
    std::cout << "<>  Threads available: " << Core_Max << std::endl;
    std::cout << "<>  Enter number of threads to use: ";
    std::cin >> OMP_Fmax;
    if (OMP_Fmax <= 0 || OMP_Fmax > Core_Max) {
      OMP_Fmax = Core_Max - 1;
    }
    WaitForEnter();
  }
  omp_set_num_threads(OMP_Fmax);
  std::cout << "<>    Using " << OMP_Fmax << " OpenMP threads" << std::endl;

  // セッション出力ディレクトリを作成
  OutputDirResult output_result =
      CreateSessionOutputDir(OUTPUT_DIR, "asteroid_orbit_intersection", cli_output_tag);
  if (!output_result.success) {
    return -1;
  }
  std::string session_output_dir = output_result.session_dir;
  std::cout << "<>    Session output directory: " << session_output_dir << std::endl;

  auto start_ofall = std::chrono::system_clock::now();

  // 各コンフィグファイルを処理
  for (const auto& configfilepath : config_file_list) {
    std::cout << "\n<>----------------------------------------------------------------"
              << std::endl;
    std::cout << "<>    Processing: " << configfilepath << std::endl;

    IntersectionConfig config;
    if (!LoadConfig(configfilepath, &config)) {
      std::cerr << "<>    Failed to load config, skipping..." << std::endl;
      continue;
    }

    PrintConfig(config);

    if (!is_continuous && !skip_wait) {
      WaitForEnter();
    }

    // 出力ディレクトリ作成
    std::string config_basename = fs::path(configfilepath).stem().string();
    std::string output_dir = session_output_dir + "/" + config_basename;
    if (!config.output_tag.empty()) {
      output_dir += "_" + config.output_tag;
    }
    fs::create_directories(output_dir);

    // コンフィグファイルをコピー
    std::string config_copy_path = output_dir + "/" + fs::path(configfilepath).filename().string();
    try {
      fs::copy_file(configfilepath, config_copy_path, fs::copy_options::overwrite_existing);
    } catch (const fs::filesystem_error& e) {
      std::cerr << "<>    Warning: Failed to copy config: " << e.what() << std::endl;
    }

    auto start = std::chrono::system_clock::now();

    // ========================================
    // フォルダ内の軌道データファイルを探索
    // ========================================
    std::cout << "<>    Discovering orbit data files..." << std::endl;
    std::vector<std::string> orbit_files =
        DiscoverOrbitDataFiles(config.orbit_data_folder, config.orbit_file_pattern);

    if (orbit_files.empty()) {
      std::cerr << "<>    No orbit files found in: " << config.orbit_data_folder << std::endl;
      std::cerr << "<>    Pattern: " << config.orbit_file_pattern << std::endl;
      continue;
    }

    std::cout << "<>        Found " << orbit_files.size() << " orbit files" << std::endl;
    for (const auto& f : orbit_files) {
      std::cout << "<>            - " << fs::path(f).filename().string() << std::endl;
    }

    // 各軌道ファイルに対して処理を実行
    for (const auto& orbit_file_path : orbit_files) {
      std::string orbit_basename = fs::path(orbit_file_path).stem().string();
      std::cout << "\n<>    ----------------------------------------" << std::endl;
      std::cout << "<>    Processing orbit file: " << orbit_basename << std::endl;

      // 軌道ファイルごとの出力ディレクトリ
      std::string file_output_dir = output_dir + "/" + orbit_basename;
      fs::create_directories(file_output_dir);

      // ========================================
      // Phase 0: 軌道データ読み込み・軌道面法線計算
      // ========================================
      std::cout << "<>    [Phase 0] Loading orbit data & calculating orbital plane..." << std::endl;

      std::vector<OrbitDataRow> orbit_data;
      if (!LoadOrbitData(orbit_file_path, config.phase_step, &orbit_data)) {
        std::cerr << "<>    Failed to load orbit data, skipping..." << std::endl;
        continue;
      }
      std::cout << "<>        Loaded " << orbit_data.size() << " data points" << std::endl;

      // J2000 -> 回転座標系に変換
      std::vector<State<double>> rotated_states;
      rotated_states.reserve(orbit_data.size());
      for (const auto& row : orbit_data) {
        State<double> rotated = ConvertInertial2RotatingV2(row.asteroid, row.earth, astro_params);
        rotated_states.push_back(rotated);
      }
      std::cout << "<>        Coordinate conversion completed" << std::endl;

      // 軌道の重心（メッシュ中心）を計算
      State3d<double> mesh_center{0.0, 0.0, 0.0};
      for (const auto& s : rotated_states) {
        mesh_center.x += s.x;
        mesh_center.y += s.y;
        mesh_center.z += s.z;
      }
      mesh_center.x /= rotated_states.size();
      mesh_center.y /= rotated_states.size();
      mesh_center.z /= rotated_states.size();
      std::cout << "<>        Mesh center: (" << mesh_center.x << ", " << mesh_center.y << ", "
                << mesh_center.z << ")" << std::endl;

      // 軌道面法線ベクトルを計算
      State3d<double> plane_normal = CalcOrbitalPlaneNormal(rotated_states, mesh_center);
      std::cout << "<>        Orbital plane normal: (" << plane_normal.x << ", " << plane_normal.y
                << ", " << plane_normal.z << ")" << std::endl;

      // 回転座標系軌道を出力
      std::string rotated_path = file_output_dir + "/asteroid_rotated.csv";
      std::ofstream ofs_rotated(rotated_path);
      ofs_rotated << "# Asteroid trajectory in rotating frame\n";
      ofs_rotated << "idx,x,y,z,vx,vy,vz,jacobi\n";
      for (size_t i = 0; i < rotated_states.size(); ++i) {
        const auto& s = rotated_states[i];
        double jacobi = calc_jacobi_integral(s, kMU);
        ofs_rotated << std::setprecision(12) << i << "," << s.x << "," << s.y << "," << s.z << ","
                    << s.vx << "," << s.vy << "," << s.vz << "," << jacobi << "\n";
      }
      ofs_rotated.close();

      // ========================================
      // Phase 1: 初期条件生成
      // ========================================
      std::cout << "<>    [Phase 1] Generating initial conditions..." << std::endl;
      std::cout << "<>        Source: " << config.ic_source << std::endl;
      std::cout << "<>        Target Jacobi: " << config.target_jacobi << std::endl;

      std::vector<State<double>> initial_conditions;
      int forbidden_count = 0;

      if (config.ic_source == "trajectory") {
        // ========================================
        // 軌道点から直接初期条件を生成（位置は元軌道と完全一致）
        // 各軌道点で複数の速度方向を生成
        // ========================================
        int n_angles = 1;
        if (config.velocity_angle_max > config.velocity_angle_min &&
            config.velocity_angle_step > 0) {
          n_angles = static_cast<int>((config.velocity_angle_max - config.velocity_angle_min) /
                                      config.velocity_angle_step) +
                     1;
        }
        std::cout << "<>        Using trajectory points (step=" << config.trajectory_step
                  << ", angles=" << n_angles << ")" << std::endl;

        for (size_t i = 0; i < rotated_states.size(); i += config.trajectory_step) {
          const auto& orig = rotated_states[i];

          // 位置は元軌道から直接使用
          double x = orig.x;
          double y = orig.y;
          double z = orig.z;

          // ヤコビ積分から速度の大きさを計算
          double v_abs = CalcVelocityFromJacobi(config.target_jacobi, x, y, z, kMU);
          if (v_abs < 0) {
            ++forbidden_count;
            continue;  // 禁止領域
          }

          // 元軌道の速度方向を基準ベクトルとして取得
          double orig_v_abs = std::sqrt(orig.vx * orig.vx + orig.vy * orig.vy + orig.vz * orig.vz);
          if (orig_v_abs < 1e-12) {
            // 速度がほぼゼロの場合はスキップ
            continue;
          }
          State3d<double> v_dir{orig.vx / orig_v_abs, orig.vy / orig_v_abs, orig.vz / orig_v_abs};

          // 角運動量方向（速度回転の軸）を計算: h = r × v / |r × v|
          State3d<double> h{y * orig.vz - z * orig.vy, z * orig.vx - x * orig.vz,
                            x * orig.vy - y * orig.vx};
          double h_norm = std::sqrt(h.x * h.x + h.y * h.y + h.z * h.z);
          if (h_norm > 1e-12) {
            h.x /= h_norm;
            h.y /= h_norm;
            h.z /= h_norm;
          }

          // 各速度角度で初期条件を生成
          for (double angle = config.velocity_angle_min; angle <= config.velocity_angle_max + 1e-9;
               angle += (config.velocity_angle_step > 0 ? config.velocity_angle_step : 360.0)) {
            // 角運動量軸周りで速度方向を回転
            double theta = angle * M_PI / 180.0;
            double cos_t = std::cos(theta);
            double sin_t = std::sin(theta);
            // Rodrigues回転公式: v_rot = v*cos(θ) + (h×v)*sin(θ) + h*(h·v)*(1-cos(θ))
            double h_dot_v = h.x * v_dir.x + h.y * v_dir.y + h.z * v_dir.z;
            State3d<double> h_cross_v{h.y * v_dir.z - h.z * v_dir.y, h.z * v_dir.x - h.x * v_dir.z,
                                      h.x * v_dir.y - h.y * v_dir.x};
            State3d<double> v_rot{
                v_dir.x * cos_t + h_cross_v.x * sin_t + h.x * h_dot_v * (1 - cos_t),
                v_dir.y * cos_t + h_cross_v.y * sin_t + h.y * h_dot_v * (1 - cos_t),
                v_dir.z * cos_t + h_cross_v.z * sin_t + h.z * h_dot_v * (1 - cos_t)};

            State<double> state;
            state.x = x;
            state.y = y;
            state.z = z;
            state.vx = v_rot.x * v_abs;
            state.vy = v_rot.y * v_abs;
            state.vz = v_rot.z * v_abs;
            initial_conditions.push_back(state);

            if (config.velocity_angle_step <= 0) break;  // 単一角度の場合
          }
        }
      } else {
        // ========================================
        // メッシュ点から初期条件を生成（従来方式）
        // ========================================
        std::vector<State3d<double>> mesh_points;
        if (config.mesh_type == "concentric") {
          mesh_points = CreateConcentricMesh(mesh_center, plane_normal, config.r_min, config.r_max,
                                             config.n_radial, config.n_angular);
        } else {
          mesh_points =
              CreateRectangularMesh(mesh_center, plane_normal, config.r_max, config.grid_size);
        }
        std::cout << "<>        Generated " << mesh_points.size() << " mesh points" << std::endl;

        for (const auto& pos : mesh_points) {
          double v_abs = CalcVelocityFromJacobi(config.target_jacobi, pos.x, pos.y, pos.z, kMU);
          if (v_abs < 0) {
            ++forbidden_count;
            continue;
          }

          // 各速度角度で初期条件を生成
          for (double angle = config.velocity_angle_min; angle <= config.velocity_angle_max + 1e-9;
               angle += (config.velocity_angle_step > 0 ? config.velocity_angle_step : 360.0)) {
            State3d<double> velocity =
                CalcVelocityOnPlane(pos, mesh_center, plane_normal, v_abs, angle);

            State<double> state;
            state.x = pos.x;
            state.y = pos.y;
            state.z = pos.z;
            state.vx = velocity.x;
            state.vy = velocity.y;
            state.vz = velocity.z;
            initial_conditions.push_back(state);

            if (config.velocity_angle_step <= 0) break;
          }
        }
      }

      std::cout << "<>        Valid initial conditions: " << initial_conditions.size() << std::endl;
      if (forbidden_count > 0) {
        std::cout << "<>        Skipped " << forbidden_count << " points in forbidden region"
                  << std::endl;
      }

      // 初期条件をCSVに出力
      std::string ic_path = file_output_dir + "/initial_conditions.csv";
      std::ofstream ofs_ic(ic_path);
      ofs_ic << "# Initial conditions (source=" << config.ic_source
             << ", target_jacobi=" << config.target_jacobi << ")\n";
      ofs_ic << "idx,x,y,z,vx,vy,vz,jacobi\n";
      for (size_t i = 0; i < initial_conditions.size(); ++i) {
        const auto& s = initial_conditions[i];
        double jacobi = calc_jacobi_integral(s, kMU);
        ofs_ic << std::setprecision(12) << i << "," << s.x << "," << s.y << "," << s.z << ","
               << s.vx << "," << s.vy << "," << s.vz << "," << jacobi << "\n";
      }
      ofs_ic.close();
      std::cout << "<>        Initial conditions saved to: " << ic_path << std::endl;

      // ========================================
      // Phase 2: 周期軌道探索（Newton-Raphson） [オプション]
      // ========================================
      std::vector<PeriodicOrbit<double>> refined_orbits;
      if (config.enable_phase2 && !initial_conditions.empty()) {
        std::cout << "<>    [Phase 2] Refining periodic orbits (N-R)..." << std::endl;

        int section_index = GetVarIndex(config.section_var);
        int converged_count = 0;

        for (size_t i = 0; i < initial_conditions.size(); ++i) {
          const auto& ic = initial_conditions[i];

          NewtonConvergenceInfo<double> convergence_info;
          PeriodicOrbit<double> orbit = RefinePeriodicOrbit(
              ic, kMU, section_index, config.section_value, config.newton_max_iterations,
              config.newton_tolerance, config.max_integration_time, config.dt, &convergence_info);

          if (convergence_info.converged) {
            refined_orbits.push_back(orbit);
            ++converged_count;
          }
        }

        std::cout << "<>        Converged: " << converged_count << " / "
                  << initial_conditions.size() << std::endl;

        // 周期軌道を出力
        std::string periodic_path = file_output_dir + "/periodic_orbits.csv";
        std::ofstream ofs_periodic(periodic_path);
        ofs_periodic << "# Refined periodic orbits\n";
        ofs_periodic << "idx,x0,y0,z0,vx0,vy0,vz0,period,jacobi\n";
        for (size_t i = 0; i < refined_orbits.size(); ++i) {
          const auto& orbit = refined_orbits[i];
          double jacobi = calc_jacobi_integral(orbit.initial_state, kMU);
          ofs_periodic << std::setprecision(12) << i << "," << orbit.initial_state.x << ","
                       << orbit.initial_state.y << "," << orbit.initial_state.z << ","
                       << orbit.initial_state.vx << "," << orbit.initial_state.vy << ","
                       << orbit.initial_state.vz << "," << orbit.period << "," << jacobi << "\n";
        }
        ofs_periodic.close();
        std::cout << "<>        Periodic orbits saved to: " << periodic_path << std::endl;
      } else if (!config.enable_phase2) {
        std::cout << "<>    [Phase 2] Skipped (N-R disabled)" << std::endl;
      }

      // ========================================
      // Phase 3: 元軌道との交点検証
      // ========================================
      std::vector<IntersectionEvent> all_intersections;
      if (config.enable_phase3) {
        std::cout << "<>    [Phase 3] Verifying intersections with reference trajectory..."
                  << std::endl;

        // 候補軌道のリストを構築（Phase2実行時は周期軌道、未実行時は初期条件を使用）
        std::vector<std::pair<int, State<double>>> candidates;
        if (!refined_orbits.empty()) {
          for (size_t i = 0; i < refined_orbits.size(); ++i) {
            candidates.push_back({static_cast<int>(i), refined_orbits[i].initial_state});
          }
          std::cout << "<>        Using " << candidates.size() << " refined periodic orbits"
                    << std::endl;
        } else {
          for (size_t i = 0; i < initial_conditions.size(); ++i) {
            candidates.push_back({static_cast<int>(i), initial_conditions[i]});
          }
          std::cout << "<>        Using " << candidates.size()
                    << " initial conditions (quasi-periodic search)" << std::endl;
        }

        int completed = 0;
        int total = static_cast<int>(candidates.size());

#pragma omp parallel for schedule(dynamic) shared(all_intersections, completed)
        for (int idx = 0; idx < total; ++idx) {
          const auto& [cand_idx, initial_state] = candidates[idx];

          std::vector<IntersectionEvent> events =
              VerifyIntersection(initial_state, cand_idx, rotated_states, config, kMU);

#pragma omp critical
          {
            all_intersections.insert(all_intersections.end(), events.begin(), events.end());
            completed++;
            if (completed % 10 == 0 || completed == total) {
              std::cout << "<>        Progress: " << completed << "/" << total << std::endl;
            }
          }
        }

        std::cout << "<>        Found " << all_intersections.size() << " intersection events"
                  << std::endl;

        // 交点イベントを出力
        std::string intersections_path = file_output_dir + "/intersection_events.csv";
        std::ofstream ofs_intersections(intersections_path);
        ofs_intersections << "# Intersection events with reference trajectory\n";
        ofs_intersections << "candidate_idx,time,cand_x,cand_y,cand_z,ref_x,ref_y,ref_z,distance\n";
        for (const auto& e : all_intersections) {
          ofs_intersections << std::setprecision(12) << e.candidate_idx << "," << e.time << ","
                            << e.candidate_state.x << "," << e.candidate_state.y << ","
                            << e.candidate_state.z << "," << e.reference_state.x << ","
                            << e.reference_state.y << "," << e.reference_state.z << ","
                            << e.distance << "\n";
        }
        ofs_intersections.close();
      }

      auto end = std::chrono::system_clock::now();
      double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
      std::cout << "<>    File processing completed in " << elapsed << " seconds" << std::endl;
      std::cout << "<>    Output directory: " << file_output_dir << std::endl;
    }  // end of orbit file loop
  }  // end of config file loop

  auto end_ofall = std::chrono::system_clock::now();
  double total_elapsed =
      std::chrono::duration_cast<std::chrono::seconds>(end_ofall - start_ofall).count();
  std::cout << "\n<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    All processing completed in " << total_elapsed << " seconds" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  return 0;
}
