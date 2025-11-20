#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "rtbp.hpp"
#include "utils.hpp"
#include "vector3d.hpp"

using ::kSecondsPerDay;
using ::Vector3d;
using crtbp::AstroConstants;
using crtbp::calc_r2;
using crtbp::CanonicalState;
using crtbp::ConvertInertial2RotatingV2;
using crtbp::ConvertToCanonical;
using crtbp::ConvertToPhysical;
using crtbp::SaliState;
using crtbp::State;
using utils::loadConstants;

namespace fs = std::filesystem;
constexpr double kPi = std::numbers::pi;
constexpr double kFilteredScore = -1.0e9;

/**
 * @brief スティッキー軌道探索の設定
 */
struct StickyConfig {
  State<double> asteroid_state_helio{};
  State<double> earth_state_helio{};

  // 自然軌道の積分設定 (days)
  double natural_duration_days = 5.0;
  double natural_dt_days = 1e-3;
  int sample_stride = 1;

  // Delta-v マップ設定
  double dv_magnitude_mps = 100.0;
  double theta_start_deg = -180.0;
  double theta_end_deg = 180.0;
  double theta_step_deg = 5.0;
  double phi_start_deg = 0.0;
  double phi_end_deg = 0.0;
  double phi_step_deg = 90.0;

  // 生存時間計測の設定
  double survival_max_days = 36500.0;  // 最大生存時間 (100年など)
  double survival_dt_days = 1e-3;
  double escape_radius_hill = 3.0;      // 離脱判定距離 (ヒル半径の倍数)
  double collision_radius_km = 6500.0;  // 衝突判定距離 (地球半径+大気)
  double sali_threshold = 1e-8;         // SALI threshold used to detect sticky decay
  double jacobi_min = -std::numeric_limits<double>::infinity();
  double jacobi_max = std::numeric_limits<double>::infinity();  // Jacobi gate window after DV

  std::string output_tag;
};

/**
 * @brief 自然単位系
 */
struct NaturalUnits {
  double mu = 0.0;
  double lu_au = 0.0;
  double tu_day = 0.0;
  double vu_au_per_day = 0.0;
  double vu_m_per_s = 0.0;
  double hill_radius = 0.0;
  double collision_radius = 0.0;  // 無次元化された衝突半径
};

/**
 * @brief 軌道サンプル
 */
struct TrajectorySample {
  double t_nd = 0.0;
  double t_days = 0.0;
  State<double> state{};
  double r2 = 0.0;
};

/**
 * @brief 探索タスク
 */
struct SearchTask {
  const TrajectorySample* sample = nullptr;
  double theta_deg = 0.0;
  double phi_deg = 0.0;
};

/**
 * @brief 終了理由
 */
enum class EndReason {
  TIME_UP,    // 指定時間まで生存
  COLLISION,  // 地球に衝突
  ESCAPE,     // ヒル球外へ離脱
  FILTERED,   // プレフィルタで棄却
  ERROR       // 数値計算エラー等
};

/**
 * @brief 探索結果
 */
struct SearchResult {
  double t_nd = 0.0;
  double t_days = 0.0;
  double theta_deg = 0.0;
  double survival_days = 0.0;
  double sali_decay_days = 0.0;  // SALIが閾値を下回るまでの時間
  double final_sali = 0.0;
  double score = 0.0;
  double jacobi_constant = 0.0;
  EndReason reason = EndReason::ERROR;
  State<double> final_state{};
  Vector3d<double> dv_direction{};
};

// Forward declarations
bool LoadStickyConfig(const std::string& filepath, StickyConfig* config);
NaturalUnits ComputeNaturalUnits(const State<double>& earth_state_helio,
                                 const AstroConstants<double>& astro, double collision_km);
std::vector<TrajectorySample> IntegrateNaturalTrajectory(const AstroConstants<double>& astro,
                                                         const State<double>& initial_state_nd,
                                                         double duration_nd, double dt_nd,
                                                         double mu, const NaturalUnits& units);
bool EvaluateStickyOrbit(const SearchTask& task, double theta_rad, double phi_rad, double dv_nd,
                         double mu, double escape_r, double collision_r, double jacobi_min,
                         double jacobi_max, double max_duration_nd, double dt_nd,
                         double sali_threshold, SearchResult* output);
void WriteResultCsv(const std::string& filepath, const NaturalUnits& units, const StickyConfig& cfg,
                    const std::vector<SearchTask>& tasks, double dv_nd);

// Level 1: Perigee Filter
bool CheckPerigeeFilter(const State<double>& state_rot, double mu, double collision_r);

std::string EndReasonToString(EndReason r) {
  switch (r) {
    case EndReason::TIME_UP:
      return "TIME_UP";
    case EndReason::COLLISION:
      return "COLLISION";
    case EndReason::ESCAPE:
      return "ESCAPE";
    case EndReason::FILTERED:
      return "FILTERED";
    default:
      return "ERROR";
  }
}

int main() {
  using namespace param;

  const std::string kConfigDir = std::string(CONFIG_DIR) + "/sticky_orbit_search/";
  const std::string kOutputDir = std::string(OUTPUT_DIR) + "/sticky_orbit_search";
  const std::string kConfigPrefix = "sticky_search";

  AstroConstants<double> astro_params =
      loadConstants<double>(kConfigDir + "/../astro_param/astro_param.txt");
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  std::cout << "<>------------------------------------------------------------\n";
  std::cout << "<>  Sticky Orbit Search (Survival Time & SALI Map)\n";
  std::cout << "<>------------------------------------------------------------\n\n";

  std::vector<std::string> config_files = utils::DiscoverConfigFiles(kConfigDir, kConfigPrefix);
  if (config_files.empty()) {
    // Fallback: try finding any .txt if prefix match fails (or just rely on DiscoverConfigFiles)
    // DiscoverConfigFiles implementation in utils.hpp handles regex matching.
    // If empty, maybe user didn't name it correctly. Let's try to list all .txt for now if empty.
    if (fs::exists(kConfigDir)) {
      for (const auto& entry : fs::directory_iterator(kConfigDir)) {
        if (entry.path().extension() == ".txt") {
          // Check if already in list
          bool found = false;
          for (const auto& f : config_files)
            if (f == fs::absolute(entry.path()).string()) found = true;
          if (!found) config_files.push_back(fs::absolute(entry.path()).string());
        }
      }
    }
  }

  if (config_files.empty()) {
    std::cerr << "No config files found in " << kConfigDir << "\n";
    return 1;
  }

  fs::create_directories(kOutputDir);

  int config_index = 0;
  for (const auto& config_path : config_files) {
    ++config_index;
    StickyConfig cfg;
    if (!LoadStickyConfig(config_path, &cfg)) {
      continue;
    }
    if (cfg.jacobi_min > cfg.jacobi_max) std::swap(cfg.jacobi_min, cfg.jacobi_max);

    NaturalUnits units{};
    try {
      units = ComputeNaturalUnits(cfg.earth_state_helio, astro_params, cfg.collision_radius_km);
    } catch (const std::exception& e) {
      std::cerr << "Unit computation failed: " << e.what() << "\n";
      continue;
    }

    std::cout << "<> Config #" << config_index << ": " << fs::path(config_path).filename() << "\n";
    std::cout << "   Collision Radius: " << cfg.collision_radius_km << " km ("
              << units.collision_radius << " nd)\n";
    std::cout << "   Escape Radius: " << cfg.escape_radius_hill << " * Hill ("
              << units.hill_radius * cfg.escape_radius_hill << " nd)\n";

    State<double> initial_state_nd{};
    try {
      initial_state_nd =
          ConvertInertial2RotatingV2(cfg.asteroid_state_helio, cfg.earth_state_helio, astro_params);
    } catch (const std::exception& e) {
      std::cerr << "Coordinate conversion failed: " << e.what() << "\n";
      continue;
    }

    // 自然軌道積分
    const double natural_dur_nd = cfg.natural_duration_days / units.tu_day;
    const double natural_dt_nd = cfg.natural_dt_days / units.tu_day;
    auto samples = IntegrateNaturalTrajectory(astro_params, initial_state_nd, natural_dur_nd,
                                              natural_dt_nd, kMU, units);

    if (samples.empty()) {
      std::cerr << "No samples generated.\n";
      continue;
    }

    // タスク生成
    std::vector<SearchTask> tasks;
    const size_t stride = std::max(1, cfg.sample_stride);
    auto compute_step = [](double start, double end, double raw_step) {
      double step = (raw_step == 0.0) ? 1.0 : std::abs(raw_step);
      if (end < start) step = -step;
      return step;
    };

    const double theta_step = compute_step(cfg.theta_start_deg, cfg.theta_end_deg, cfg.theta_step_deg);
    const double phi_step = compute_step(cfg.phi_start_deg, cfg.phi_end_deg, cfg.phi_step_deg);

    for (size_t i = 0; i < samples.size(); i += stride) {
      for (double th = cfg.theta_start_deg;
           (theta_step >= 0.0 ? (th <= cfg.theta_end_deg + 1e-9) : (th >= cfg.theta_end_deg - 1e-9));
           th += theta_step) {
        for (double ph = cfg.phi_start_deg;
             (phi_step >= 0.0 ? (ph <= cfg.phi_end_deg + 1e-9) : (ph >= cfg.phi_end_deg - 1e-9));
             ph += phi_step) {
          tasks.push_back(SearchTask{&samples[i], th, ph});
        }
      }
    }

    std::cout << "   Tasks: " << tasks.size() << "\n";

    const double dv_nd = cfg.dv_magnitude_mps / units.vu_m_per_s;
    WriteResultCsv(kOutputDir + "/" + (cfg.output_tag.empty() ? "result" : cfg.output_tag) + ".csv",
                   units, cfg, tasks, dv_nd);
  }

  return 0;
}

bool LoadStickyConfig(const std::string& filepath, StickyConfig* config) {
  std::ifstream ifs(filepath);
  if (!ifs) return false;

  auto trim = [](const std::string& s) {
    auto start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return std::string();
    auto end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
  };

  auto parse_state = [](const std::string& s, State<double>* out) {
    std::stringstream ss(s);
    double v[6];
    for (int i = 0; i < 6; ++i) ss >> v[i];
    *out = State<double>{v[0], v[1], v[2], v[3], v[4], v[5]};
  };

  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    auto pos = line.find('=');
    if (pos == std::string::npos) continue;
    std::string key = trim(line.substr(0, pos));
    std::string val = trim(line.substr(pos + 1));

    if (key == "ASTEROID_STATE_HELIO_AU")
      parse_state(val, &config->asteroid_state_helio);
    else if (key == "EARTH_STATE_HELIO_AU")
      parse_state(val, &config->earth_state_helio);
    else if (key == "NATURAL_DURATION_DAYS")
      config->natural_duration_days = std::stod(val);
    else if (key == "NATURAL_DT_DAYS")
      config->natural_dt_days = std::stod(val);
    else if (key == "SAMPLE_STRIDE")
      config->sample_stride = std::stoi(val);
    else if (key == "DV_MPS")
      config->dv_magnitude_mps = std::stod(val);
    else if (key == "THETA_RANGE_DEG") {
      std::stringstream ss(val);
      ss >> config->theta_start_deg >> config->theta_end_deg >> config->theta_step_deg;
    } else if (key == "PHI_RANGE_DEG") {
      std::stringstream ss(val);
      ss >> config->phi_start_deg >> config->phi_end_deg >> config->phi_step_deg;
    } else if (key == "SURVIVAL_MAX_DAYS")
      config->survival_max_days = std::stod(val);
    else if (key == "SURVIVAL_DT_DAYS")
      config->survival_dt_days = std::stod(val);
    else if (key == "ESCAPE_RADIUS_HILL")
      config->escape_radius_hill = std::stod(val);
    else if (key == "COLLISION_RADIUS_KM")
      config->collision_radius_km = std::stod(val);
    else if (key == "SALI_THRESHOLD")
      config->sali_threshold = std::stod(val);
    else if (key == "CJ_GATE_RANGE" || key == "JACOBI_RANGE") {
      std::stringstream ss(val);
      ss >> config->jacobi_min >> config->jacobi_max;
    }
    else if (key == "OUTPUT_TAG")
      config->output_tag = val;
  }
  return true;
}

NaturalUnits ComputeNaturalUnits(const State<double>& earth_state_helio,
                                 const AstroConstants<double>& astro, double collision_km) {
  NaturalUnits units;
  units.mu = astro.gm_earth / (astro.gm_earth + astro.gm_sun);

  // AU, Day -> ND conversion
  double kM3s2_to_AU3day2 = (kSecondsPerDay * kSecondsPerDay) / (astro.au * astro.au * astro.au);
  double gm_total_AD = (astro.gm_sun + astro.gm_earth) * kM3s2_to_AU3day2;

  Vector3d<double> r_vec{earth_state_helio.x, earth_state_helio.y, earth_state_helio.z};
  double r_mag = r_vec.magnitude();
  double n = std::sqrt(gm_total_AD / std::pow(r_mag, 3));

  units.lu_au = r_mag;
  units.tu_day = 1.0 / n;
  units.vu_au_per_day = units.lu_au / units.tu_day;
  units.vu_m_per_s = units.vu_au_per_day * astro.au / kSecondsPerDay;
  units.hill_radius = std::cbrt(units.mu / 3.0);

  // 衝突半径の無次元化 (km -> AU -> ND)
  // astro.au is usually meters.
  units.collision_radius = (collision_km * 1000.0 / astro.au) / units.lu_au;

  return units;
}

std::vector<TrajectorySample> IntegrateNaturalTrajectory(const AstroConstants<double>& astro,
                                                         const State<double>& initial_state_nd,
                                                         double duration_nd, double dt_nd,
                                                         double mu, const NaturalUnits& units) {
  std::vector<TrajectorySample> samples;
  int steps = static_cast<int>(std::ceil(duration_nd / dt_nd));
  samples.reserve(steps + 1);

  auto observer = [&](const State<double>& s, double t) {
    double r2 = calc_r2(s.x, s.y, s.z, mu);
    samples.push_back({t, t * units.tu_day, s, r2});
  };

  auto integrator = [&](const State<double>& s, double t, double dt) {
    return crtbp::SymplecticStep6thOrder(mu, s, dt);
  };

  State<double> state = initial_state_nd;
  crtbp::Integrate(state, integrator, observer, 0.0, dt_nd, steps);
  return samples;
}

bool CheckPerigeeFilter(const State<double>& state_rot, double mu, double collision_r) {
  // 回転座標系での位置・速度
  Vector3d<double> r{state_rot.x - (1.0 - mu), state_rot.y, state_rot.z};  // 地球中心
  Vector3d<double> v_rot{state_rot.vx, state_rot.vy, state_rot.vz};

  // 慣性系への速度変換 (簡易版: 地球中心慣性系とみなす)
  // v_in = v_rot + omega x r
  // omega = (0, 0, 1)
  Vector3d<double> omega{0.0, 0.0, 1.0};
  Vector3d<double> v_in = v_rot + omega.gaiseki(r);

  double r_mag = r.magnitude();
  double v_mag2 = v_in.magnitude() * v_in.magnitude();

  double specific_energy = 0.5 * v_mag2 - mu / r_mag;

  if (specific_energy >= 0) {
    // 双曲線軌道 -> 衝突チェックは近地点距離で行う
    // a < 0 となる。
    double a = -mu / (2.0 * specific_energy);
    // 角運動量 h = r x v
    Vector3d<double> h_vec = r.gaiseki(v_in);
    double h2 = h_vec.magnitude() * h_vec.magnitude();
    // e = sqrt(1 + 2 E h^2 / mu^2)
    double e = std::sqrt(1.0 + 2.0 * specific_energy * h2 / (mu * mu));
    double r_p = (h2 / mu) / (1.0 + e);

    if (r_p < collision_r) return false;  // 衝突
    return true;                          // 衝突はしないが離脱する
  }

  double a = -mu / (2.0 * specific_energy);

  // 角運動量
  Vector3d<double> h_vec = r.gaiseki(v_in);
  double h2 = h_vec.magnitude() * h_vec.magnitude();

  // 離心率 e = sqrt(1 + 2Eh^2/mu^2)
  double term = 1.0 + 2.0 * specific_energy * h2 / (mu * mu);
  double e = (term > 0) ? std::sqrt(term) : 0.0;

  double r_p = a * (1.0 - e);

  if (r_p < collision_r) {
    return false;  // 衝突コース
  }
  return true;
}

bool EvaluateStickyOrbit(const SearchTask& task, double theta_rad, double phi_rad, double dv_nd,
                         double mu, double escape_r, double collision_r, double jacobi_min,
                         double jacobi_max, double max_duration_nd, double dt_nd,
                         double sali_threshold, SearchResult* output) {
  if (!output || !task.sample) return false;

  output->t_nd = task.sample->t_nd;
  output->t_days = task.sample->t_days;
  output->theta_deg = task.theta_deg;
  output->phi_deg = task.phi_deg;

  const auto& base = task.sample->state;
  Vector3d<double> v{base.vx, base.vy, base.vz};
  Vector3d<double> r_rel{base.x - (1.0 - mu), base.y, base.z};

  Vector3d<double> tangent = v;
  if (tangent.magnitude() < 1e-12) {
    tangent = Vector3d<double>(1.0, 0.0, 0.0);
  } else {
    tangent = tangent.normalise();
  }

  Vector3d<double> h = r_rel.gaiseki(v);
  if (h.magnitude() < 1e-12) {
    h = Vector3d<double>(0.0, 0.0, 1.0);
  }
  Vector3d<double> normal = h.normalise();
  Vector3d<double> radial = normal.gaiseki(tangent);
  if (radial.magnitude() < 1e-12) {
    radial = normal.gaiseki(Vector3d<double>(1.0, 0.0, 0.0));
  }
  if (radial.magnitude() < 1e-12) {
    radial = normal.gaiseki(Vector3d<double>(0.0, 1.0, 0.0));
  }
  radial = radial.normalise();

  Vector3d<double> in_plane = tangent * std::cos(theta_rad) + radial * std::sin(theta_rad);
  if (in_plane.magnitude() < 1e-12) in_plane = tangent;
  Vector3d<double> dv_dir = in_plane.normalise() * std::cos(phi_rad) + normal * std::sin(phi_rad);
  if (dv_dir.magnitude() < 1e-12) dv_dir = tangent;
  dv_dir = dv_dir.normalise();

  State<double> state_rot = base;
  state_rot.vx += dv_nd * dv_dir.x();
  state_rot.vy += dv_nd * dv_dir.y();
  state_rot.vz += dv_nd * dv_dir.z();
  output->dv_direction = dv_dir;

  auto mark_filtered = [&](double jacobi_value) {
    output->survival_days = 0.0;
    output->sali_decay_days = 0.0;
    output->final_sali = 0.0;
    output->score = kFilteredScore;
    output->reason = EndReason::FILTERED;
    output->final_state = state_rot;
    output->jacobi_constant = jacobi_value;
  };

  if (!CheckPerigeeFilter(state_rot, mu, collision_r)) {
    double jacobi_val = crtbp::calc_jacobi_integral(state_rot, mu);
    mark_filtered(jacobi_val);
    return true;
  }

  double jacobi_val = crtbp::calc_jacobi_integral(state_rot, mu);
  output->jacobi_constant = jacobi_val;
  const bool jacobi_gate_enabled = std::isfinite(jacobi_min) || std::isfinite(jacobi_max);
  if (jacobi_gate_enabled && (jacobi_val < jacobi_min || jacobi_val > jacobi_max)) {
    mark_filtered(jacobi_val);
    return true;
  }

  SaliState<double> sali_state;
  sali_state.state = ConvertToCanonical(state_rot);
  sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  double t = 0.0;
  EndReason reason = EndReason::TIME_UP;
  double sali_decay_time = max_duration_nd;
  bool sali_dropped = false;
  double current_sali = 2.0;

  while (t < max_duration_nd) {
    crtbp::SymplecticStep6thOrderSALI(mu, &sali_state, dt_nd);
    t += dt_nd;

    const double norm_plus = (sali_state.w1 + sali_state.w2).Norm();
    const double norm_minus = (sali_state.w1 - sali_state.w2).Norm();
    current_sali = std::min(norm_plus, norm_minus);

    if (!sali_dropped && current_sali < sali_threshold) {
      sali_decay_time = t;
      sali_dropped = true;
    }

    State<double> phys = ConvertToPhysical(sali_state.state);
    double r2 = calc_r2(phys.x, phys.y, phys.z, mu);

    if (r2 < collision_r) {
      reason = EndReason::COLLISION;
      break;
    }
    if (r2 > escape_r) {
      reason = EndReason::ESCAPE;
      break;
    }
  }

  const double clamped_time = std::min(t, max_duration_nd);
  output->survival_days = clamped_time;
  output->sali_decay_days =
      sali_dropped ? std::min(sali_decay_time, max_duration_nd) : clamped_time;
  output->final_sali = current_sali;
  output->reason = reason;
  output->final_state = ConvertToPhysical(sali_state.state);

  auto score_term = [](double value) { return std::log10(1.0 + std::max(0.0, value)); };
  double survive_term = score_term(output->survival_days);
  double sali_term = score_term(output->sali_decay_days);
  double reason_bonus = 0.0;
  if (reason == EndReason::TIME_UP) {
    reason_bonus = 0.5;
  } else if (reason == EndReason::ESCAPE) {
    reason_bonus = -0.5;
  } else if (reason == EndReason::COLLISION) {
    reason_bonus = -1.0;
  }
  output->score = survive_term + sali_term + reason_bonus;

  return true;
}void WriteResultCsv(const std::string& filepath, const NaturalUnits& units, const StickyConfig& cfg,
                    const std::vector<SearchTask>& tasks, double dv_nd) {
  std::ofstream ofs(filepath);
  ofs << std::fixed << std::setprecision(10);
  auto format_bound = [](double value) {
    if (std::isfinite(value)) return std::to_string(value);
    return value > 0 ? std::string("+inf") : std::string("-inf");
  };

  ofs << "# Sticky Orbit Search Results\n";
  ofs << "# DV=" << cfg.dv_magnitude_mps << " m/s\n";
  ofs << "# Collision=" << cfg.collision_radius_km << " km, Escape=" << cfg.escape_radius_hill
      << " Hill\n";
  ofs << "# Jacobi gate: [" << format_bound(cfg.jacobi_min) << ", " << format_bound(cfg.jacobi_max)
      << "]\n";
  ofs << "time_days,theta_deg,phi_deg,survival_days,sali_decay_days,final_sali,score,end_reason,"
         "jacobi_constant,final_r2,final_x,final_y,final_z,dv_dir_x,dv_dir_y,dv_dir_z\n";

  const double escape_r_nd = units.hill_radius * cfg.escape_radius_hill;
  const double max_dur_nd = cfg.survival_max_days / units.tu_day;
  const double dt_nd = cfg.survival_dt_days / units.tu_day;

  // 結果をバッファリングするためのベクタ
  std::vector<SearchResult> results(tasks.size());
  std::vector<const TrajectorySample*> sample_order;
  sample_order.reserve(tasks.size());
  std::unordered_set<const TrajectorySample*> seen_samples;
  for (const auto& task : tasks) {
    if (task.sample && seen_samples.insert(task.sample).second) {
      sample_order.push_back(task.sample);
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < static_cast<int>(tasks.size()); ++i) {
    double th_rad = tasks[i].theta_deg * kPi / 180.0;
    double ph_rad = tasks[i].phi_deg * kPi / 180.0;
    EvaluateStickyOrbit(tasks[i], th_rad, ph_rad, dv_nd, units.mu, escape_r_nd,
                        units.collision_radius, cfg.jacobi_min, cfg.jacobi_max, max_dur_nd, dt_nd,
                        cfg.sali_threshold, &results[i]);

    // ND -> Days conversion for output
    if (results[i].reason != EndReason::FILTERED) {
      results[i].survival_days *= units.tu_day;
      results[i].sali_decay_days *= units.tu_day;
    }
  }

  // シリアル書き込み (Gnuplot pm3d用に空行を入れる)
  std::unordered_map<const TrajectorySample*, SearchResult> best_per_time;
  best_per_time.reserve(sample_order.size());
  for (size_t i = 0; i < tasks.size(); ++i) {
    const TrajectorySample* sample = tasks[i].sample;
    if (!sample) continue;
    const auto& candidate = results[i];
    auto it = best_per_time.find(sample);
    if (it == best_per_time.end() || candidate.score > it->second.score) {
      best_per_time[sample] = candidate;
    }
  }
  double prev_t = -1.0;
  for (const auto& res : results) {
    if (prev_t >= 0.0 && std::abs(res.t_days - prev_t) > 1e-6) {
      ofs << "\n";  // Block separator for pm3d
    }
    prev_t = res.t_days;

    ofs << res.t_days << "," << res.theta_deg << "," << res.phi_deg << "," << res.survival_days
        << "," << res.sali_decay_days << "," << res.final_sali << "," << res.score << ","
        << EndReasonToString(res.reason) << "," << res.jacobi_constant << ","
        << calc_r2(res.final_state.x, res.final_state.y, res.final_state.z, units.mu) << ","
        << res.final_state.x << "," << res.final_state.y << "," << res.final_state.z << ","
        << res.dv_direction.x() << "," << res.dv_direction.y() << "," << res.dv_direction.z()
        << "\n";
  }

  std::string best_filepath = filepath;
  const auto dot_pos = best_filepath.find_last_of('.');
  if (dot_pos == std::string::npos) {
    best_filepath.append("_best.csv");
  } else {
    best_filepath.insert(dot_pos, "_best");
  }

  std::ofstream ofs_best(best_filepath);
  ofs_best << std::fixed << std::setprecision(10);
  ofs_best << "# Sticky Orbit Search | Best dV per strip sample\n";
  ofs_best << "# DV=" << cfg.dv_magnitude_mps << " m/s\n";
  ofs_best << "# Collision=" << cfg.collision_radius_km << " km, Escape=" << cfg.escape_radius_hill
           << " Hill\n";
  ofs_best << "time_days,theta_deg,phi_deg,score,survival_days,sali_decay_days,end_reason,"
              "jacobi_constant,dv_dir_x,dv_dir_y,dv_dir_z\n";
  for (const auto* sample : sample_order) {
    auto it = best_per_time.find(sample);
    if (it == best_per_time.end()) continue;
    const auto& best = it->second;
    ofs_best << best.t_days << "," << best.theta_deg << "," << best.phi_deg << ","
             << best.score << "," << best.survival_days << "," << best.sali_decay_days << ","
             << EndReasonToString(best.reason) << "," << best.jacobi_constant << ","
             << best.dv_direction.x() << "," << best.dv_direction.y() << ","
             << best.dv_direction.z() << "\n";
  }
}
