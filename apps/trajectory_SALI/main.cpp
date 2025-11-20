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
using crtbp::IntegrateDopri5Orbit;
using crtbp::IntegrateDopri5SALI;
using crtbp::SaliState;
using crtbp::State;
using crtbp::State3d;
using utils::loadConstants;

namespace fs = std::filesystem;
constexpr double kPi = std::numbers::pi;

struct TrajectorySaliConfig {
  State<double> asteroid_state_helio{};
  State<double> earth_state_helio{};

  // Natural trajectory settings (days in inertial time).
  double natural_duration_days = 5.0;
  double natural_dt_days = 1e-3;
  int sample_stride = 1;

  // Delta-v map settings.
  double dv_magnitude_mps = 100.0;
  double theta_start_deg = -180.0;
  double theta_end_deg = 180.0;
  double theta_step_deg = 5.0;

  // SALI integration after the impulse (still in CRTBP non-dimensional time).
  double sali_duration_days = 5.0;
  double sali_dt_days = 1e-4;
  double forbidden_radius = 1e-8;  // dimensionless distance to avoid singularity
  bool limit_to_hill = true;

  std::string output_tag;
};

struct NaturalUnits {
  double mu = 0.0;
  double lu_au = 0.0;
  double tu_day = 0.0;
  double vu_au_per_day = 0.0;
  double vu_m_per_s = 0.0;
  double hill_radius = 0.0;  // dimensionless Hill radius around P2
};

struct TrajectorySample {
  double t_nd = 0.0;
  double t_days = 0.0;
  State<double> state{};
  double r2 = 0.0;  // distance from Earth (dimensionless)
};

struct MapTask {
  const TrajectorySample* sample = nullptr;
  double theta_deg = 0.0;
};

struct MapResult {
  double t_nd = 0.0;
  double t_days = 0.0;
  double theta_deg = 0.0;
  double base_r2 = 0.0;
  double final_r2 = 0.0;
  double sali = std::numeric_limits<double>::quiet_NaN();
  State<double> base_state{};
  State<double> post_state{};
};

// Forward declarations
bool LoadTrajectoryConfig(const std::string& filepath, TrajectorySaliConfig* config);

NaturalUnits ComputeNaturalUnits(const State<double>& earth_state_helio,
                                 const AstroConstants<double>& astro);

std::vector<TrajectorySample> IntegrateNaturalTrajectory(const AstroConstants<double>& astro,
                                                         const State<double>& initial_state_nd,
                                                         double duration_nd, double dt_nd,
                                                         double mu, const NaturalUnits& units);

bool EvaluateSali(const MapTask& task, double theta_rad, double dv_nd, double mu,
                  double hill_radius, double forbidden_radius, double sali_duration_nd,
                  double sali_dt_nd, const AstroConstants<double>& astro, MapResult* output);

void WriteTrajectoryCsv(const std::string& filepath, const std::vector<TrajectorySample>& samples,
                        double hill_radius);

void WriteMapCsv(const std::string& filepath, const NaturalUnits& units,
                 const TrajectorySaliConfig& cfg, double sali_duration_nd, double sali_dt_nd,
                 const std::vector<MapTask>& tasks, double dv_nd, double hill_radius,
                 double forbidden_radius, const AstroConstants<double>& astro);

std::vector<std::string> DiscoverConfigFiles(const std::string& directory,
                                             const std::string& prefix);

int main() {
  using namespace param;

  const std::string kConfigDir = std::string(CONFIG_DIR) + "/trajectory_SALI/";
  const std::string kOutputDir = std::string(OUTPUT_DIR) + "/trajectory_SALI";
  const std::string kConfigPrefix = "trajectorySALI";

  AstroConstants<double> astro_params =
      loadConstants<double>(kConfigDir + "/../astro_param/astro_param.txt");
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  std::cout << "<>------------------------------------------------------------\n";
  std::cout << "<>  CRTBP SALI map along natural trajectory (Earth Hill sphere)\n";
  std::cout << "<>------------------------------------------------------------\n\n";

  std::vector<std::string> config_files = DiscoverConfigFiles(kConfigDir, kConfigPrefix);
  if (config_files.empty()) {
    std::cerr << "No config files found under " << kConfigDir << " with prefix " << kConfigPrefix
              << "\n";
    return 1;
  }

  fs::create_directories(kOutputDir);

  int config_index = 0;
  for (const auto& config_path : config_files) {
    ++config_index;
    TrajectorySaliConfig cfg;
    if (!LoadTrajectoryConfig(config_path, &cfg)) {
      std::cerr << "Skipping config due to parse errors: " << config_path << "\n";
      continue;
    }

    NaturalUnits units{};
    try {
      units = ComputeNaturalUnits(cfg.earth_state_helio, astro_params);
    } catch (const std::exception& e) {
      std::cerr << "Failed to compute natural units for " << config_path << ": " << e.what()
                << "\n";
      continue;
    }
    std::cout << "<> Config #" << config_index << ": " << config_path << "\n";
    std::cout << "    mu=" << std::setprecision(8) << units.mu
              << ", Hill radius=" << std::setprecision(6) << units.hill_radius
              << ", LU(AU)=" << units.lu_au << ", TU(days)=" << units.tu_day
              << ", VU(m/s)=" << units.vu_m_per_s << "\n";

    State<double> initial_state_nd{};
    try {
      initial_state_nd = crtbp::ConvertInertial2Rotating____(cfg.asteroid_state_helio,
                                                             cfg.earth_state_helio, astro_params);
    } catch (const std::exception& e) {
      std::cerr << "Failed to convert heliocentric states to CRTBP frame: " << e.what() << "\n";
      continue;
    }

    std::cout << "  Initial state (before conberting): "
              << cfg.asteroid_state_helio.x - cfg.earth_state_helio.x << " "
              << cfg.asteroid_state_helio.y - cfg.earth_state_helio.y << " "
              << cfg.asteroid_state_helio.z - cfg.earth_state_helio.z << " "
              << cfg.asteroid_state_helio.vx - cfg.earth_state_helio.vx << " "
              << cfg.asteroid_state_helio.vy - cfg.earth_state_helio.vy << " "
              << cfg.asteroid_state_helio.vz - cfg.earth_state_helio.vz << "\n";
    std::cout << "  Initial state (nd): " << initial_state_nd.x << " " << initial_state_nd.y << " "
              << initial_state_nd.z << " " << initial_state_nd.vx << " " << initial_state_nd.vy
              << " " << initial_state_nd.vz << "\n";
    std::cout << " Jacobi constant: " << crtbp::calc_jacobi_integral(initial_state_nd, units.mu)
              << "\n";
    double target_jacobi = 3.0;
    std::cout << " Target Jacobi constant: " << target_jacobi << "\n";
    State3d<double> state(initial_state_nd.x, initial_state_nd.y, initial_state_nd.z);
    double needed_v_abs =
        crtbp::calc_v_abs(
            State3d<double>(initial_state_nd.x, initial_state_nd.y, initial_state_nd.z),
            target_jacobi, units.mu) -
        crtbp::calc_v_abs(
            State3d<double>(initial_state_nd.x, initial_state_nd.y, initial_state_nd.z),
            crtbp::calc_jacobi_integral(initial_state_nd, units.mu), units.mu);
    std::cout << " Needed delta-v to reach target Jacobi: " << needed_v_abs * units.vu_m_per_s
              << " m/s\n";
    utils::WaitForEnter();

    const double natural_duration_nd = cfg.natural_duration_days / units.tu_day;
    const double natural_dt_nd = cfg.natural_dt_days / units.tu_day;
    if (natural_duration_nd <= 0 || natural_dt_nd <= 0) {
      std::cerr << "Invalid natural trajectory time settings. duration_nd=" << natural_duration_nd
                << " dt_nd=" << natural_dt_nd << "\n";
      continue;
    }

    auto samples = IntegrateNaturalTrajectory(astro_params, initial_state_nd, natural_duration_nd,
                                              natural_dt_nd, kMU, units);
    if (samples.empty()) {
      std::cerr << "Natural trajectory integration produced no samples.\n";
      continue;
    }

    // Identify the segment inside the Hill sphere.
    size_t start_idx = 0;
    size_t end_idx = samples.size() - 1;
    if (cfg.limit_to_hill) {
      bool found = false;
      for (size_t i = 0; i < samples.size(); ++i) {
        if (samples[i].r2 <= 3.0 * units.hill_radius) {
          start_idx = i;
          found = true;
          break;
        }
      }
      for (size_t i = samples.size(); i-- > 0;) {
        if (samples[i].r2 <= 3.0 * units.hill_radius) {
          end_idx = i;
          found = true;
          break;
        }
      }
      if (!found) {
        std::cout << "  Warning: orbit never enters Hill sphere. Using full trajectory.\n";
        start_idx = 0;
        end_idx = samples.size() - 1;
      }
    }

    // Select samples along the arc (with stride).
    std::vector<MapTask> tasks;
    const size_t stride = static_cast<size_t>(std::max(1, cfg.sample_stride));
    const double theta_start = cfg.theta_start_deg;
    const double theta_end = cfg.theta_end_deg;
    const double theta_step = cfg.theta_step_deg > 0 ? cfg.theta_step_deg : 1.0;

    for (size_t idx = start_idx; idx <= end_idx; idx += stride) {
      for (double th = theta_start; th <= theta_end + 1e-9; th += theta_step) {
        tasks.push_back(MapTask{&samples[idx], th});
      }
    }

    std::cout << "  Using samples [" << start_idx << ", " << end_idx << "] (total "
              << samples.size() << ")\n";
    std::cout << "  Map tasks: " << tasks.size() << " (theta step " << theta_step << " deg)\n";

    // Dump natural trajectory and SALI map.
    const std::string tag =
        cfg.output_tag.empty() ? fs::path(config_path).stem().string() : cfg.output_tag;
    const std::string trajectory_file = kOutputDir + "/" + tag + "_trajectory.csv";
    const std::string map_file = kOutputDir + "/" + tag + "_sali_map.csv";

    WriteTrajectoryCsv(trajectory_file, samples, units.hill_radius);

    const double sali_duration_nd = cfg.sali_duration_days / units.tu_day;
    const double sali_dt_nd = cfg.sali_dt_days / units.tu_day;
    const double dv_nd = cfg.dv_magnitude_mps / units.vu_m_per_s;
    WriteMapCsv(map_file, units, cfg, sali_duration_nd, sali_dt_nd, tasks, dv_nd, units.hill_radius,
                cfg.forbidden_radius, astro_params);

    std::cout << "  -> trajectory: " << trajectory_file << "\n";
    std::cout << "  -> SALI map : " << map_file << "\n";
  }

  std::cout << "\nAll configs processed.\n";
  return 0;
}

bool LoadTrajectoryConfig(const std::string& filepath, TrajectorySaliConfig* config) {
  if (config == nullptr) return false;

  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "Cannot open config: " << filepath << "\n";
    return false;
  }

  std::string line;
  auto trim = [](const std::string& s) -> std::string {
    const auto start = s.find_first_not_of(" \t");
    const auto end = s.find_last_not_of(" \t");
    if (start == std::string::npos) {
      std::cout << "  Warning: empty or whitespace-only line encountered during config parsing.\n";
      return "";
    };

    return s.substr(start, end - start + 1);
  };

  auto parse_state = [](const std::string& rhs, State<double>* out) -> bool {
    if (out == nullptr) return false;
    std::stringstream ss(rhs);
    std::array<double, 6> values{};
    for (int i = 0; i < 6; ++i) {
      if (!(ss >> values[i])) return false;
    }
    *out = State<double>{values[0], values[1], values[2], values[3], values[4], values[5]};
    return true;
  };

  bool earth_found = false;
  bool asteroid_found = false;

  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    const auto pos = line.find('=');
    if (pos == std::string::npos) continue;
    const std::string key = trim(line.substr(0, pos));
    const std::string value = trim(line.substr(pos + 1));

    if (key == "ASTEROID_STATE_HELIO_AU" || key == "ASTEROID_STATE_AU") {
      if (!parse_state(value, &config->asteroid_state_helio)) {
        std::cerr << "  Failed to parse ASTEROID_STATE_HELIO_AU\n";
        return false;
      }
      asteroid_found = true;
    } else if (key == "EARTH_STATE_HELIO_AU" || key == "EARTH_STATE_AU") {
      if (!parse_state(value, &config->earth_state_helio)) {
        std::cerr << "  Failed to parse EARTH_STATE_HELIO_AU\n";
        return false;
      }
      earth_found = true;
    } else if (key == "NATURAL_DURATION_DAYS") {
      config->natural_duration_days = std::stod(value);
    } else if (key == "NATURAL_DT_DAYS") {
      config->natural_dt_days = std::stod(value);
    } else if (key == "MAP_SAMPLE_STRIDE") {
      config->sample_stride = std::stoi(value);
    } else if (key == "DV_MPS") {
      config->dv_magnitude_mps = std::stod(value);
    } else if (key == "THETA_RANGE_DEG") {
      std::stringstream ss(value);
      ss >> config->theta_start_deg >> config->theta_end_deg >> config->theta_step_deg;
    } else if (key == "SALI_DURATION_DAYS") {
      config->sali_duration_days = std::stod(value);
    } else if (key == "SALI_DT_DAYS") {
      config->sali_dt_days = std::stod(value);
    } else if (key == "FORBIDDEN_RADIUS") {
      config->forbidden_radius = std::stod(value);
    } else if (key == "LIMIT_TO_HILL") {
      config->limit_to_hill = (std::stoi(value) != 0);
    } else if (key == "OUTPUT_TAG") {
      config->output_tag = value;
    }
  }

  if (!earth_found || !asteroid_found) {
    std::cerr << "  Missing Earth/Asteroid state definitions in " << filepath << "\n";
    return false;
  }

  return true;
}

NaturalUnits ComputeNaturalUnits(const State<double>& earth_state_helio,
                                 const AstroConstants<double>& astro) {
  NaturalUnits units{};
  units.mu = astro.gm_earth / (astro.gm_earth + astro.gm_sun);

  const double kM3s2_to_AU3day2 =
      (kSecondsPerDay * kSecondsPerDay) / (astro.au * astro.au * astro.au);
  const double gm_sun_AD = astro.gm_sun * kM3s2_to_AU3day2;
  const double gm_earth_AD = astro.gm_earth * kM3s2_to_AU3day2;
  const double gm_total_AD = gm_sun_AD + gm_earth_AD;

  const Vector3d<double> r_p2_G{earth_state_helio.x, earth_state_helio.y, earth_state_helio.z};
  if (r_p2_G.magnitude() <= 0.0) {
    throw std::runtime_error("Earth position magnitude must be positive to build CRTBP units.");
  }
  const double mean_motion_AD = std::sqrt(gm_total_AD);  // LU = 1 AU -> r = 1
  const double tu_day = 1.0 / mean_motion_AD;
  const double vu_au_day = 1.0 / tu_day;  // LU = 1 AU

  units.lu_au = 1.0;
  units.tu_day = tu_day;
  units.vu_au_per_day = vu_au_day;
  units.vu_m_per_s = vu_au_day * astro.au / kSecondsPerDay;
  units.hill_radius = std::cbrt(units.mu / 3.0);  // normalized Earth Hill sphere radius
  return units;
}

std::vector<TrajectorySample> IntegrateNaturalTrajectory(const AstroConstants<double>& astro,
                                                         const State<double>& initial_state_nd,
                                                         double duration_nd, double dt_nd,
                                                         double mu, const NaturalUnits& units) {
  std::vector<TrajectorySample> samples;
  const int steps = static_cast<int>(std::ceil(duration_nd / dt_nd));
  samples.reserve(static_cast<size_t>(steps) + 1);

  auto observer = [&](const State<double>& state, double t) {
    const double r2 = calc_r2(state.x, state.y, state.z, mu);
    samples.push_back(TrajectorySample{t, t * units.tu_day, state, r2});
  };

  State<double> work_state = initial_state_nd;
  IntegrateDopri5Orbit(astro, work_state, 0.0, duration_nd, dt_nd, observer);
  return samples;
}

class SaliRecorder {
 public:
  void operator()(const SaliState<double>& state, double /*t*/) {
    const double norm_plus = (state.w1 + state.w2).Norm();
    const double norm_minus = (state.w1 - state.w2).Norm();
    last_sali_ = std::min(norm_plus, norm_minus);
  }

  double last_sali() const { return last_sali_; }

 private:
  double last_sali_ = std::numeric_limits<double>::quiet_NaN();
};

bool EvaluateSali(const MapTask& task, double theta_rad, double dv_nd, double mu,
                  double hill_radius, double forbidden_radius, double sali_duration_nd,
                  double sali_dt_nd, const AstroConstants<double>& astro, MapResult* output) {
  if (output == nullptr || task.sample == nullptr) return false;

  const State<double>& base_state = task.sample->state;
  Vector3d<double> v_vec{base_state.vx, base_state.vy, base_state.vz};
  const double speed = v_vec.magnitude();
  if (speed < 1e-10) {
    return false;
  }
  Vector3d<double> tangent = v_vec.normalise();

  // Build an in-plane basis (tangent, radial).
  Vector3d<double> r_to_earth{base_state.x - (1.0 - mu), base_state.y, base_state.z};
  Vector3d<double> h_vec = r_to_earth.gaiseki(v_vec);
  if (h_vec.magnitude() < 1e-12) {
    // Degenerate: pick an arbitrary normal perpendicular to tangent.
    if (std::abs(tangent.x()) < 0.9) {
      h_vec = Vector3d<double>(0.0, 0.0, 1.0);
    } else {
      h_vec = Vector3d<double>(0.0, 1.0, 0.0);
    }
  }
  Vector3d<double> normal = h_vec.normalise();
  Vector3d<double> radial = normal.gaiseki(tangent).normalise();

  const double c = std::cos(theta_rad);
  const double s = std::sin(theta_rad);
  const Vector3d<double> dv_dir = tangent * c + radial * s;

  State<double> kicked_state = base_state;
  kicked_state.vx += dv_nd * dv_dir.x();
  kicked_state.vy += dv_nd * dv_dir.y();
  kicked_state.vz += dv_nd * dv_dir.z();

  SaliState<double> sali_state;
  sali_state.state = ConvertToCanonical(kicked_state);
  sali_state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  sali_state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  SaliRecorder recorder;
  IntegrateDopri5SALI(astro, sali_state, 0.0, sali_duration_nd, sali_dt_nd, recorder, true);

  const State<double> final_state = ConvertToPhysical(sali_state.state);
  const double final_r2 = calc_r2(final_state.x, final_state.y, final_state.z, mu);

  output->t_nd = task.sample->t_nd;
  output->t_days = task.sample->t_days;
  output->theta_deg = task.theta_deg;
  output->base_r2 = task.sample->r2;
  output->final_r2 = final_r2;
  output->sali = recorder.last_sali();
  output->base_state = base_state;
  output->post_state = final_state;

  // If the orbit after the impulse escapes or collides, flag with NaN.
  if (final_r2 > hill_radius || final_r2 < forbidden_radius) {
    output->sali = std::numeric_limits<double>::quiet_NaN();
  }
  return true;
}

void WriteTrajectoryCsv(const std::string& filepath, const std::vector<TrajectorySample>& samples,
                        double hill_radius) {
  std::ofstream ofs(filepath);
  if (!ofs) {
    std::cerr << "Failed to open trajectory output: " << filepath << "\n";
    return;
  }
  ofs << std::fixed << std::setprecision(15);
  ofs << "# time_days,time_nd,x,y,z,vx,vy,vz,r2,inside_hill\n";
  for (const auto& s : samples) {
    const bool inside = (s.r2 <= hill_radius);
    ofs << s.t_days << "," << s.t_nd << "," << s.state.x << "," << s.state.y << "," << s.state.z
        << "," << s.state.vx << "," << s.state.vy << "," << s.state.vz << "," << s.r2 << ","
        << (inside ? 1 : 0) << "\n";
  }
}

void WriteMapCsv(const std::string& filepath, const NaturalUnits& units,
                 const TrajectorySaliConfig& cfg, double sali_duration_nd, double sali_dt_nd,
                 const std::vector<MapTask>& tasks, double dv_nd, double hill_radius,
                 double forbidden_radius, const AstroConstants<double>& astro) {
  std::ofstream ofs(filepath);
  if (!ofs) {
    std::cerr << "Failed to open SALI map output: " << filepath << "\n";
    return;
  }

  ofs << std::fixed << std::setprecision(15);
  ofs << "# SALI map along natural trajectory\n";
  ofs << "# DV_mps=" << cfg.dv_magnitude_mps << " (nd=" << dv_nd << ")\n";
  ofs << "# theta_range_deg=" << cfg.theta_start_deg << " " << cfg.theta_end_deg
      << " step=" << cfg.theta_step_deg << "\n";
  ofs << "# sali_duration_nd=" << sali_duration_nd << " sali_dt_nd=" << sali_dt_nd << "\n";
  ofs << "# hill_radius=" << hill_radius << " forbidden_radius=" << forbidden_radius << "\n";
  ofs << "# LU(AU)=" << units.lu_au << " TU(days)=" << units.tu_day
      << " VU(m/s)=" << units.vu_m_per_s << "\n";
  ofs << "time_days,time_nd,theta_deg,base_r2,final_r2,inside_hill,x,y,z,vx,vy,vz,"
         "post_x,post_y,post_z,post_vx,post_vy,post_vz,sali\n";

  std::ostringstream buffer;
  buffer << std::fixed << std::setprecision(15);

#pragma omp parallel
  {
    std::ostringstream local;
    local << std::fixed << std::setprecision(15);

#pragma omp for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(tasks.size()); ++i) {
      const MapTask& task = tasks[static_cast<size_t>(i)];
      MapResult result;
      const double theta_rad = task.theta_deg * kPi / 180.0;
      if (!EvaluateSali(task, theta_rad, dv_nd, units.mu, 3.0 * hill_radius, forbidden_radius,
                        sali_duration_nd, sali_dt_nd, astro, &result)) {
        continue;
      }

      const bool inside_final =
          result.final_r2 <= hill_radius && result.final_r2 >= forbidden_radius;
      local << result.t_days << "," << result.t_nd << "," << result.theta_deg << ","
            << result.base_r2 << "," << result.final_r2 << "," << (inside_final ? 1 : 0) << ","
            << result.base_state.x << "," << result.base_state.y << "," << result.base_state.z
            << "," << result.base_state.vx << "," << result.base_state.vy << ","
            << result.base_state.vz << "," << result.post_state.x << "," << result.post_state.y
            << "," << result.post_state.z << "," << result.post_state.vx << ","
            << result.post_state.vy << "," << result.post_state.vz << "," << result.sali << "\n";
    }

#pragma omp critical
    buffer << local.str();
  }

  ofs << buffer.str();
}

std::vector<std::string> DiscoverConfigFiles(const std::string& directory,
                                             const std::string& prefix) {
  std::vector<std::string> files;
  if (!fs::exists(directory)) {
    std::cerr << "Config directory does not exist: " << directory << "\n";
    return files;
  }
  const std::regex pattern("^" + prefix + "(?:_\\d+)?\\.txt$");
  try {
    for (const auto& entry : fs::directory_iterator(directory)) {
      if (!entry.is_regular_file()) continue;
      const auto name = entry.path().filename().string();
      if (std::regex_match(name, pattern)) {
        files.push_back(fs::absolute(entry.path()).string());
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "Failed to read config directory: " << e.what() << "\n";
    return files;
  }

  auto sorter = [](const std::string& a, const std::string& b) {
    const std::string stem_a = fs::path(a).stem().string();
    const std::string stem_b = fs::path(b).stem().string();
    auto number_from_stem = [](const std::string& stem) -> int {
      const auto pos = stem.find_last_of('_');
      if (pos == std::string::npos) return 0;
      return std::stoi(stem.substr(pos + 1));
    };
    return number_from_stem(stem_a) < number_from_stem(stem_b);
  };
  std::sort(files.begin(), files.end(), sorter);
  return files;
}
