#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rtbp.hpp"
#include "utils.hpp"

namespace fs = std::filesystem;

using crtbp::AstroConstants;
using crtbp::Dopri5OrbitSystem;
using crtbp::State;
using utils::DiscoverConfigFiles;
using utils::getcurrent_date;
using utils::loadConstants;
using utils::trim;

struct PoincareConfig {
  double jacobi_constant = 3.0009;
  double x_min = 0.988;
  double x_max = 0.999;
  int x_samples = 80;
  double x_step = 0.0;
  double initial_vx = 0.0;
  double initial_vz = 0.0;
  double initial_z = 0.0;
  double max_time = 300.0;
  double time_step = 1e-3;
  int max_crossings_per_traj = 200;
  double plane_velocity_threshold = 1e-7;
  double crossing_tolerance = 1e-10;
  double earth_collision_radius = 1e-4;
  double sun_collision_radius = 1e-4;
  double escape_radius = 2.0;
  int plane_direction = 1;
  bool record_initial_crossing = true;
  std::string output_basename = "poincare_map";
  bool append_timestamp = true;
};

struct CrossingRecord {
  int sample_index = 0;
  int crossing_index = 0;
  double time = 0.0;
  double initial_x = 0.0;
  double initial_vx = 0.0;
  double initial_vy = 0.0;
  double jacobi = 0.0;
  double jacobi_error = 0.0;
  State<double> state{};
};

enum class TerminationReason {
  kMaxTime,
  kMaxCrossings,
  kCollisionEarth,
  kCollisionSun,
  kEscape,
  kEnergyInvalid,
  kCompleted
};

struct CrossingEvent {
  State<double> state{};
  double time = 0.0;
  double jacobi = 0.0;
  double jacobi_error = 0.0;
};

bool LoadConfig(const std::string& filepath, PoincareConfig* cfg) {
  std::ifstream ifs(filepath);
  if (!ifs.is_open()) {
    std::cerr << "Failed to open config: " << filepath << "\n";
    return false;
  }

  std::string line;
  while (std::getline(ifs, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    if (line.empty()) continue;
    const auto eq = line.find('=');
    if (eq == std::string::npos) continue;

    std::string key = trim(line.substr(0, eq));
    std::string value = trim(line.substr(eq + 1));
    if (key.empty() || value.empty()) continue;

    try {
      if (key == "jacobi_constant") {
        cfg->jacobi_constant = std::stod(value);
      } else if (key == "x_min") {
        cfg->x_min = std::stod(value);
      } else if (key == "x_max") {
        cfg->x_max = std::stod(value);
      } else if (key == "x_samples") {
        cfg->x_samples = std::stoi(value);
      } else if (key == "x_step") {
        cfg->x_step = std::stod(value);
      } else if (key == "initial_vx") {
        cfg->initial_vx = std::stod(value);
      } else if (key == "initial_vz") {
        cfg->initial_vz = std::stod(value);
      } else if (key == "initial_z") {
        cfg->initial_z = std::stod(value);
      } else if (key == "max_time") {
        cfg->max_time = std::stod(value);
      } else if (key == "time_step") {
        cfg->time_step = std::stod(value);
      } else if (key == "max_crossings_per_traj") {
        cfg->max_crossings_per_traj = std::stoi(value);
      } else if (key == "plane_velocity_threshold") {
        cfg->plane_velocity_threshold = std::stod(value);
      } else if (key == "crossing_tolerance") {
        cfg->crossing_tolerance = std::stod(value);
      } else if (key == "earth_collision_radius") {
        cfg->earth_collision_radius = std::stod(value);
      } else if (key == "sun_collision_radius") {
        cfg->sun_collision_radius = std::stod(value);
      } else if (key == "escape_radius") {
        cfg->escape_radius = std::stod(value);
      } else if (key == "plane_direction") {
        cfg->plane_direction = std::stoi(value);
      } else if (key == "record_initial_crossing") {
        cfg->record_initial_crossing = (value == "1" || value == "true" || value == "TRUE");
      } else if (key == "output_basename") {
        cfg->output_basename = value;
      } else if (key == "append_timestamp") {
        cfg->append_timestamp = (value == "1" || value == "true" || value == "TRUE");
      }
    } catch (const std::exception& e) {
      std::cerr << "Failed to parse key '" << key << "': " << e.what() << "\n";
      return false;
    }
  }

  if (cfg->x_max < cfg->x_min) {
    std::swap(cfg->x_max, cfg->x_min);
  }
  if (cfg->x_samples < 1 && cfg->x_step <= 0.0) {
    std::cerr << "Either x_samples >= 1 or x_step > 0 must be provided.\n";
    return false;
  }
  if (cfg->max_time <= 0.0 || cfg->time_step <= 0.0) {
    std::cerr << "max_time and time_step must be positive.\n";
    return false;
  }
  if (cfg->max_crossings_per_traj <= 0) {
    std::cerr << "max_crossings_per_traj must be positive.\n";
    return false;
  }
  if (cfg->escape_radius <= 0.0) {
    std::cerr << "escape_radius must be positive.\n";
    return false;
  }
  if (cfg->plane_direction == 0) {
    cfg->plane_direction = 1;
  } else if (cfg->plane_direction > 0) {
    cfg->plane_direction = 1;
  } else {
    cfg->plane_direction = -1;
  }
  if (cfg->output_basename.empty()) {
    cfg->output_basename = "poincare_map";
  }
  return true;
}

std::vector<double> BuildInitialXGrid(const PoincareConfig& cfg) {
  std::vector<double> xs;
  if (cfg.x_step > 0.0) {
    for (double x = cfg.x_min; x <= cfg.x_max + 0.5 * cfg.x_step; x += cfg.x_step) {
      xs.push_back(x);
    }
  } else if (cfg.x_samples <= 1) {
    xs.push_back(0.5 * (cfg.x_min + cfg.x_max));
  } else {
    const double step = (cfg.x_max - cfg.x_min) / static_cast<double>(cfg.x_samples - 1);
    for (int i = 0; i < cfg.x_samples; ++i) {
      xs.push_back(cfg.x_min + step * static_cast<double>(i));
    }
  }
  return xs;
}

std::optional<double> ComputeInitialVy(double x, double y, double z, double vx, double vz,
                                       double jacobi_constant, double mu, int direction) {
  const double v_sq =
      2.0 * crtbp::calc_potential_U(x, y, z, mu) - jacobi_constant - vx * vx - vz * vz;
  if (v_sq < -1e-13) {
    return std::nullopt;
  }
  const double vy = std::sqrt(std::max(0.0, v_sq));
  return direction >= 0 ? vy : -vy;
}

bool DetectCrossing(const State<double>& prev_state, const State<double>& curr_state,
                    double prev_time, double curr_time, const PoincareConfig& cfg, double mu,
                    double jacobi_target, CrossingEvent* event) {
  const double y0 = prev_state.y;
  const double y1 = curr_state.y;
  if (cfg.plane_direction > 0) {
    if (!(y0 <= 0.0 && y1 >= 0.0)) {
      return false;
    }
  } else {
    if (!(y0 >= 0.0 && y1 <= 0.0)) {
      return false;
    }
  }

  const double denom = y1 - y0;
  double alpha = 0.5;
  if (std::abs(denom) >= cfg.crossing_tolerance) {
    alpha = -y0 / denom;
  }
  alpha = std::clamp(alpha, 0.0, 1.0);

  State<double> cross_state{};
  cross_state.x = prev_state.x + alpha * (curr_state.x - prev_state.x);
  cross_state.y = 0.0;
  cross_state.z = prev_state.z + alpha * (curr_state.z - prev_state.z);
  cross_state.vx = prev_state.vx + alpha * (curr_state.vx - prev_state.vx);
  cross_state.vy = prev_state.vy + alpha * (curr_state.vy - prev_state.vy);
  cross_state.vz = prev_state.vz + alpha * (curr_state.vz - prev_state.vz);

  if (cfg.plane_direction > 0) {
    if (cross_state.vy <= cfg.plane_velocity_threshold) {
      return false;
    }
  } else {
    if (cross_state.vy >= -cfg.plane_velocity_threshold) {
      return false;
    }
  }

  event->state = cross_state;
  event->time = prev_time + alpha * (curr_time - prev_time);
  event->jacobi = crtbp::calc_jacobi_integral(cross_state, mu);
  event->jacobi_error = jacobi_target - event->jacobi;
  return true;
}

std::string_view ToString(TerminationReason reason) {
  switch (reason) {
    case TerminationReason::kMaxTime:
      return "max_time";
    case TerminationReason::kMaxCrossings:
      return "max_crossings";
    case TerminationReason::kCollisionEarth:
      return "earth_collision";
    case TerminationReason::kCollisionSun:
      return "sun_collision";
    case TerminationReason::kEscape:
      return "escape";
    case TerminationReason::kEnergyInvalid:
      return "invalid_initial_condition";
    case TerminationReason::kCompleted:
    default:
      return "completed";
  }
}

TerminationReason SimulateTrajectory(int sample_index, double x0, double vy0,
                                     const PoincareConfig& cfg,
                                     const AstroConstants<double>& astro_params, double mu,
                                     std::vector<CrossingRecord>* global_records) {
  using Stepper = boost::numeric::odeint::runge_kutta_dopri5<
      State<double>, double, State<double>, double,
      boost::numeric::odeint::vector_space_algebra>;

  State<double> current{};
  current.x = x0;
  current.y = 0.0;
  current.z = cfg.initial_z;
  current.vx = cfg.initial_vx;
  current.vy = vy0;
  current.vz = cfg.initial_vz;

  Stepper stepper;
  Dopri5OrbitSystem<double> system(astro_params);

  double current_time = 0.0;
  State<double> prev_state = current;
  double prev_time = current_time;
  bool has_prev = false;
  int recorded_crossings = 0;

  const auto record_event = [&](const CrossingEvent& event) {
    CrossingRecord rec;
    rec.sample_index = sample_index;
    rec.crossing_index = recorded_crossings;
    rec.time = event.time;
    rec.initial_x = x0;
    rec.initial_vx = cfg.initial_vx;
    rec.initial_vy = vy0;
    rec.jacobi = event.jacobi;
    rec.jacobi_error = event.jacobi_error;
    rec.state = event.state;
    global_records->push_back(rec);
  };

  if (cfg.record_initial_crossing && std::abs(current.y) <= cfg.crossing_tolerance) {
    const bool direction_ok =
        (cfg.plane_direction > 0 && current.vy > cfg.plane_velocity_threshold) ||
        (cfg.plane_direction < 0 && current.vy < -cfg.plane_velocity_threshold);
    if (direction_ok) {
      CrossingEvent event;
      event.state = current;
      event.time = current_time;
      event.jacobi = crtbp::calc_jacobi_integral(current, mu);
      event.jacobi_error = cfg.jacobi_constant - event.jacobi;
      record_event(event);
      ++recorded_crossings;
    }
  }

  while (current_time < cfg.max_time) {
    if (recorded_crossings >= cfg.max_crossings_per_traj) {
      return TerminationReason::kMaxCrossings;
    }

    const double dt = std::min(cfg.time_step, cfg.max_time - current_time);
    if (!has_prev) {
      prev_state = current;
      prev_time = current_time;
      has_prev = true;
    }
    stepper.do_step(system, current, current_time, dt);
    current_time += dt;

    const double r_earth = crtbp::calc_r2(current.x, current.y, current.z, mu);
    if (r_earth < cfg.earth_collision_radius) {
      return TerminationReason::kCollisionEarth;
    }
    const double r_sun = crtbp::calc_r1(current.x, current.y, current.z, mu);
    if (r_sun < cfg.sun_collision_radius) {
      return TerminationReason::kCollisionSun;
    }
    const double r_total =
        std::sqrt(current.x * current.x + current.y * current.y + current.z * current.z);
    if (r_total > cfg.escape_radius) {
      return TerminationReason::kEscape;
    }

    CrossingEvent event;
    if (DetectCrossing(prev_state, current, prev_time, current_time, cfg, mu,
                       cfg.jacobi_constant, &event)) {
      record_event(event);
      ++recorded_crossings;
      if (recorded_crossings >= cfg.max_crossings_per_traj) {
        return TerminationReason::kMaxCrossings;
      }
    }

    prev_state = current;
    prev_time = current_time;
  }

  return TerminationReason::kMaxTime;
}

void WriteCsv(const fs::path& output_path, const PoincareConfig& cfg,
              const std::vector<CrossingRecord>& records) {
  std::ofstream ofs(output_path);
  if (!ofs.is_open()) {
    throw std::runtime_error("Failed to open output file: " + output_path.string());
  }
  ofs << std::setprecision(15);
  ofs << "# JacobiConstant=" << cfg.jacobi_constant << "\n";
  ofs << "# XRange=[" << cfg.x_min << "," << cfg.x_max << "]\n";
  ofs << "# InitialVX=" << cfg.initial_vx << ", InitialVZ=" << cfg.initial_vz << "\n";
  ofs << "# MaxTime=" << cfg.max_time << ", TimeStep=" << cfg.time_step << "\n";
  ofs << "# MaxCrossingsPerTrajectory=" << cfg.max_crossings_per_traj << "\n";
  ofs << "sample_index,crossing_index,time,x,y,z,vx,vy,vz,initial_x,initial_vx,initial_vy,"
         "jacobi,jacobi_error\n";

  for (const auto& rec : records) {
    ofs << rec.sample_index << "," << rec.crossing_index << "," << rec.time << ","
        << rec.state.x << "," << rec.state.y << "," << rec.state.z << "," << rec.state.vx << ","
        << rec.state.vy << "," << rec.state.vz << "," << rec.initial_x << "," << rec.initial_vx
        << "," << rec.initial_vy << "," << rec.jacobi << "," << rec.jacobi_error << "\n";
  }
}

int main(int argc, char** argv) {
  try {
    const fs::path config_dir = fs::path(CONFIG_DIR) / "poincare_map";
    const fs::path output_dir = fs::path(OUTPUT_DIR) / "poincare_map";
    fs::create_directories(output_dir);

    std::string config_path;
    if (argc > 1) {
      config_path = argv[1];
    } else {
      const auto files = DiscoverConfigFiles(config_dir.string(), "poincare");
      if (files.empty()) {
        std::cerr << "No config files found under " << config_dir << "\n";
        return 1;
      }
      config_path = files.front();
      std::cout << "Using config: " << config_path << "\n";
    }

    PoincareConfig config;
    if (!LoadConfig(config_path, &config)) {
      return 1;
    }

    const auto astro = loadConstants<double>(std::string(CONFIG_DIR) + "/astro_param/astro_param.txt");
    const double mu = astro.gm_earth / (astro.gm_earth + astro.gm_sun);

    const auto x_values = BuildInitialXGrid(config);
    std::cout << "Sweeping " << x_values.size() << " initial points.\n";

    std::vector<CrossingRecord> records;
    records.reserve(static_cast<size_t>(x_values.size()) * 32);
    std::map<TerminationReason, int> reason_counts;

    for (size_t idx = 0; idx < x_values.size(); ++idx) {
      const double x0 = x_values[idx];
      auto vy0 = ComputeInitialVy(x0, 0.0, config.initial_z, config.initial_vx, config.initial_vz,
                                  config.jacobi_constant, mu, config.plane_direction);

      if (!vy0.has_value()) {
        reason_counts[TerminationReason::kEnergyInvalid] += 1;
        std::cout << "[skip] x0=" << x0 << " violates Jacobi constraint.\n";
        continue;
      }

      const TerminationReason reason =
          SimulateTrajectory(static_cast<int>(idx), x0, *vy0, config, astro, mu, &records);
      reason_counts[reason] += 1;

      std::cout << "[" << idx + 1 << "/" << x_values.size() << "] x0=" << x0
                << " -> " << ToString(reason) << "\n";
    }

    std::string filename = config.output_basename;
    if (config.append_timestamp) {
      filename += "_" + getcurrent_date();
    }
    if (fs::path(filename).extension().empty()) {
      filename += ".csv";
    }
    const fs::path output_path = output_dir / filename;
    WriteCsv(output_path, config, records);

    std::cout << "Generated " << records.size() << " intersections -> " << output_path << "\n";
    std::cout << "Termination summary:\n";
    for (const auto& entry : reason_counts) {
      std::cout << "  - " << ToString(entry.first) << ": " << entry.second << "\n";
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Fatal error: " << e.what() << "\n";
    return 1;
  }
}
