/**
 * @file main_english.cpp
 * @brief Periodic Orbit Analysis Application (English version)
 * @details Detect, analyze stability, and compute invariant manifolds for periodic orbits in CR3BP
 * @date 2025-12-24
 */

#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "periodic_orbit.hpp"
#include "rtbp.hpp"

// Use toml++
#include <toml++/toml.hpp>

using namespace my_type;
using namespace periodic_orbit;

// ---------------------------------------------------------------------------
// Config Structure
// ---------------------------------------------------------------------------

struct PeriodicOrbitConfig {
  // Mode settings
  std::string execution_mode = "all";  // "detect", "stability", "manifold", "all"

  // Newton-Raphson settings
  int newton_max_iterations = 50;
  double newton_tolerance = 1e-12;
  std::string initial_guess_source = "file";  // "file" or "poincare_map"
  std::string initial_guess_file;
  std::string poincare_map_output_dir;

  // Poincare section settings
  std::string section_var = "y";
  double section_value = 0.0;
  std::string crossing_direction = "positive";
  double max_integration_time = 500.0;

  // Stability analysis settings
  bool compute_eigenvectors = true;
  double eigenvalue_threshold = 1.0;
  double stm_epsilon = 1e-8;

  // Manifold settings
  ManifoldConfig<double> manifold_config;

  // Integrator settings
  std::string integrator_type = "rk45";
  double timestep = 0.0001;
  double rk45_atol = 1e-12;
  double rk45_rtol = 1e-12;
  double rk45_initial_step = 0.0001;
  double rk45_max_step = 0.01;

  // System parameters
  double mu = 3.0404233870218e-06;
  double jacobi_constant = 3.0007;

  // Output settings
  std::string output_dir = "data/periodic_orbit";
  std::string tag = "test_orbit";
  bool save_trajectory = true;
  bool save_monodromy_matrix = true;

  // Collision detection
  double collision_radius_earth = 1e-12;
  double collision_radius_sun = 1e-12;
};

// ---------------------------------------------------------------------------
// Function Declarations
// ---------------------------------------------------------------------------

bool LoadConfig(const std::string& filepath, PeriodicOrbitConfig* config);
void PrintConfig(const PeriodicOrbitConfig& config);
std::string CreateOutputDirectory(const PeriodicOrbitConfig& config);
std::vector<State<double>> LoadInitialGuesses(const PeriodicOrbitConfig& config);
void SavePeriodicOrbit(const PeriodicOrbit<double>& orbit, const std::string& output_dir,
                       const PeriodicOrbitConfig& config);
void SaveManifolds(const std::vector<ManifoldTrajectory<double>>& manifolds,
                   const std::string& output_dir);

// ---------------------------------------------------------------------------
// Load Configuration File
// ---------------------------------------------------------------------------

bool LoadConfig(const std::string& filepath, PeriodicOrbitConfig* config) {
  try {
    auto data = toml::parse_file(filepath);

    // Mode settings
    if (auto mode = data["mode"]["execution_mode"].value<std::string>()) {
      config->execution_mode = *mode;
    }

    // Newton-Raphson settings
    if (auto val = data["newton_raphson"]["max_iterations"].value<int>()) {
      config->newton_max_iterations = *val;
    }
    if (auto val = data["newton_raphson"]["tolerance"].value<double>()) {
      config->newton_tolerance = *val;
    }
    if (auto val = data["newton_raphson"]["initial_guess_source"].value<std::string>()) {
      config->initial_guess_source = *val;
    }
    if (auto val = data["newton_raphson"]["initial_guess_file"].value<std::string>()) {
      config->initial_guess_file = *val;
    }
    if (auto val = data["newton_raphson"]["poincare_map_output_dir"].value<std::string>()) {
      config->poincare_map_output_dir = *val;
    }

    // Poincare section settings
    if (auto val = data["poincare_section"]["section_var"].value<std::string>()) {
      config->section_var = *val;
    }
    if (auto val = data["poincare_section"]["section_value"].value<double>()) {
      config->section_value = *val;
    }
    if (auto val = data["poincare_section"]["crossing_direction"].value<std::string>()) {
      config->crossing_direction = *val;
    }
    if (auto val = data["poincare_section"]["max_integration_time"].value<double>()) {
      config->max_integration_time = *val;
    }

    // Stability analysis settings
    if (auto val = data["stability"]["compute_eigenvectors"].value<bool>()) {
      config->compute_eigenvectors = *val;
    }
    if (auto val = data["stability"]["eigenvalue_threshold"].value<double>()) {
      config->eigenvalue_threshold = *val;
    }
    if (auto val = data["stability"]["stm_finite_diff_epsilon"].value<double>()) {
      config->stm_epsilon = *val;
    }

    // Manifold settings
    if (auto val = data["manifold"]["epsilon"].value<double>()) {
      config->manifold_config.epsilon = *val;
    }
    if (auto val = data["manifold"]["forward_time"].value<double>()) {
      config->manifold_config.forward_time = *val;
    }
    if (auto val = data["manifold"]["backward_time"].value<double>()) {
      config->manifold_config.backward_time = *val;
    }
    if (auto val = data["manifold"]["num_initial_points"].value<int>()) {
      config->manifold_config.num_initial_points = *val;
    }
    if (auto val = data["manifold"]["compute_stable"].value<bool>()) {
      config->manifold_config.compute_stable = *val;
    }
    if (auto val = data["manifold"]["compute_unstable"].value<bool>()) {
      config->manifold_config.compute_unstable = *val;
    }

    // Integrator settings
    if (auto val = data["integrator"]["type"].value<std::string>()) {
      config->integrator_type = *val;
    }
    if (auto val = data["integrator"]["timestep"].value<double>()) {
      config->timestep = *val;
    }
    if (auto val = data["integrator"]["rk45_atol"].value<double>()) {
      config->rk45_atol = *val;
    }
    if (auto val = data["integrator"]["rk45_rtol"].value<double>()) {
      config->rk45_rtol = *val;
    }

    // System parameters
    if (auto val = data["system"]["mu"].value<double>()) {
      config->mu = *val;
    }
    if (auto val = data["system"]["jacobi_constant"].value<double>()) {
      config->jacobi_constant = *val;
    }

    // Output settings
    if (auto val = data["output"]["output_dir"].value<std::string>()) {
      config->output_dir = *val;
    }
    if (auto val = data["output"]["tag"].value<std::string>()) {
      config->tag = *val;
    }
    if (auto val = data["output"]["save_trajectory"].value<bool>()) {
      config->save_trajectory = *val;
    }
    if (auto val = data["output"]["save_monodromy_matrix"].value<bool>()) {
      config->save_monodromy_matrix = *val;
    }

    // Collision detection
    if (auto val = data["collision"]["collision_radius_earth"].value<double>()) {
      config->collision_radius_earth = *val;
    }
    if (auto val = data["collision"]["collision_radius_sun"].value<double>()) {
      config->collision_radius_sun = *val;
    }

    return true;
  } catch (const toml::parse_error& err) {
    std::cerr << "Config parse error: " << err << std::endl;
    return false;
  } catch (const std::exception& e) {
    std::cerr << "Config load error: " << e.what() << std::endl;
    return false;
  }
}

// ---------------------------------------------------------------------------
// Print Configuration
// ---------------------------------------------------------------------------

void PrintConfig(const PeriodicOrbitConfig& config) {
  std::cout << "\n========================================\n";
  std::cout << "  Configuration Summary\n";
  std::cout << "========================================\n\n";

  std::cout << "[Execution Mode]\n";
  std::cout << "  Mode: " << config.execution_mode << "\n\n";

  std::cout << "[Newton-Raphson Settings]\n";
  std::cout << "  Max iterations: " << config.newton_max_iterations << "\n";
  std::cout << "  Tolerance: " << config.newton_tolerance << "\n";
  std::cout << "  Initial guess source: " << config.initial_guess_source << "\n";
  if (config.initial_guess_source == "file") {
    std::cout << "  Initial guess file: " << config.initial_guess_file << "\n";
  }
  std::cout << "\n";

  std::cout << "[Poincare Section]\n";
  std::cout << "  Section variable: " << config.section_var << " = " << config.section_value
            << "\n";
  std::cout << "  Crossing direction: " << config.crossing_direction << "\n\n";

  std::cout << "[Stability Analysis]\n";
  std::cout << "  Compute eigenvectors: " << (config.compute_eigenvectors ? "enabled" : "disabled")
            << "\n";
  std::cout << "  Eigenvalue threshold: " << config.eigenvalue_threshold << "\n\n";

  std::cout << "[Invariant Manifolds]\n";
  std::cout << "  Initial displacement: " << config.manifold_config.epsilon << "\n";
  std::cout << "  Forward time: " << config.manifold_config.forward_time << "\n";
  std::cout << "  Backward time: " << config.manifold_config.backward_time << "\n";
  std::cout << "  Number of initial points: " << config.manifold_config.num_initial_points
            << "\n\n";

  std::cout << "[Integrator]\n";
  std::cout << "  Type: " << config.integrator_type << "\n";
  std::cout << "  Timestep: " << config.timestep << "\n\n";

  std::cout << "[System Parameters]\n";
  std::cout << "  Mass ratio mu: " << std::scientific << config.mu << "\n";
  std::cout << "  Jacobi constant: " << config.jacobi_constant << "\n\n";

  std::cout << "[Output Settings]\n";
  std::cout << "  Output directory: " << config.output_dir << "\n";
  std::cout << "  Tag: " << config.tag << "\n";
  std::cout << "========================================\n\n";
}

// ---------------------------------------------------------------------------
// Create Output Directory
// ---------------------------------------------------------------------------

std::string CreateOutputDirectory(const PeriodicOrbitConfig& config) {
  namespace fs = std::filesystem;

  // Generate timestamp
  auto now = std::chrono::system_clock::now();
  auto time_t = std::chrono::system_clock::to_time_t(now);
  std::tm tm;
  localtime_s(&tm, &time_t);

  std::ostringstream oss;
  oss << std::put_time(&tm, "%y_%m%d_%H%M");
  std::string timestamp = oss.str();

  // Output directory path
  fs::path output_path = fs::path(config.output_dir) / timestamp / config.tag;

  // Create directory
  fs::create_directories(output_path);

  std::cout << "Output directory created: " << output_path.string() << "\n\n";

  return output_path.string();
}

// ---------------------------------------------------------------------------
// Load Initial Guesses
// ---------------------------------------------------------------------------

std::vector<State<double>> LoadInitialGuesses(const PeriodicOrbitConfig& config) {
  std::vector<State<double>> guesses;

  if (config.initial_guess_source == "file") {
    std::ifstream file(config.initial_guess_file);
    if (!file.is_open()) {
      std::cerr << "Cannot open initial guess file: " << config.initial_guess_file << std::endl;
      return guesses;
    }

    std::string line;
    while (std::getline(file, line)) {
      // Skip comments and empty lines
      if (line.empty() || line[0] == '#') continue;

      std::istringstream iss(line);
      State<double> state;
      char comma;
      if (iss >> state.x >> comma >> state.y >> comma >> state.z >> comma >> state.vx >> comma >>
          state.vy >> comma >> state.vz) {
        guesses.push_back(state);
      }
    }

    std::cout << "Loaded initial guesses: " << guesses.size() << " points\n";
  } else if (config.initial_guess_source == "poincare_map") {
    // TODO: Load from PoincareMap output
    std::cerr << "Loading from PoincareMap output not implemented yet\n";
  }

  return guesses;
}

// ---------------------------------------------------------------------------
// Save Periodic Orbit
// ---------------------------------------------------------------------------

void SavePeriodicOrbit(const PeriodicOrbit<double>& orbit, const std::string& output_dir,
                       const PeriodicOrbitConfig& config) {
  namespace fs = std::filesystem;

  // Save basic information
  std::ofstream info_file(fs::path(output_dir) / "orbit_info.txt");
  info_file << std::setprecision(16);
  info_file << "Periodic Orbit Information\n";
  info_file << "==========================\n\n";
  info_file << "Initial State:\n";
  info_file << "  x  = " << orbit.initial_state.x << "\n";
  info_file << "  y  = " << orbit.initial_state.y << "\n";
  info_file << "  z  = " << orbit.initial_state.z << "\n";
  info_file << "  vx = " << orbit.initial_state.vx << "\n";
  info_file << "  vy = " << orbit.initial_state.vy << "\n";
  info_file << "  vz = " << orbit.initial_state.vz << "\n\n";
  info_file << "Period: " << orbit.period << "\n";
  info_file << "Jacobi Constant: " << orbit.jacobi_constant << "\n\n";

  if (orbit.stability_computed) {
    info_file << "Stability Analysis:\n";
    info_file << "  Stable: " << (orbit.is_stable ? "Yes (SPO)" : "No (UPO)") << "\n";
    info_file << "  Stability Index: " << orbit.stability_index << "\n";
    info_file << "  Eigenvalues:\n";
    for (size_t i = 0; i < orbit.eigenvalues.size(); ++i) {
      info_file << "    lambda" << i + 1 << " = " << orbit.eigenvalues[i]
                << " (|lambda| = " << std::abs(orbit.eigenvalues[i]) << ")\n";
    }
  }

  info_file.close();
  std::cout << "Saved orbit info: orbit_info.txt\n";

  // Save monodromy matrix
  if (config.save_monodromy_matrix && orbit.stability_computed) {
    std::ofstream matrix_file(fs::path(output_dir) / "monodromy_matrix.csv");
    matrix_file << std::setprecision(16);
    for (const auto& row : orbit.monodromy_matrix) {
      for (size_t j = 0; j < 6; ++j) {
        matrix_file << row[j];
        if (j < 5) matrix_file << ",";
      }
      matrix_file << "\n";
    }
    matrix_file.close();
    std::cout << "Saved monodromy matrix: monodromy_matrix.csv\n";
  }
}

// ---------------------------------------------------------------------------
// Save Manifold Data
// ---------------------------------------------------------------------------

void SaveManifolds(const std::vector<ManifoldTrajectory<double>>& manifolds,
                   const std::string& output_dir) {
  namespace fs = std::filesystem;

  for (size_t i = 0; i < manifolds.size(); ++i) {
    const auto& manifold = manifolds[i];
    std::string type_str =
        (manifold.type == ManifoldTrajectory<double>::Type::STABLE) ? "stable" : "unstable";

    std::ostringstream filename;
    filename << "manifold_" << type_str << "_" << std::setfill('0') << std::setw(3) << i << ".csv";

    std::ofstream file(fs::path(output_dir) / filename.str());
    file << std::setprecision(16);
    file << "# Manifold Type: " << type_str << "\n";
    file << "# time,x,y,z,vx,vy,vz\n";

    for (size_t j = 0; j < manifold.trajectory.size(); ++j) {
      const auto& state = manifold.trajectory[j];
      file << manifold.times[j] << "," << state.x << "," << state.y << "," << state.z << ","
           << state.vx << "," << state.vy << "," << state.vz << "\n";
    }
    file.close();
  }

  std::cout << "Saved manifold data: " << manifolds.size() << " files\n";
}

// ---------------------------------------------------------------------------
// Main Function
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  std::cout << "========================================\n";
  std::cout << "  Periodic Orbit Analysis\n";
  std::cout << "  CR3BP Periodic Orbit Detection\n";
  std::cout << "========================================\n\n";

  // Check command line arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <config_file.toml>\n";
    return 1;
  }

  std::string config_file = argv[1];

  // Load configuration file
  PeriodicOrbitConfig config;
  if (!LoadConfig(config_file, &config)) {
    std::cerr << "Failed to load configuration file\n";
    return 1;
  }

  PrintConfig(config);

  // Create output directory
  std::string output_dir = CreateOutputDirectory(config);

  // Load initial guesses
  auto guesses = LoadInitialGuesses(config);
  if (guesses.empty()) {
    std::cerr << "No initial guesses loaded\n";
    return 1;
  }

  std::cout << "\n========================================\n";
  std::cout << "  Starting Periodic Orbit Detection\n";
  std::cout << "========================================\n\n";

  // Get section index
  int section_index = -1;
  if (config.section_var == "x")
    section_index = 0;
  else if (config.section_var == "y")
    section_index = 1;
  else if (config.section_var == "z")
    section_index = 2;
  else if (config.section_var == "vx")
    section_index = 3;
  else if (config.section_var == "vy")
    section_index = 4;
  else if (config.section_var == "vz")
    section_index = 5;

  if (section_index < 0) {
    std::cerr << "Invalid section variable: " << config.section_var << "\n";
    return 1;
  }

  // Detect periodic orbits for each initial guess
  std::vector<PeriodicOrbit<double>> detected_orbits;

  for (size_t i = 0; i < guesses.size(); ++i) {
    std::cout << "\nInitial guess " << (i + 1) << "/" << guesses.size() << ":\n";
    std::cout << "  (" << guesses[i].x << ", " << guesses[i].y << ", " << guesses[i].z << ", "
              << guesses[i].vx << ", " << guesses[i].vy << ", " << guesses[i].vz << ")\n";

    try {
      NewtonConvergenceInfo<double> conv_info;
      PeriodicOrbit<double> orbit = RefinePeriodicOrbit(
          guesses[i], config.mu, section_index, config.section_value, config.newton_max_iterations,
          config.newton_tolerance, config.max_integration_time, config.timestep, &conv_info);

      std::cout << "  Success! Converged\n";
      std::cout << "    Iterations: " << conv_info.iterations << "\n";
      std::cout << "    Final residual: " << std::scientific << conv_info.final_residual << "\n";
      std::cout << "    Period: " << std::fixed << std::setprecision(6) << orbit.period << "\n";
      std::cout << "    Jacobi constant: " << orbit.jacobi_constant << "\n";

      // Stability analysis
      if (config.execution_mode == "all" || config.execution_mode == "stability") {
        std::cout << "    Computing stability...\n";

        try {
          // Compute monodromy matrix
          orbit.monodromy_matrix = ComputeMonodromyMatrix(orbit, config.mu, config.timestep);

          // Analyze stability
          AnalyzeStability(&orbit, config.eigenvalue_threshold);

          std::cout << "    Stability: " << (orbit.is_stable ? "SPO (Stable)" : "UPO (Unstable)")
                    << "\n";
          std::cout << "    Stability index: " << orbit.stability_index << "\n";

        } catch (const std::exception& e) {
          std::cerr << "    Stability analysis failed: " << e.what() << "\n";
        }
      }

      detected_orbits.push_back(orbit);

    } catch (const std::exception& e) {
      std::cerr << "  Failed: " << e.what() << "\n";
    }
  }

  std::cout << "\n========================================\n";
  std::cout << "  Detected periodic orbits: " << detected_orbits.size() << "\n";
  std::cout << "========================================\n\n";

  if (detected_orbits.empty()) {
    std::cout << "No periodic orbits detected\n";
    return 0;
  }

  // Save detected orbits
  std::cout << "Saving periodic orbit information...\n";
  for (size_t i = 0; i < detected_orbits.size(); ++i) {
    std::string orbit_dir = output_dir + "/orbit_" + std::to_string(i + 1);
    namespace fs = std::filesystem;
    fs::create_directories(orbit_dir);
    SavePeriodicOrbit(detected_orbits[i], orbit_dir, config);
  }

  std::cout << "\nProcessing complete\n";
  return 0;
}
