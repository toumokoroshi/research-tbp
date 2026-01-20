/**
 * @file test_halo_continuation.cpp
 * @brief Halo軌道族の継続テスト（分岐点から生成）
 * @details L2 Lyapunov軌道族のピッチフォーク分岐点からHalo軌道族を生成
 */

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "continuation.hpp"
#include "periodic_orbit_stability.hpp"
#include "rtbp.hpp"

using namespace continuation;
using namespace periodic_orbit;
using namespace lyapunov;

namespace fs = std::filesystem;

// 太陽-地球系の質量パラメータ
constexpr double kMuSunEarth = 3.003480593e-06;

/**
 * @brief CSVファイルからLyapunov軌道族を読み込む
 * @param filepath CSVファイルパス
 * @return 軌道のベクター（index順）
 */
std::vector<PeriodicOrbit<double>> LoadLyapunovFamily(const std::string& filepath) {
  std::vector<PeriodicOrbit<double>> orbits;
  std::ifstream file(filepath);
  
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filepath);
  }
  
  std::string line;
  while (std::getline(file, line)) {
    // コメント行をスキップ
    if (line.empty() || line[0] == '#') continue;
    
    std::istringstream iss(line);
    std::string token;
    std::vector<double> values;
    
    while (std::getline(iss, token, ',')) {
      try {
        values.push_back(std::stod(token));
      } catch (...) {
        break;
      }
    }
    
    // index,x0,y0,z0,vx0,vy0,vz0,period,jacobi,stable,stability_index
    if (values.size() >= 9) {
      PeriodicOrbit<double> orbit;
      orbit.initial_state.x = values[1];
      orbit.initial_state.y = values[2];
      orbit.initial_state.z = values[3];
      orbit.initial_state.vx = values[4];
      orbit.initial_state.vy = values[5];
      orbit.initial_state.vz = values[6];
      orbit.period = values[7];
      orbit.jacobi_constant = values[8];
      if (values.size() >= 10) {
        orbit.is_stable = (values[9] != 0);
      }
      if (values.size() >= 11) {
        orbit.stability_index = values[10];
      }
      orbits.push_back(orbit);
    }
  }
  
  return orbits;
}

/**
 * @brief Halo軌道族をCSVファイルに出力
 */
void SaveHaloFamily(const PeriodicOrbitFamily<double>& family, const std::string& filepath) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filepath);
  }
  
  file << std::setprecision(16);
  file << "# " << family.family_name << "\n";
  file << "# index,x0,y0,z0,vx0,vy0,vz0,period,jacobi,stable,stability_index\n";
  
  for (size_t i = 0; i < family.orbits.size(); ++i) {
    const auto& orb = family.orbits[i];
    file << i << ","
         << orb.initial_state.x << ","
         << orb.initial_state.y << ","
         << orb.initial_state.z << ","
         << orb.initial_state.vx << ","
         << orb.initial_state.vy << ","
         << orb.initial_state.vz << ","
         << orb.period << ","
         << orb.jacobi_constant << ","
         << (orb.is_stable ? 1 : 0) << ","
         << orb.stability_index << "\n";
  }
  
  file.close();
  std::cout << "  Saved: " << filepath << " (" << family.orbits.size() << " orbits)\n";
}

int main() {
  std::cout << std::setprecision(12);
  std::cout << "\n========================================\n";
  std::cout << "  Halo Family Continuation from Bifurcation\n";
  std::cout << "========================================\n";
  std::cout << "  mu = " << kMuSunEarth << " (Sun-Earth system)\n";

  // ============================================
  // 1. Lyapunov軌道族を読み込み
  // ============================================
  std::cout << "\n[Loading Lyapunov Family]\n";
  
  std::string lyapunov_file = "data/lyapunov_continuation/l2_lyapunov_family.csv";
  std::vector<PeriodicOrbit<double>> lyapunov_orbits;
  
  try {
    lyapunov_orbits = LoadLyapunovFamily(lyapunov_file);
    std::cout << "  Loaded " << lyapunov_orbits.size() << " orbits from " << lyapunov_file << "\n";
  } catch (const std::exception& e) {
    std::cerr << "  Error loading Lyapunov family: " << e.what() << std::endl;
    return 1;
  }

  // ============================================
  // 2. 分岐点軌道を取得（orbit #116 - より大きな振幅で安定）
  // ============================================
  // Note: orbit #2 (C=3.0009) は分岐点に近すぎてNewtonが不安定
  //       orbit #116 (C=3.0003) の方がHalo/Lyapunov分離が明確
  constexpr int kBifurcationOrbitIndex = 116;
  
  if (lyapunov_orbits.size() <= kBifurcationOrbitIndex) {
    std::cerr << "  Error: Not enough orbits in family (need at least " 
              << (kBifurcationOrbitIndex + 1) << ")\n";
    return 1;
  }
  
  PeriodicOrbit<double> bifurcation_orbit = lyapunov_orbits[kBifurcationOrbitIndex];
  
  std::cout << "\n[Bifurcation Orbit #" << kBifurcationOrbitIndex << "]\n";
  std::cout << "  x0  = " << bifurcation_orbit.initial_state.x << "\n";
  std::cout << "  vy0 = " << bifurcation_orbit.initial_state.vy << "\n";
  std::cout << "  T   = " << bifurcation_orbit.period << "\n";
  std::cout << "  C   = " << bifurcation_orbit.jacobi_constant << "\n";

  // ============================================
  // 3. モノドロミー行列を計算
  // ============================================
  std::cout << "\n[Computing Monodromy Matrix]\n";
  
  double integration_dt = 0.001;
  bifurcation_orbit.monodromy_matrix = 
      ComputeMonodromyMatrix(bifurcation_orbit, kMuSunEarth, integration_dt);
  AnalyzeStability(&bifurcation_orbit, 1.5);
  
  std::cout << "  Stability index: " << bifurcation_orbit.stability_index << "\n";
  std::cout << "  Is stable: " << (bifurcation_orbit.is_stable ? "yes" : "no") << "\n";
  std::cout << "  Eigenvalues:\n";
  for (const auto& ev : bifurcation_orbit.eigenvalues) {
    std::cout << "    " << ev.real() << " + " << ev.imag() << "i\n";
  }

  // ============================================
  // 4. 継続法の設定
  // ============================================
  ContinuationConfig<double> config;
  config.method = ContinuationMethod::PSEUDO_ARCLENGTH;
  config.initial_step_size = 0.0005;
  config.min_step_size = 1e-7;
  config.max_step_size = 0.01;
  config.max_steps = 300;  // Halo軌道族は短めに
  
  config.jacobi_min = 2.85;
  config.jacobi_max = 3.02;
  
  config.newton_max_iterations = 50;
  config.newton_tolerance = 1e-9;
  config.max_integration_time = 15.0;  // 周期の約5倍
  config.integration_timestep = 0.001;
  
  config.detect_bifurcations = true;
  config.eigenvalue_threshold = 1e-4;
  
  config.section_index = 1;  // y = 0
  config.section_value = 0.0;

  OrbitContinuator<double> continuator(kMuSunEarth, config);

  // ============================================
  // 5. Halo族の生成（North）
  // ============================================
  std::cout << "\n========================================\n";
  std::cout << "  Generating North Halo Family\n";
  std::cout << "========================================\n";
  
  PeriodicOrbitFamily<double> halo_north;
  try {
    halo_north = continuator.ContinueHaloFromBifurcation(bifurcation_orbit, true);
    std::cout << "\n  North Halo family generated: " << halo_north.orbits.size() << " orbits\n";
    std::cout << "  Termination: " << halo_north.termination_reason << "\n";
  } catch (const std::exception& e) {
    std::cerr << "  Error generating North Halo: " << e.what() << std::endl;
  }

  // ============================================
  // 6. Halo族の生成（South）
  // ============================================
  std::cout << "\n========================================\n";
  std::cout << "  Generating South Halo Family\n";
  std::cout << "========================================\n";
  
  PeriodicOrbitFamily<double> halo_south;
  try {
    halo_south = continuator.ContinueHaloFromBifurcation(bifurcation_orbit, false);
    std::cout << "\n  South Halo family generated: " << halo_south.orbits.size() << " orbits\n";
    std::cout << "  Termination: " << halo_south.termination_reason << "\n";
  } catch (const std::exception& e) {
    std::cerr << "  Error generating South Halo: " << e.what() << std::endl;
  }

  // ============================================
  // 7. 結果の出力
  // ============================================
  std::string output_dir = "data/halo_continuation";
  fs::create_directories(output_dir);
  
  std::cout << "\n[Saving Results]\n";
  
  if (!halo_north.orbits.empty()) {
    SaveHaloFamily(halo_north, output_dir + "/l2_halo_north.csv");
    
    // 最初と最後の軌道情報を表示
    const auto& first = halo_north.orbits.front();
    const auto& last = halo_north.orbits.back();
    std::cout << "  North Halo range:\n";
    std::cout << "    z: " << first.initial_state.z << " -> " << last.initial_state.z << "\n";
    std::cout << "    C: " << first.jacobi_constant << " -> " << last.jacobi_constant << "\n";
    std::cout << "    T: " << first.period << " -> " << last.period << "\n";
  }
  
  if (!halo_south.orbits.empty()) {
    SaveHaloFamily(halo_south, output_dir + "/l2_halo_south.csv");
    
    const auto& first = halo_south.orbits.front();
    const auto& last = halo_south.orbits.back();
    std::cout << "  South Halo range:\n";
    std::cout << "    z: " << first.initial_state.z << " -> " << last.initial_state.z << "\n";
    std::cout << "    C: " << first.jacobi_constant << " -> " << last.jacobi_constant << "\n";
    std::cout << "    T: " << first.period << " -> " << last.period << "\n";
  }

  // ============================================
  // 8. テスト結果サマリー
  // ============================================
  std::cout << "\n========================================\n";
  std::cout << "  Test Summary\n";
  std::cout << "========================================\n";
  
  bool north_success = !halo_north.orbits.empty() && halo_north.orbits.size() > 1;
  bool south_success = !halo_south.orbits.empty() && halo_south.orbits.size() > 1;
  
  // z座標の符号チェック
  bool north_z_positive = north_success && halo_north.orbits.back().initial_state.z > 0;
  bool south_z_negative = south_success && halo_south.orbits.back().initial_state.z < 0;
  
  std::cout << "  North Halo: " << (north_success ? "[PASS]" : "[FAIL]") 
            << " (" << halo_north.orbits.size() << " orbits)\n";
  if (north_success) {
    std::cout << "    z > 0: " << (north_z_positive ? "[OK]" : "[WARN]") << "\n";
  }
  
  std::cout << "  South Halo: " << (south_success ? "[PASS]" : "[FAIL]") 
            << " (" << halo_south.orbits.size() << " orbits)\n";
  if (south_success) {
    std::cout << "    z < 0: " << (south_z_negative ? "[OK]" : "[WARN]") << "\n";
  }
  
  std::cout << "========================================\n";
  
  // 少なくとも1つのHalo族が生成されれば成功
  return (north_success || south_success) ? 0 : 1;
}
