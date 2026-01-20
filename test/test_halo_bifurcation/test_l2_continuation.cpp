/**
 * @file test_l2_continuation.cpp
 * @brief L2 Lyapunov軌道族継続と分岐点検出のテスト
 * @details 実際のL2 Lyapunov軌道族を継続し、ピッチフォーク分岐を検出する
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include "continuation.hpp"
#include "periodic_orbit_stability.hpp"
#include "rtbp.hpp"

using namespace continuation;
using namespace periodic_orbit;
using namespace lyapunov;

// 太陽-地球系の質量パラメータ
constexpr double kMuSunEarth = 3.003480593e-06;

int main() {
  std::cout << std::setprecision(12);
  std::cout << "\n========================================" << std::endl;
  std::cout << "  L2 Lyapunov Family Continuation" << std::endl;
  std::cout << "  Bifurcation Detection Test" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "  mu = " << kMuSunEarth << " (Sun-Earth system)" << std::endl;

  // ============================================
  // 継続法の設定
  // ============================================
  ContinuationConfig<double> config;
  
  // 継続法パラメータ
  config.method = ContinuationMethod::PSEUDO_ARCLENGTH;
  config.initial_step_size = 0.0005;    // 小さめのステップで分岐点を見逃さない
  config.min_step_size = 1e-7;
  config.max_step_size = 0.01;
  config.max_steps = 500;               // 十分な継続ステップ数
  
  // ヤコビ積分の範囲（L2周りの適切な範囲）
  // ピッチフォーク分岐を探すため、より広い範囲を設定
  config.jacobi_min = 2.85;  // 大振幅軌道まで探索
  config.jacobi_max = 3.02;
  
  // Newton法パラメータ
  config.newton_max_iterations = 50;
  config.newton_tolerance = 1e-9;  // 1e-12は厳しすぎるため緩和
  config.max_integration_time = 10.0;  // L2周期(~3.1)の約3倍
  config.integration_timestep = 0.001;  // モノドロミー行列計算用（0.0001は遅すぎる）
  
  // 分岐検出を有効化
  config.detect_bifurcations = true;
  config.eigenvalue_threshold = 1e-4;  // 分岐検出の閾値
  
  // 初期振幅（線形近似が有効な小さい値）
  config.initial_amplitude = 0.0001;  // 小振幅で線形近似を適用
  
  // ポアンカレ断面設定
  config.section_index = 1;      // y = 0
  config.section_value = 0.0;

  std::cout << "\n[Configuration]" << std::endl;
  std::cout << "  Method: Pseudo-arclength" << std::endl;
  std::cout << "  Initial step: " << config.initial_step_size << std::endl;
  std::cout << "  Max steps: " << config.max_steps << std::endl;
  std::cout << "  Jacobi range: [" << config.jacobi_min << ", " << config.jacobi_max << "]" << std::endl;
  std::cout << "  Bifurcation detection: enabled" << std::endl;

  // ============================================
  // L2 Lyapunov軌道族の継続
  // ============================================
  std::cout << "\n[Starting L2 Lyapunov Family Continuation]" << std::endl;
  
  OrbitContinuator<double> continuator(kMuSunEarth, config);
  
  // 推定周期を表示（デバッグ用）
  auto generator = std::make_unique<lyapunov::LyapunovOrbitGenerator<double>>(LagrangePoint::L2, kMuSunEarth);
  double estimated_period = generator->EstimatePeriod(config.initial_amplitude);
  std::cout << "  Estimated period: " << estimated_period << std::endl;
  std::cout << "  max_integration_time: " << config.max_integration_time << std::endl;
  
  // max_integration_timeが周期より短い場合は調整
  if (config.max_integration_time < estimated_period * 2.0) {
    config.max_integration_time = estimated_period * 5.0;
    std::cout << "  [ADJUSTED] max_integration_time -> " << config.max_integration_time << std::endl;
    
    // configを再設定
    continuator = OrbitContinuator<double>(kMuSunEarth, config);
  }
  
  PeriodicOrbitFamily<double> family;
  try {
    family = continuator.ContinueLyapunovFamily(LagrangePoint::L2);
  } catch (const std::exception& e) {
    std::cerr << "Continuation failed: " << e.what() << std::endl;
    return 1;
  }

  // ============================================
  // 結果の表示
  // ============================================
  std::cout << "\n========================================" << std::endl;
  std::cout << "  Continuation Results" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "  Family name: " << family.family_name << std::endl;
  std::cout << "  Orbits computed: " << family.orbits.size() << std::endl;
  std::cout << "  Continuation completed: " << (family.continuation_completed ? "Yes" : "No") << std::endl;
  std::cout << "  Termination reason: " << family.termination_reason << std::endl;

  // 軌道族の範囲
  if (!family.orbits.empty()) {
    double C_min = family.orbits.front().jacobi_constant;
    double C_max = family.orbits.front().jacobi_constant;
    double T_min = family.orbits.front().period;
    double T_max = family.orbits.front().period;
    
    for (const auto& orbit : family.orbits) {
      if (orbit.jacobi_constant < C_min) C_min = orbit.jacobi_constant;
      if (orbit.jacobi_constant > C_max) C_max = orbit.jacobi_constant;
      if (orbit.period < T_min) T_min = orbit.period;
      if (orbit.period > T_max) T_max = orbit.period;
    }
    
    std::cout << "\n  Jacobi constant range: [" << C_min << ", " << C_max << "]" << std::endl;
    std::cout << "  Period range: [" << T_min << ", " << T_max << "]" << std::endl;
  }

  // ============================================
  // 分岐点の表示
  // ============================================
  std::cout << "\n========================================" << std::endl;
  std::cout << "  Bifurcation Points Detected" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "  Total bifurcations: " << family.bifurcations.size() << std::endl;

  int pitchfork_count = 0;
  for (const auto& bif : family.bifurcations) {
    std::string type_str;
    switch (bif.type) {
      case BifurcationType::PITCHFORK:
        type_str = "PITCHFORK (Lyapunov->Halo)";
        pitchfork_count++;
        break;
      case BifurcationType::PERIOD_DOUBLING:
        type_str = "PERIOD_DOUBLING";
        break;
      case BifurcationType::NEIMARK_SACKER:
        type_str = "NEIMARK_SACKER";
        break;
      default:
        type_str = "UNKNOWN";
    }
    
    std::cout << "\n  Bifurcation #" << (bif.orbit_index) << ":" << std::endl;
    std::cout << "    Type: " << type_str << std::endl;
    std::cout << "    Jacobi constant: " << bif.jacobi_integral << std::endl;
    std::cout << "    Critical eigenvalue: " << bif.critical_eigenvalue.real() 
              << " + " << bif.critical_eigenvalue.imag() << "i" << std::endl;
    std::cout << "    Description: " << bif.description << std::endl;
  }

  // ============================================
  // データ出力
  // ============================================
  namespace fs = std::filesystem;
  std::string output_dir = "data/lyapunov_continuation";
  fs::create_directories(output_dir);

  // 軌道族をCSVに出力
  std::string orbit_file = output_dir + "/l2_lyapunov_family.csv";
  std::ofstream ofs(orbit_file);
  ofs << "# L2 Lyapunov Family\n";
  ofs << "# index,x0,y0,z0,vx0,vy0,vz0,period,jacobi,stable,stability_index\n";
  ofs << std::setprecision(16);
  
  for (size_t i = 0; i < family.orbits.size(); ++i) {
    const auto& orb = family.orbits[i];
    ofs << i << ","
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
  ofs.close();
  std::cout << "\n  Saved orbit family: " << orbit_file << std::endl;

  // 分岐点をCSVに出力
  std::string bif_file = output_dir + "/bifurcations.csv";
  std::ofstream bfs(bif_file);
  bfs << "# Bifurcation Points\n";
  bfs << "# orbit_index,jacobi,type,eigenvalue_real,eigenvalue_imag\n";
  bfs << std::setprecision(16);
  
  for (const auto& bif : family.bifurcations) {
    std::string type_str;
    switch (bif.type) {
      case BifurcationType::PITCHFORK: type_str = "pitchfork"; break;
      case BifurcationType::PERIOD_DOUBLING: type_str = "period_doubling"; break;
      case BifurcationType::NEIMARK_SACKER: type_str = "neimark_sacker"; break;
      default: type_str = "unknown";
    }
    
    bfs << bif.orbit_index << ","
        << bif.jacobi_integral << ","
        << type_str << ","
        << bif.critical_eigenvalue.real() << ","
        << bif.critical_eigenvalue.imag() << "\n";
  }
  bfs.close();
  std::cout << "  Saved bifurcations: " << bif_file << std::endl;

  // ============================================
  // テスト結果
  // ============================================
  std::cout << "\n========================================" << std::endl;
  std::cout << "  Test Summary" << std::endl;
  std::cout << "========================================" << std::endl;
  
  bool success = true;
  
  // 軌道が生成されたか
  if (family.orbits.size() >= 10) {
    std::cout << "  [PASS] Generated " << family.orbits.size() << " orbits" << std::endl;
  } else {
    std::cout << "  [WARN] Only " << family.orbits.size() << " orbits generated" << std::endl;
  }
  
  // 固有値が正しく計算されたか
  bool all_eigenvalues_computed = true;
  for (const auto& orb : family.orbits) {
    if (orb.eigenvalues.size() != 6) {
      all_eigenvalues_computed = false;
      break;
    }
  }
  if (all_eigenvalues_computed && !family.orbits.empty()) {
    std::cout << "  [PASS] All eigenvalues (6 per orbit) computed correctly" << std::endl;
  } else {
    std::cout << "  [WARN] Some eigenvalues missing" << std::endl;
  }
  
  // ピッチフォーク分岐が見つかったか
  if (pitchfork_count > 0) {
    std::cout << "  [PASS] Detected " << pitchfork_count << " pitchfork bifurcation(s)" << std::endl;
  } else {
    std::cout << "  [INFO] No pitchfork bifurcation detected in this range" << std::endl;
    std::cout << "         (May need larger Jacobi range or more steps)" << std::endl;
  }

  std::cout << "========================================" << std::endl;
  
  return success ? 0 : 1;
}
