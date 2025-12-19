/**
 * @brief TOML設定ファイルのパース検証テスト
 *
 * オリジナルの.txtファイルとTOML形式のファイルから同じ情報を
 * 取得できることを確認する
 */
#include <cmath>
#include <iomanip>
#include <iostream>

#include "utils.hpp"


int main() {
  std::cout << "=== TOML Config Parser Verification Test ===" << std::endl;

  bool all_passed = true;

  try {
    // TOML設定ファイルをパース
    utils::TomlConfigParser toml_config(
        "C:/WORKSPACE/Research_git/TBP/configs/trajectory_calc/trajectory_calc_config_1.toml");

    // 期待値 (オリジナルの.txtファイルから)
    const double expected_timestep = 0.000001;
    const double expected_threshold = 6.28;
    const std::string expected_chaos = "NONE";
    const std::string expected_integrator = "SYMPLECTIC6";
    const size_t expected_coords_count = 4;

    // 検証
    std::cout << std::fixed << std::setprecision(10);

    // calc_timestep
    double actual_timestep = toml_config.GetDouble("simulation.calc_timestep", 0.0);
    std::cout << "\n[1] simulation.calc_timestep" << std::endl;
    std::cout << "    Expected: " << expected_timestep << std::endl;
    std::cout << "    Actual:   " << actual_timestep << std::endl;
    if (std::abs(actual_timestep - expected_timestep) < 1e-12) {
      std::cout << "    => PASS" << std::endl;
    } else {
      std::cout << "    => FAIL" << std::endl;
      all_passed = false;
    }

    // time_threshold
    double actual_threshold = toml_config.GetDouble("simulation.time_threshold", 0.0);
    std::cout << "\n[2] simulation.time_threshold" << std::endl;
    std::cout << "    Expected: " << expected_threshold << std::endl;
    std::cout << "    Actual:   " << actual_threshold << std::endl;
    if (std::abs(actual_threshold - expected_threshold) < 1e-12) {
      std::cout << "    => PASS" << std::endl;
    } else {
      std::cout << "    => FAIL" << std::endl;
      all_passed = false;
    }

    // chaos.index_type
    std::string actual_chaos = toml_config.GetString("chaos.index_type", "");
    std::cout << "\n[3] chaos.index_type" << std::endl;
    std::cout << "    Expected: " << expected_chaos << std::endl;
    std::cout << "    Actual:   " << actual_chaos << std::endl;
    if (actual_chaos == expected_chaos) {
      std::cout << "    => PASS" << std::endl;
    } else {
      std::cout << "    => FAIL" << std::endl;
      all_passed = false;
    }

    // integrator.type
    std::string actual_integrator = toml_config.GetString("integrator.type", "");
    std::cout << "\n[4] integrator.type" << std::endl;
    std::cout << "    Expected: " << expected_integrator << std::endl;
    std::cout << "    Actual:   " << actual_integrator << std::endl;
    if (actual_integrator == expected_integrator) {
      std::cout << "    => PASS" << std::endl;
    } else {
      std::cout << "    => FAIL" << std::endl;
      all_passed = false;
    }

    // coords array
    auto coords = toml_config.GetCoordsArray("coords");
    std::cout << "\n[5] coords array count" << std::endl;
    std::cout << "    Expected: " << expected_coords_count << std::endl;
    std::cout << "    Actual:   " << coords.size() << std::endl;
    if (coords.size() == expected_coords_count) {
      std::cout << "    => PASS" << std::endl;
    } else {
      std::cout << "    => FAIL" << std::endl;
      all_passed = false;
    }

    // 最初の座標を詳細検証
    if (!coords.empty()) {
      std::cout << "\n[6] First coord values" << std::endl;
      std::cout << "    x:  " << coords[0].x << " (expected: 0.989293478931466)" << std::endl;
      std::cout << "    y:  " << coords[0].y << " (expected: 0.000777777777778)" << std::endl;
      std::cout << "    z:  " << coords[0].z << " (expected: 0.0)" << std::endl;
      std::cout << "    vx: " << coords[0].vx << " (expected: 0.001443829870393)" << std::endl;
      std::cout << "    vy: " << coords[0].vy << " (expected: 0.019869503672238)" << std::endl;
      std::cout << "    vz: " << coords[0].vz << " (expected: 0.0)" << std::endl;

      if (std::abs(coords[0].x - 0.989293478931466) < 1e-12 &&
          std::abs(coords[0].vx - 0.001443829870393) < 1e-12) {
        std::cout << "    => PASS" << std::endl;
      } else {
        std::cout << "    => FAIL" << std::endl;
        all_passed = false;
      }
    }

  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    all_passed = false;
  }

  std::cout << "\n========================================" << std::endl;
  if (all_passed) {
    std::cout << "All tests PASSED!" << std::endl;
    return 0;
  } else {
    std::cout << "Some tests FAILED!" << std::endl;
    return 1;
  }
}
