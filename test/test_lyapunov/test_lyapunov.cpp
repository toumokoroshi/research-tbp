/**
 * @file test_lyapunov.cpp
 * @brief Lyapunov軌道初期条件生成のテスト
 * @details 物理的・プログラム的健全性を検証するテストケース
 * @date 2026-01-01
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

#include "lyapunov_initial.hpp"
#include "periodic_orbit.hpp"
#include "rtbp.hpp"

using namespace lyapunov;
using namespace my_type;
using namespace crtbp;

// テスト結果カウンタ
int tests_passed = 0;
int tests_failed = 0;

// テスト用マクロ
#define TEST_CASE(name) std::cout << "\n=== Test: " << name << " ===" << std::endl;

#define ASSERT_NEAR(actual, expected, tolerance, msg)                                    \
  do {                                                                                   \
    double a = (actual);                                                                 \
    double e = (expected);                                                               \
    double t = (tolerance);                                                              \
    if (std::abs(a - e) < t) {                                                           \
      std::cout << "  PASS: " << msg << " (actual=" << a << ", expected=" << e << ")\n"; \
      tests_passed++;                                                                    \
    } else {                                                                             \
      std::cout << "  FAIL: " << msg << " (actual=" << a << ", expected=" << e           \
                << ", diff=" << std::abs(a - e) << ")\n";                                \
      tests_failed++;                                                                    \
    }                                                                                    \
  } while (0)

#define ASSERT_TRUE(condition, msg)           \
  do {                                        \
    if (condition) {                          \
      std::cout << "  PASS: " << msg << "\n"; \
      tests_passed++;                         \
    } else {                                  \
      std::cout << "  FAIL: " << msg << "\n"; \
      tests_failed++;                         \
    }                                         \
  } while (0)

// ---------------------------------------------------------------------------
// Test 1: L点位置の計算精度
// ---------------------------------------------------------------------------
void TestLagrangePointPositions() {
  TEST_CASE("Lagrange Point Positions (Sun-Earth System)");

  // Sun-Earth系の質量比 (μ = M_Earth / (M_Sun + M_Earth))
  const double mu = 3.0404233870218e-06;

  // 文献値 (Koon et al., "Dynamical Systems, the Three-Body Problem and Space Mission Design")
  // L1: x ≈ 0.9900 (地球-太陽間)
  // L2: x ≈ 1.0100 (地球の外側)
  // L3: x ≈ -1.0000 (太陽の反対側)

  double x_L1 = ComputeL1Position(mu);
  double x_L2 = ComputeL2Position(mu);
  double x_L3 = ComputeL3Position(mu);

  std::cout << std::setprecision(12);
  std::cout << "  L1 position: x = " << x_L1 << "\n";
  std::cout << "  L2 position: x = " << x_L2 << "\n";
  std::cout << "  L3 position: x = " << x_L3 << "\n";

  // L1: 地球と太陽の間、約0.990 (1-mu から約0.01離れた位置)
  ASSERT_TRUE(x_L1 > 0.98 && x_L1 < 1.0, "L1 is between Sun and Earth");
  ASSERT_TRUE(x_L1 < 1.0 - mu, "L1 is closer to Sun than Earth");

  // L2: 地球の外側、約1.010
  ASSERT_TRUE(x_L2 > 1.0 - mu, "L2 is beyond Earth");
  ASSERT_TRUE(x_L2 < 1.02, "L2 is within expected range");

  // L3: 太陽の反対側、約-1.0
  ASSERT_TRUE(x_L3 < -0.99, "L3 is on opposite side of Sun");
  ASSERT_TRUE(x_L3 > -1.01, "L3 is within expected range");

  // Hill球半径との整合性チェック
  // r_Hill ≈ (μ/3)^(1/3) ≈ 0.01 for Sun-Earth
  double r_hill = std::pow(mu / 3.0, 1.0 / 3.0);
  std::cout << "  Hill radius: " << r_hill << "\n";

  // L1とL2は地球からHill球半径程度離れているはず
  double dist_L1_to_Earth = std::abs((1.0 - mu) - x_L1);
  double dist_L2_to_Earth = std::abs(x_L2 - (1.0 - mu));
  ASSERT_NEAR(dist_L1_to_Earth, r_hill, 0.005, "L1 distance from Earth ~ Hill radius");
  ASSERT_NEAR(dist_L2_to_Earth, r_hill, 0.005, "L2 distance from Earth ~ Hill radius");
}

// ---------------------------------------------------------------------------
// Test 2: 平衡点条件の検証
// ---------------------------------------------------------------------------
void TestEquilibriumCondition() {
  TEST_CASE("Equilibrium Point Condition (Force Balance)");

  const double mu = 3.0404233870218e-06;

  // 各L点で力の釣り合いを確認
  // CR3BPの平衡点条件: x - (1-μ)(x+μ)/r1³ - μ(x-(1-μ))/r2³ = 0
  // y方向と z方向は y=0, z=0 で自動的に満たされる

  for (auto point : {LagrangePoint::L1, LagrangePoint::L2, LagrangePoint::L3}) {
    double x = ComputeLagrangePointPosition(point, mu);
    double y = 0.0;

    double x1 = x + mu;
    double x2 = x - (1.0 - mu);
    double r1 = std::abs(x1);
    double r2 = std::abs(x2);

    // 有効ポテンシャルの勾配（平衡点では0になるはず）
    double grad_x = x - (1.0 - mu) * x1 / (r1 * r1 * r1) - mu * x2 / (r2 * r2 * r2);

    std::string point_name = (point == LagrangePoint::L1)   ? "L1"
                             : (point == LagrangePoint::L2) ? "L2"
                                                            : "L3";
    ASSERT_NEAR(grad_x, 0.0, 1e-10, point_name + " x-gradient = 0");
  }
}

// ---------------------------------------------------------------------------
// Test 3: Lyapunov係数の物理的妥当性
// ---------------------------------------------------------------------------
void TestLyapunovCoefficients() {
  TEST_CASE("Lyapunov Coefficients Physical Validity");

  const double mu = 3.0404233870218e-06;

  for (auto point : {LagrangePoint::L1, LagrangePoint::L2}) {
    auto coeff = ComputeLyapunovCoefficients<double>(point, mu);

    std::string point_name = (point == LagrangePoint::L1) ? "L1" : "L2";
    std::cout << "  " << point_name << ":\n";
    std::cout << "    c2 = " << coeff.c2 << "\n";
    std::cout << "    omega_xy = " << coeff.omega_xy << "\n";
    std::cout << "    omega_z = " << coeff.omega_z << "\n";
    std::cout << "    lambda = " << coeff.lambda << "\n";
    std::cout << "    kappa = " << coeff.kappa << "\n";

    // c2 > 0 (attraction to L-point in some directions)
    ASSERT_TRUE(coeff.c2 > 0, point_name + " c2 > 0");

    // omega_xy: 面内振動周波数（実数、正）
    ASSERT_TRUE(coeff.omega_xy > 0, point_name + " omega_xy > 0");
    ASSERT_TRUE(coeff.omega_xy < 5.0, point_name + " omega_xy reasonable range");

    // lambda: サドル型固有値（正の実数）
    ASSERT_TRUE(coeff.lambda > 0, point_name + " lambda > 0 (saddle type)");

    // kappa: 振幅比（通常2〜4程度）
    ASSERT_TRUE(coeff.kappa > 1.0, point_name + " kappa > 1");
    ASSERT_TRUE(coeff.kappa < 10.0, point_name + " kappa < 10");

    // 周期の妥当性（L1/L2周りの周期は約3〜4無次元時間）
    double period = 2.0 * std::numbers::pi / coeff.omega_xy;
    std::cout << "    Estimated period = " << period << "\n";
    ASSERT_TRUE(period > 2.0, point_name + " period > 2 (dimensionless)");
    ASSERT_TRUE(period < 5.0, point_name + " period < 5 (dimensionless)");
  }
}

// ---------------------------------------------------------------------------
// Test 4: 初期条件の対称性
// ---------------------------------------------------------------------------
void TestInitialConditionSymmetry() {
  TEST_CASE("Lyapunov Initial Condition Symmetry");

  const double mu = 3.0404233870218e-06;
  const double amplitude = 0.001;

  auto state = GenerateLyapunovInitialGuess<double>(LagrangePoint::L1, mu, amplitude);

  std::cout << "  Initial state:\n";
  std::cout << "    x  = " << state.x << "\n";
  std::cout << "    y  = " << state.y << "\n";
  std::cout << "    z  = " << state.z << "\n";
  std::cout << "    vx = " << state.vx << "\n";
  std::cout << "    vy = " << state.vy << "\n";
  std::cout << "    vz = " << state.vz << "\n";

  // Lyapunov軌道の対称性条件（平面軌道、x軸対称）
  // t=0 で: y=0, vx=0, z=0, vz=0
  ASSERT_NEAR(state.y, 0.0, 1e-15, "y = 0 at t=0");
  ASSERT_NEAR(state.vx, 0.0, 1e-15, "vx = 0 at t=0");
  ASSERT_NEAR(state.z, 0.0, 1e-15, "z = 0 (planar orbit)");
  ASSERT_NEAR(state.vz, 0.0, 1e-15, "vz = 0 (planar orbit)");

  // x座標はL1点 + amplitude
  auto coeff = ComputeLyapunovCoefficients<double>(LagrangePoint::L1, mu);
  ASSERT_NEAR(state.x, coeff.x_L + amplitude, 1e-12, "x = x_L + amplitude");

  // vy の符号と大きさ（κ * ω * amplitude）
  double expected_vy = -coeff.kappa * coeff.omega_xy * amplitude;
  ASSERT_NEAR(state.vy, expected_vy, 1e-12, "vy = -kappa * omega * amplitude");
}

// ---------------------------------------------------------------------------
// Test 5: ヤコビ積分の保存
// ---------------------------------------------------------------------------
void TestJacobiIntegralConsistency() {
  TEST_CASE("Jacobi Integral Consistency");

  const double mu = 3.0404233870218e-06;
  const double amplitude = 0.001;

  auto state = GenerateLyapunovInitialGuess<double>(LagrangePoint::L1, mu, amplitude);

  // ヤコビ積分を計算
  double C = calc_jacobi_integral(state, mu);
  std::cout << "  Jacobi integral C = " << C << "\n";

  // L1点でのヤコビ積分（速度=0）と比較
  auto coeff = ComputeLyapunovCoefficients<double>(LagrangePoint::L1, mu);
  State<double> L1_state = {coeff.x_L, 0, 0, 0, 0, 0};
  double C_L1 = calc_jacobi_integral(L1_state, mu);
  std::cout << "  Jacobi integral at L1 = " << C_L1 << "\n";

  // Lyapunov軌道はL1点より低いエネルギー（高いヤコビ積分）を持つ...
  // 実際には振幅があると運動エネルギーが増えるのでCは減少する
  // C = 2U* - v² なので、v²が増えるとCは減る
  ASSERT_TRUE(C < C_L1, "C(Lyapunov orbit) < C(L1 point)");

  // Sun-Earth系のL1ではC ≈ 3.00
  ASSERT_TRUE(C > 2.9, "C > 2.9 (reasonable for Sun-Earth L1)");
  ASSERT_TRUE(C < 3.1, "C < 3.1 (reasonable for Sun-Earth L1)");
}

// ---------------------------------------------------------------------------
// Test 6: 異なる振幅での線形スケーリング
// ---------------------------------------------------------------------------
void TestAmplitudeScaling() {
  TEST_CASE("Amplitude Scaling");

  const double mu = 3.0404233870218e-06;

  auto state1 = GenerateLyapunovInitialGuess<double>(LagrangePoint::L1, mu, 0.001);
  auto state2 = GenerateLyapunovInitialGuess<double>(LagrangePoint::L1, mu, 0.002);

  // x座標の変化は振幅に比例
  auto coeff = ComputeLyapunovCoefficients<double>(LagrangePoint::L1, mu);
  double dx1 = state1.x - coeff.x_L;
  double dx2 = state2.x - coeff.x_L;
  ASSERT_NEAR(dx2 / dx1, 2.0, 1e-10, "x displacement scales linearly with amplitude");

  // vy の変化も振幅に比例
  ASSERT_NEAR(state2.vy / state1.vy, 2.0, 1e-10, "vy scales linearly with amplitude");
}

// ---------------------------------------------------------------------------
// Test 7: オービットファミリータイプの一貫性
// ---------------------------------------------------------------------------
void TestOrbitFamilyTypeConsistency() {
  TEST_CASE("Orbit Family Type Consistency");

  const double mu = 3.0404233870218e-06;

  {
    LyapunovOrbitGenerator<double> gen(LagrangePoint::L1, mu);
    ASSERT_TRUE(gen.GetFamilyType() == OrbitFamilyType::LYAPUNOV_L1, "L1 generator -> LYAPUNOV_L1");
  }

  {
    LyapunovOrbitGenerator<double> gen(LagrangePoint::L2, mu);
    ASSERT_TRUE(gen.GetFamilyType() == OrbitFamilyType::LYAPUNOV_L2, "L2 generator -> LYAPUNOV_L2");
  }

  {
    LyapunovOrbitGenerator<double> gen(LagrangePoint::L3, mu);
    ASSERT_TRUE(gen.GetFamilyType() == OrbitFamilyType::LYAPUNOV_L3, "L3 generator -> LYAPUNOV_L3");
  }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main() {
  std::cout << "========================================\n";
  std::cout << "  Lyapunov Initial Condition Tests\n";
  std::cout << "  Physical and Programmatic Validation\n";
  std::cout << "========================================\n";

  TestLagrangePointPositions();
  TestEquilibriumCondition();
  TestLyapunovCoefficients();
  TestInitialConditionSymmetry();
  TestJacobiIntegralConsistency();
  TestAmplitudeScaling();
  TestOrbitFamilyTypeConsistency();

  std::cout << "\n========================================\n";
  std::cout << "  Test Summary\n";
  std::cout << "========================================\n";
  std::cout << "  Passed: " << tests_passed << "\n";
  std::cout << "  Failed: " << tests_failed << "\n";
  std::cout << "  Total:  " << (tests_passed + tests_failed) << "\n";
  std::cout << "========================================\n";

  return (tests_failed == 0) ? 0 : 1;
}
