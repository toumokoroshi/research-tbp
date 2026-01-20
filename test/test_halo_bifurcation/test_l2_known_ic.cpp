/**
 * @file test_l2_known_ic.cpp
 * @brief 既知のL2 Lyapunov軌道初期条件でPoincaré map・積分器の動作確認
 * @details 文献に基づいた初期条件を使用して、軌道が正しくy=0断面に戻るか検証
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "continuation.hpp"
#include "periodic_orbit_stability.hpp"
#include "periodic_orbit_impl.hpp"
#include "rtbp.hpp"
#include "lyapunov_initial.hpp"
#include "utils.hpp"

using namespace continuation;
using namespace periodic_orbit;
using namespace lyapunov;
using namespace my_type;
using namespace param;

// 太陽-地球系の質量パラメータ（astro_param.txtから計算）
// mu = GM_earth / (GM_earth + GM_sun)
double kMuSunEarth = 0.0;  // main()で初期化

// テスト結果カウンタ
int tests_passed = 0;
int tests_failed = 0;

#define TEST_CASE(name) std::cout << "\n=== Test: " << name << " ===" << std::endl;

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

/**
 * @brief CR3BP運動方程式で状態を積分
 * @param state 初期状態
 * @param mu 質量パラメータ
 * @param t_end 積分終了時刻
 * @param dt 時刻刻み
 * @return 最終状態
 */
State<double> IntegrateState(const State<double>& state, double mu, double t_end, double dt) {
  using namespace boost::numeric::odeint;
  using Stepper = runge_kutta_dopri5<State<double>, double, State<double>, double, vector_space_algebra>;
  using ControlledStepper = controlled_runge_kutta<Stepper>;

  auto eom = [mu](const State<double>& s, State<double>& dsdt, double /*t*/) {
    double x1 = s.x + mu;
    double x2 = s.x - (1.0 - mu);
    double r1_sq = x1 * x1 + s.y * s.y + s.z * s.z;
    double r2_sq = x2 * x2 + s.y * s.y + s.z * s.z;
    double r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
    double r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

    dsdt.x = s.vx;
    dsdt.y = s.vy;
    dsdt.z = s.vz;
    dsdt.vx = 2.0 * s.vy + s.x - (1.0 - mu) * x1 * r1_inv3 - mu * x2 * r2_inv3;
    dsdt.vy = -2.0 * s.vx + s.y - (1.0 - mu) * s.y * r1_inv3 - mu * s.y * r2_inv3;
    dsdt.vz = -(1.0 - mu) * s.z * r1_inv3 - mu * s.z * r2_inv3;
  };

  ControlledStepper stepper = make_controlled(1e-14, 1e-14, Stepper());
  State<double> current = state;
  double t = 0.0;

  while (t < t_end) {
    double step = dt;
    auto result = stepper.try_step(eom, current, t, step);
    if (result != controlled_step_result::success) {
      continue;
    }
  }

  return current;
}

/**
 * @brief 半周期積分してy=0交差を検出
 * @return 交差点の状態と時刻
 */
std::pair<State<double>, double> IntegrateToHalfPeriod(
    const State<double>& state, double mu, double max_time, double dt) {
  using namespace boost::numeric::odeint;
  using Stepper = runge_kutta_dopri5<State<double>, double, State<double>, double, vector_space_algebra>;
  using ControlledStepper = controlled_runge_kutta<Stepper>;

  auto eom = [mu](const State<double>& s, State<double>& dsdt, double /*t*/) {
    double x1 = s.x + mu;
    double x2 = s.x - (1.0 - mu);
    double r1_sq = x1 * x1 + s.y * s.y + s.z * s.z;
    double r2_sq = x2 * x2 + s.y * s.y + s.z * s.z;
    double r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
    double r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

    dsdt.x = s.vx;
    dsdt.y = s.vy;
    dsdt.z = s.vz;
    dsdt.vx = 2.0 * s.vy + s.x - (1.0 - mu) * x1 * r1_inv3 - mu * x2 * r2_inv3;
    dsdt.vy = -2.0 * s.vx + s.y - (1.0 - mu) * s.y * r1_inv3 - mu * s.y * r2_inv3;
    dsdt.vz = -(1.0 - mu) * s.z * r1_inv3 - mu * s.z * r2_inv3;
  };

  ControlledStepper stepper = make_controlled(1e-14, 1e-14, Stepper());
  State<double> current = state;
  State<double> prev = state;
  double t = 0.0;
  double prev_t = 0.0;
  double min_time = 0.5;  // 最初のy=0通過を避けるため

  while (t < max_time) {
    prev = current;
    prev_t = t;
    double step = dt;
    auto result = stepper.try_step(eom, current, t, step);
    if (result != controlled_step_result::success) {
      continue;
    }

    // y=0交差を検出（初期状態がy=0から開始なので、逆向きに戻るところを検出）
    // 負のyから正のyへ交差（vy_0 < 0 で開始した場合）
    if (t > min_time && prev.y < 0 && current.y >= 0) {
      // 線形補間で交差点を求める
      double alpha = -prev.y / (current.y - prev.y);
      State<double> crossing;
      crossing.x = prev.x + alpha * (current.x - prev.x);
      crossing.y = 0.0;
      crossing.z = prev.z + alpha * (current.z - prev.z);
      crossing.vx = prev.vx + alpha * (current.vx - prev.vx);
      crossing.vy = prev.vy + alpha * (current.vy - prev.vy);
      crossing.vz = prev.vz + alpha * (current.vz - prev.vz);
      double crossing_t = prev_t + alpha * (t - prev_t);
      return {crossing, crossing_t};
    }
  }

  // 交差が見つからなかった場合
  return {current, -1.0};
}

/**
 * @brief テスト1: L2点の位置計算
 */
void TestL2Position() {
  TEST_CASE("L2 Point Position");

  double x_L2 = ComputeL2Position(kMuSunEarth);
  std::cout << "  L2 position: x = " << std::setprecision(12) << x_L2 << std::endl;

  // L2は地球の外側（x > 1 - mu）
  ASSERT_TRUE(x_L2 > 1.0 - kMuSunEarth, "L2 is beyond Earth");
  // 太陽-地球系では x_L2 ~ 1.01 程度
  ASSERT_TRUE(x_L2 > 1.009 && x_L2 < 1.012, "L2 is approximately at x = 1.01");
}

/**
 * @brief テスト2: Lyapunov係数の計算
 */
void TestLyapunovCoefficients() {
  TEST_CASE("Lyapunov Coefficients for L2");

  auto coeff = ComputeLyapunovCoefficients(LagrangePoint::L2, kMuSunEarth);

  std::cout << "  x_L = " << coeff.x_L << std::endl;
  std::cout << "  c2 = " << coeff.c2 << std::endl;
  std::cout << "  omega_xy = " << coeff.omega_xy << std::endl;
  std::cout << "  kappa = " << coeff.kappa << std::endl;
  std::cout << "  lambda (unstable) = " << coeff.lambda << std::endl;

  double estimated_period = 2.0 * M_PI / coeff.omega_xy;
  std::cout << "  Estimated period = " << estimated_period << std::endl;

  ASSERT_TRUE(coeff.omega_xy > 0, "omega_xy is positive");
  ASSERT_TRUE(coeff.kappa > 0, "kappa is positive");
  ASSERT_TRUE(coeff.lambda > 0, "lambda (unstable eigenvalue) is positive");
  // 太陽-地球L2の周期は約3.1無次元時間
  ASSERT_TRUE(estimated_period > 3.0 && estimated_period < 3.3, "Period is approximately 3.1");
}

/**
 * @brief テスト3: 線形近似初期条件の生成
 */
void TestLinearApproximationIC() {
  TEST_CASE("Linear Approximation Initial Condition");

  double amplitude = 0.001;  // 小振幅
  State<double> ic = GenerateLyapunovInitialGuess(LagrangePoint::L2, kMuSunEarth, amplitude);

  std::cout << "  Initial condition (amplitude = " << amplitude << "):" << std::endl;
  std::cout << "    x  = " << std::setprecision(12) << ic.x << std::endl;
  std::cout << "    y  = " << ic.y << std::endl;
  std::cout << "    z  = " << ic.z << std::endl;
  std::cout << "    vx = " << ic.vx << std::endl;
  std::cout << "    vy = " << ic.vy << std::endl;
  std::cout << "    vz = " << ic.vz << std::endl;

  ASSERT_TRUE(std::abs(ic.y) < 1e-14, "Initial y = 0");
  ASSERT_TRUE(std::abs(ic.z) < 1e-14, "Initial z = 0 (planar)");
  ASSERT_TRUE(std::abs(ic.vx) < 1e-14, "Initial vx = 0");
  ASSERT_TRUE(ic.vy < 0, "Initial vy < 0 (moving in -y direction)");
}

/**
 * @brief テスト4: 半周期積分のテスト（小振幅）
 */
void TestHalfPeriodIntegration() {
  TEST_CASE("Half-Period Integration (Small Amplitude)");

  double amplitude = 0.0001;  // 非常に小振幅
  State<double> ic = GenerateLyapunovInitialGuess(LagrangePoint::L2, kMuSunEarth, amplitude);

  std::cout << "  Initial: x=" << ic.x << ", y=" << ic.y << ", vy=" << ic.vy << std::endl;

  double estimated_period = EstimateLyapunovPeriod(LagrangePoint::L2, kMuSunEarth);
  double max_time = estimated_period * 2.0;

  auto [crossing, half_period] = IntegrateToHalfPeriod(ic, kMuSunEarth, max_time, 0.0001);

  if (half_period > 0) {
    std::cout << "  Half-period crossing at t = " << half_period << std::endl;
    std::cout << "  Crossing state:" << std::endl;
    std::cout << "    x  = " << crossing.x << std::endl;
    std::cout << "    y  = " << crossing.y << std::endl;
    std::cout << "    vx = " << crossing.vx << std::endl;
    std::cout << "    vz = " << crossing.vz << std::endl;

    double vx_error = std::abs(crossing.vx);
    std::cout << "  |vx| at crossing = " << vx_error << std::endl;

    // 周期軌道なら vx ~ 0 になるはず
    ASSERT_TRUE(half_period > 0, "Found y=0 crossing");
    ASSERT_TRUE(half_period > estimated_period * 0.4 && half_period < estimated_period * 0.6,
                "Half-period is approximately T/2");
  } else {
    std::cout << "  WARN: No y=0 crossing found within max_time" << std::endl;
    tests_failed++;
  }
}

/**
 * @brief テスト5: ヤコビ積分の保存
 */
void TestJacobiIntegralConservation() {
  TEST_CASE("Jacobi Integral Conservation");

  double amplitude = 0.0001;
  State<double> ic = GenerateLyapunovInitialGuess(LagrangePoint::L2, kMuSunEarth, amplitude);

  double C_initial = crtbp::calc_jacobi_integral(ic, kMuSunEarth);
  std::cout << "  Initial Jacobi integral: C = " << C_initial << std::endl;

  // 1周期分積分
  double estimated_period = EstimateLyapunovPeriod(LagrangePoint::L2, kMuSunEarth);
  State<double> final_state = IntegrateState(ic, kMuSunEarth, estimated_period, 0.0001);

  double C_final = crtbp::calc_jacobi_integral(final_state, kMuSunEarth);
  std::cout << "  Final Jacobi integral: C = " << C_final << std::endl;
  std::cout << "  Relative error: " << std::abs(C_final - C_initial) / std::abs(C_initial) << std::endl;

  ASSERT_TRUE(std::abs(C_final - C_initial) / std::abs(C_initial) < 1e-10,
              "Jacobi integral conserved to 1e-10");
}

/**
 * @brief テスト6: Poincaré写像のテスト
 */
void TestPoincareMap() {
  TEST_CASE("Poincare Map Application");

  double amplitude = 0.0001;
  State<double> ic = GenerateLyapunovInitialGuess(LagrangePoint::L2, kMuSunEarth, amplitude);

  double estimated_period = EstimateLyapunovPeriod(LagrangePoint::L2, kMuSunEarth);
  double max_time = estimated_period * 3.0;

  std::cout << "  Testing ApplyPoincareMapSafe..." << std::endl;
  std::cout << "  Initial: x=" << ic.x << ", vy=" << ic.vy << std::endl;
  std::cout << "  Estimated period: " << estimated_period << std::endl;

  double period_out = 0.0;
  auto result = ApplyPoincareMapSafe(ic, kMuSunEarth, 1, 0.0, max_time, 0.0001, &period_out);

  if (result) {
    std::cout << "  Poincare map succeeded!" << std::endl;
    std::cout << "  Period: " << period_out << std::endl;
    std::cout << "  Return state:" << std::endl;
    std::cout << "    x  = " << result->x << std::endl;
    std::cout << "    y  = " << result->y << std::endl;
    std::cout << "    vx = " << result->vx << std::endl;
    std::cout << "    vy = " << result->vy << std::endl;

    // 残差を計算
    double residual_x = std::abs(result->x - ic.x);
    double residual_vy = std::abs(result->vy - ic.vy);
    std::cout << "  Residual |x - x0| = " << residual_x << std::endl;
    std::cout << "  Residual |vy - vy0| = " << residual_vy << std::endl;

    ASSERT_TRUE(period_out > estimated_period * 0.8 && period_out < estimated_period * 1.2,
                "Period is close to estimate");
  } else {
    std::cout << "  WARN: Poincare map timed out" << std::endl;
    tests_failed++;
  }
}

/**
 * @brief テスト7: 対称性を利用した微分修正
 */
void TestSymmetricCorrection() {
  TEST_CASE("Symmetric Differential Correction");

  double amplitude = 0.0001;  // 小振幅で開始
  State<double> ic = GenerateLyapunovInitialGuess(LagrangePoint::L2, kMuSunEarth, amplitude);

  std::cout << "  Linear approximation IC:" << std::endl;
  std::cout << "    x0  = " << ic.x << std::endl;
  std::cout << "    vy0 = " << ic.vy << std::endl;

  double estimated_period = EstimateLyapunovPeriod(LagrangePoint::L2, kMuSunEarth);
  double max_time = estimated_period * 2.0;

  try {
    NewtonConvergenceInfo<double> conv_info;
    auto orbit = RefinePeriodicOrbitSymmetric(
        ic, kMuSunEarth,
        20,      // max_iterations
        1e-10,   // tolerance
        max_time,
        0.0001,  // dt
        &conv_info
    );

    std::cout << "  Refined orbit:" << std::endl;
    std::cout << "    x0  = " << orbit.initial_state.x << std::endl;
    std::cout << "    vy0 = " << orbit.initial_state.vy << std::endl;
    std::cout << "    Period = " << orbit.period << std::endl;
    std::cout << "    Jacobi = " << orbit.jacobi_constant << std::endl;
    std::cout << "    Iterations = " << conv_info.iterations << std::endl;
    std::cout << "    Final |vx| = " << conv_info.final_residual << std::endl;

    ASSERT_TRUE(conv_info.converged, "Newton converged");
    ASSERT_TRUE(conv_info.iterations < 15, "Converged in < 15 iterations");
    ASSERT_TRUE(conv_info.final_residual < 1e-9, "Final |vx| < 1e-9");
    ASSERT_TRUE(orbit.period > estimated_period * 0.9 && orbit.period < estimated_period * 1.1,
                "Period is close to linear estimate");

    // 全周期での戻りを検証
    State<double> final_state = IntegrateState(orbit.initial_state, kMuSunEarth, orbit.period, 0.0001);
    double x_error = std::abs(final_state.x - orbit.initial_state.x);
    double vy_error = std::abs(final_state.vy - orbit.initial_state.vy);
    std::cout << "  Full period verification:" << std::endl;
    std::cout << "    |x - x0| = " << x_error << std::endl;
    std::cout << "    |vy - vy0| = " << vy_error << std::endl;

    ASSERT_TRUE(x_error < 1e-6, "Full period x error < 1e-6");
    ASSERT_TRUE(vy_error < 1e-6, "Full period vy error < 1e-6");

  } catch (const std::exception& e) {
    std::cout << "  FAILED: " << e.what() << std::endl;
    tests_failed++;
  }
}

int main() {
  std::cout << std::setprecision(12);
  std::cout << "\n========================================" << std::endl;
  std::cout << "  L2 Lyapunov Known IC Test Suite" << std::endl;
  std::cout << "========================================" << std::endl;

  // 設定ファイルからmu値をロード
  std::string astro_param_file = "../../test/configs/astro_param/astro_param.txt";
  try {
    AstroConstants<double> astro_params = utils::loadConstants<double>(astro_param_file);
    kMuSunEarth = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);
  } catch (const std::exception& e) {
    std::cerr << "Failed to load astro_param.txt: " << e.what() << std::endl;
    // フォールバック値を使用
    kMuSunEarth = 3.003480593e-06;
    std::cout << "  Using fallback mu value" << std::endl;
  }

  std::cout << "  mu = " << kMuSunEarth << " (Sun-Earth system)" << std::endl;

  TestL2Position();
  TestLyapunovCoefficients();
  TestLinearApproximationIC();
  TestHalfPeriodIntegration();
  TestJacobiIntegralConservation();
  TestPoincareMap();
  TestSymmetricCorrection();

  std::cout << "\n========================================" << std::endl;
  std::cout << "  Test Summary" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "  Passed: " << tests_passed << std::endl;
  std::cout << "  Failed: " << tests_failed << std::endl;
  std::cout << "  Total:  " << (tests_passed + tests_failed) << std::endl;
  std::cout << "========================================" << std::endl;

  return (tests_failed == 0) ? 0 : 1;
}
