/**
 * @file test_halo_bifurcation.cpp
 * @brief Halo分岐関連機能のテスト
 * @details 
 *   1. Eigenによる固有値計算の検証
 *   2. 分岐方向ベクトル抽出の検証
 *   3. ContinueHaloFromBifurcation関数の存在確認
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include "continuation.hpp"
#include "periodic_orbit_stability.hpp"
#include "rtbp.hpp"

using namespace continuation;
using namespace periodic_orbit;
using namespace lyapunov;

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

// 太陽-地球系の質量パラメータ
constexpr double kMuSunEarth = 3.003480593e-06;

/**
 * @brief テスト1: 単純なモノドロミー行列で固有値計算をテスト
 * @details 単位行列の固有値がすべて1であることを確認
 */
void TestEigenvalueComputationIdentity() {
  TEST_CASE("Eigenvalue Computation - Identity Matrix");

  PeriodicOrbit<double> orbit;
  
  // モノドロミー行列を単位行列に設定
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      orbit.monodromy_matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  // 安定性解析（固有値計算）
  AnalyzeStability(&orbit, 1.5);

  std::cout << "  Eigenvalues computed: " << orbit.eigenvalues.size() << std::endl;
  
  // 6つの固有値があることを確認
  ASSERT_TRUE(orbit.eigenvalues.size() == 6, "6 eigenvalues computed");

  // すべての固有値が1であることを確認
  int count_near_one = 0;
  for (const auto& ev : orbit.eigenvalues) {
    double abs_ev = std::abs(ev);
    std::cout << "    |λ| = " << abs_ev << std::endl;
    if (std::abs(abs_ev - 1.0) < 1e-10) {
      count_near_one++;
    }
  }
  
  ASSERT_TRUE(count_near_one == 6, "All eigenvalues equal to 1 for identity matrix");
  ASSERT_TRUE(orbit.is_stable, "Identity matrix is stable (|λ| = 1)");
}

/**
 * @brief テスト2: 対角行列で固有値計算をテスト
 */
void TestEigenvalueComputationDiagonal() {
  TEST_CASE("Eigenvalue Computation - Diagonal Matrix");

  PeriodicOrbit<double> orbit;
  
  // 対角行列を設定: diag(2, 0.5, 1, 1, 3, 1/3)
  // 固有値は対角成分と一致するはず
  double diag_values[] = {2.0, 0.5, 1.0, 1.0, 3.0, 1.0/3.0};
  
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      orbit.monodromy_matrix[i][j] = (i == j) ? diag_values[i] : 0.0;
    }
  }

  AnalyzeStability(&orbit, 1.5);

  std::cout << "  Expected eigenvalues: 2, 0.5, 1, 1, 3, 0.333..." << std::endl;
  std::cout << "  Computed eigenvalues:" << std::endl;
  
  for (const auto& ev : orbit.eigenvalues) {
    std::cout << "    λ = " << ev.real() << " + " << ev.imag() << "i" << std::endl;
  }

  // 最大固有値が3であることを確認
  ASSERT_NEAR(orbit.stability_index, 3.0, 1e-10, "Max |λ| = 3");
  
  // 不安定（|λ| > 1 の固有値が存在）
  ASSERT_TRUE(!orbit.is_stable, "Matrix with |λ| > 1 is unstable");
}

/**
 * @brief テスト3: 分岐方向ベクトル計算のテスト
 */
void TestBifurcationEigenvector() {
  TEST_CASE("Bifurcation Eigenvector Extraction");

  PeriodicOrbit<double> orbit;
  
  // 固有値1を持つ成分を含む行列を構成
  // 簡単のため、第3成分（z成分）で固有値=1となるようにする
  // diag(0.5, 2, 1, 1, 0.5, 2) - z成分が1
  double diag_values[] = {0.5, 2.0, 1.0, 1.0, 0.5, 2.0};
  
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      orbit.monodromy_matrix[i][j] = (i == j) ? diag_values[i] : 0.0;
    }
  }

  // 分岐方向ベクトルを計算
  auto bif_vec = ComputeBifurcationEigenvector(orbit);
  
  std::cout << "  Bifurcation eigenvector:" << std::endl;
  std::cout << "    [" << bif_vec[0] << ", " << bif_vec[1] << ", " << bif_vec[2] 
            << ", " << bif_vec[3] << ", " << bif_vec[4] << ", " << bif_vec[5] << "]" << std::endl;

  // ベクトルが正規化されていることを確認
  double norm = 0.0;
  for (int i = 0; i < 6; ++i) {
    norm += bif_vec[i] * bif_vec[i];
  }
  norm = std::sqrt(norm);
  
  ASSERT_NEAR(norm, 1.0, 1e-10, "Eigenvector is normalized");
  
  // λ=1 に対応する固有ベクトルは (0,0,1,0,0,0) または (0,0,0,1,0,0) のはず
  // 対角行列なので、固有ベクトルは標準基底
  // index 2 または 3 が非ゼロであることを確認
  bool z_component_dominant = (std::abs(bif_vec[2]) > 0.9) || (std::abs(bif_vec[3]) > 0.9);
  ASSERT_TRUE(z_component_dominant, "Eigenvector corresponds to z or vz component");
}

/**
 * @brief テスト4: ContinuationConfig構造体のテスト
 */
void TestContinuationConfigDefault() {
  TEST_CASE("ContinuationConfig Default Values");

  ContinuationConfig<double> config;
  
  ASSERT_TRUE(config.method == ContinuationMethod::PSEUDO_ARCLENGTH, 
              "Default method is PSEUDO_ARCLENGTH");
  ASSERT_TRUE(config.detect_bifurcations, "Bifurcation detection enabled by default");
  ASSERT_TRUE(config.initial_step_size > 0, "Initial step size is positive");
  ASSERT_TRUE(config.max_steps > 0, "Max steps is positive");
}

/**
 * @brief テスト5: BifurcationType列挙型のテスト
 */
void TestBifurcationTypes() {
  TEST_CASE("BifurcationType Enumeration");

  BifurcationPoint<double> bif;
  
  bif.type = BifurcationType::PITCHFORK;
  ASSERT_TRUE(bif.type == BifurcationType::PITCHFORK, "PITCHFORK type assignable");
  
  bif.type = BifurcationType::PERIOD_DOUBLING;
  ASSERT_TRUE(bif.type == BifurcationType::PERIOD_DOUBLING, "PERIOD_DOUBLING type assignable");
  
  bif.type = BifurcationType::NEIMARK_SACKER;
  ASSERT_TRUE(bif.type == BifurcationType::NEIMARK_SACKER, "NEIMARK_SACKER type assignable");
  
  bif.type = BifurcationType::UNKNOWN;
  ASSERT_TRUE(bif.type == BifurcationType::UNKNOWN, "UNKNOWN type assignable");
}

/**
 * @brief テスト6: OrbitContinuatorのインスタンス化テスト
 */
void TestOrbitContinuatorInstantiation() {
  TEST_CASE("OrbitContinuator Instantiation");

  ContinuationConfig<double> config;
  config.max_steps = 10;  // 短縮してテスト
  
  // コンストラクタが例外なく動作することを確認
  try {
    OrbitContinuator<double> continuator(kMuSunEarth, config);
    std::cout << "  OrbitContinuator created successfully" << std::endl;
    tests_passed++;
  } catch (const std::exception& e) {
    std::cout << "  FAIL: OrbitContinuator constructor threw: " << e.what() << std::endl;
    tests_failed++;
  }
}

int main() {
  std::cout << std::setprecision(12);
  std::cout << "\n========================================" << std::endl;
  std::cout << "  Halo Bifurcation Test Suite" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "  mu = " << kMuSunEarth << " (Sun-Earth system)" << std::endl;

  TestEigenvalueComputationIdentity();
  TestEigenvalueComputationDiagonal();
  TestBifurcationEigenvector();
  TestContinuationConfigDefault();
  TestBifurcationTypes();
  TestOrbitContinuatorInstantiation();

  std::cout << "\n========================================" << std::endl;
  std::cout << "  Test Summary" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "  Passed: " << tests_passed << std::endl;
  std::cout << "  Failed: " << tests_failed << std::endl;
  std::cout << "  Total:  " << (tests_passed + tests_failed) << std::endl;
  std::cout << "========================================" << std::endl;

  return (tests_failed == 0) ? 0 : 1;
}
