/**
 * @file periodic_orbit.hpp
 * @brief 周期軌道の検出・安定性解析・不変多様体計算のためのライブラリ
 * @details CR3BP（円制限三体問題）における周期軌道の解析ツール。
 *
 * @par 理論的背景
 * 周期軌道はポアンカレ写像 \f$ P: \Sigma \to \Sigma \f$ の不動点として
 * 定義される。すなわち、断面 \f$ \Sigma \f$ を出発し、流れに沿って
 * 再び \f$ \Sigma \f$ に戻ったときに同一状態となる点が周期軌道の初期条件。
 * Newton-Raphson法で \f$ P(\mathbf{x}) - \mathbf{x} = 0 \f$ を解くことで精密化する。
 *
 * @par 参考文献
 * - [Koon2011] Koon, W.S. et al. "Dynamical Systems, the Three-Body Problem
 *   and Space Mission Design" (Marsden Books, 2011) - Ch.6: 周期軌道と不変多様体
 * - [Parker2006] Parker, J.S., Lo, M.W. "Practical Methods for Libration Point
 *   Orbit Analysis" - Ch.4: 微分修正法
 * - [Meyer2009] Meyer, K.R. et al. "Introduction to Hamiltonian Dynamical
 *   Systems and the N-Body Problem" (Springer, 2009)
 *
 * @date 2025-12-24
 */

#ifndef PERIODIC_ORBIT_HPP
#define PERIODIC_ORBIT_HPP

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "rtbp.hpp"
#include "vector3d.hpp"

namespace periodic_orbit {

using namespace my_type;

// ---------------------------------------------------------------------------
// 構造体定義
// ---------------------------------------------------------------------------

/**
 * @brief 周期軌道の情報を保持する構造体
 *
 * @par モノドロミー行列
 * モノドロミー行列 \f$ M = \Phi(T, 0) \f$ は、周期軌道の初期点からの
 * 微小摂動が1周期後にどう変化するかを記述する状態遷移行列である:
 * \f[ \delta\mathbf{x}(T) = M \cdot \delta\mathbf{x}(0) \f]
 *
 * @par ハミルトン系の制約
 * CR3BPはハミルトン系であるため、\f$ M \f$ はシンプレクティック行列となり、
 * \f$ \det(M) = 1 \f$ かつ固有値は必ず逆数ペア \f$ (\lambda, 1/\lambda) \f$ で現れる。
 * 周期軌道自体の方向には常に \f$ \lambda = 1 \f$ の固有値のペアが存在する。
 *
 * @see [Meyer2009] シンプレクティック構造と固有値の制約
 * @see [Koon2011] Ch.6: モノドロミー行列と安定性
 */
template <typename ScalarType>
struct PeriodicOrbit {
  State<ScalarType> initial_state;   ///< 初期状態 \f$ (x, y, z, \dot{x}, \dot{y}, \dot{z}) \f$ [無次元、回転座標系]
  ScalarType period = 0.0;           ///< 周期 \f$ T \f$ [無次元時間]
  ScalarType jacobi_constant = 0.0;  ///< ヤコビ積分の値 \f$ C_J \f$

  // 安定性解析結果
  bool stability_computed = false;                    ///< 安定性解析が実行されたか
  std::vector<std::complex<ScalarType>> eigenvalues;  ///< モノドロミー行列の固有値 \f$ \lambda_i \f$
  bool is_stable = false;                             ///< 全 \f$ |\lambda_i| \leq 1 \f$ なら安定 (SPO)
  ScalarType stability_index = 0.0;                   ///< 安定性指数 \f$ \max_i |\lambda_i| \f$

  // モノドロミー行列 (6x6)
  std::array<std::array<ScalarType, 6>, 6> monodromy_matrix;

  PeriodicOrbit() {
    // モノドロミー行列をゼロ初期化
    for (auto& row : monodromy_matrix) {
      row.fill(0.0);
    }
  }
};

/**
 * @brief Newton-Raphson法の収束情報を保持する構造体
 */
template <typename ScalarType>
struct NewtonConvergenceInfo {
  bool converged = false;           ///< 収束したか
  int iterations = 0;               ///< 反復回数
  ScalarType final_residual = 0.0;  ///< 最終残差
  std::string message;              ///< 収束状態のメッセージ
};

/**
 * @brief 不変多様体の設定を保持する構造体
 *
 * @par 理論的背景
 * 不安定周期軌道 (UPO) の安定/不安定多様体は、モノドロミー行列の
 * \f$ |\lambda| < 1 \f$ / \f$ |\lambda| > 1 \f$ の固有ベクトル方向に沿って存在する。
 *
 * @par 安定多様体定理
 * 周期軌道 \f$ \gamma \f$ の各点 \f$ x_0 \f$ において、安定多様体 \f$ W^s(\gamma) \f$ は
 * 安定固有空間 \f$ E^s \f$ に接する滑らかな曲面である:
 * \f[ T_{x_0} W^s(\gamma) = E^s(x_0) \f]
 * したがって、十分小さい \f$ \varepsilon \f$ に対して
 * \f$ x_0 + \varepsilon \mathbf{v}_s \f$ は近似的に \f$ W^s \f$ 上にあり、
 * 線形近似の誤差は \f$ O(\varepsilon^2) \f$ である。
 *
 * @see [Koon2011] Ch.6: 不変多様体の計算手法
 * @see [Parker2006] 数値的不変多様体計算の実装
 */
template <typename ScalarType>
struct ManifoldConfig {
  ScalarType epsilon = 1e-6;          ///< 固有ベクトル方向への初期変位量 \f$ \varepsilon \f$
  ScalarType forward_time = 100.0;    ///< 正方向積分時間 [無次元時間]
  ScalarType backward_time = -100.0;  ///< 負方向積分時間 [無次元時間]
  int num_initial_points = 20;        ///< 周期軌道上の初期点数
  bool compute_stable = true;         ///< 安定多様体 \f$ W^s \f$ を計算するか
  bool compute_unstable = true;       ///< 不安定多様体 \f$ W^u \f$ を計算するか
  std::string integrator_type = "rk4";     ///< "rk4" or "rk45"
  ScalarType abs_tolerance = 1e-12;        ///< RK45用 abs tol
  ScalarType rel_tolerance = 1e-12;        ///< RK45用 rel tol
};

/**
 * @brief 不変多様体上の軌道データを保持する構造体
 */
template <typename ScalarType>
struct ManifoldTrajectory {
  enum class Type { STABLE, UNSTABLE };

  Type type;                                  ///< 多様体の種類
  State<ScalarType> initial_displacement;     ///< 初期変位
  std::vector<State<ScalarType>> trajectory;  ///< 軌道データ
  std::vector<ScalarType> times;              ///< 時刻データ
};

// ---------------------------------------------------------------------------
// 関数宣言
// ---------------------------------------------------------------------------

/**
 * @brief ポアンカレ写像の不動点をNewton-Raphson法で精密化
 * @tparam ScalarType スカラー型 (double等)
 * @param initial_guess 初期推定値 (ポアンカレ断面上の状態)
 * @param mu CR3BPの質量比パラメータ
 * @param section_index ポアンカレ断面の変数インデックス (0:x, 1:y, 2:z)
 * @param section_value 断面の値
 * @param max_iterations 最大反復回数
 * @param tolerance 収束判定の許容誤差
 * @param max_time 積分最大時間
 * @param dt 積分時間刻み
 * @param convergence_info 収束情報 (出力)
 * @return 精密化された周期軌道
 */
template <typename ScalarType>
PeriodicOrbit<ScalarType> RefinePeriodicOrbit(
    const State<ScalarType>& initial_guess, ScalarType mu, int section_index,
    ScalarType section_value, int max_iterations, ScalarType tolerance, ScalarType max_time,
    ScalarType dt, NewtonConvergenceInfo<ScalarType>* convergence_info = nullptr);

/**
 * @brief モノドロミー行列を数値計算 (STM: State Transition Matrix)
 * @tparam ScalarType スカラー型 (double等)
 * @param orbit 周期軌道
 * @param mu CR3BPの質量比パラメータ
 * @param dt 積分の時間刻み
 * @return モノドロミー行列 (6x6)
 */
template <typename ScalarType>
std::array<std::array<ScalarType, 6>, 6> ComputeMonodromyMatrix(
    const PeriodicOrbit<ScalarType>& orbit, ScalarType mu, ScalarType dt = 0.0001);

/**
 * @brief 固有値解析による安定性判定
 * @tparam ScalarType スカラー型 (double等)
 * @param orbit 周期軌道（モノドロミー行列が計算済みであること）
 * @param eigenvalue_threshold 不安定判定のしきい値 (|λ| > threshold で不安定)
 */
template <typename ScalarType>
void AnalyzeStability(PeriodicOrbit<ScalarType>* orbit, ScalarType eigenvalue_threshold = 1.0);

/**
 * @brief 不変多様体を計算
 * @tparam ScalarType スカラー型 (double等)
 * @param orbit 周期軌道（安定性解析が完了していること）
 * @param mu CR3BPの質量比パラメータ
 * @param config 多様体計算の設定
 * @param dt 積分の時間刻み
 * @return 多様体上の軌道データのリスト
 */
template <typename ScalarType>
std::vector<ManifoldTrajectory<ScalarType>> ComputeInvariantManifolds(
    const PeriodicOrbit<ScalarType>& orbit, ScalarType mu, const ManifoldConfig<ScalarType>& config,
    ScalarType dt = 0.0001);

// ---------------------------------------------------------------------------
// ヘルパー関数宣言
// ---------------------------------------------------------------------------

/**
 * @brief ポアンカレ写像を1回適用（初期状態から次の断面交差まで積分）
 * @tparam ScalarType スカラー型
 * @param state 初期状態
 * @param mu CR3BPの質量比
 * @param section_index 断面変数のインデックス
 * @param section_value 断面の値
 * @param max_time 最大積分時間
 * @param dt 時間刻み
 * @param period 周期（出力）
 * @return 次の断面交差点の状態
 */
template <typename ScalarType>
State<ScalarType> ApplyPoincareMap(const State<ScalarType>& state, ScalarType mu, int section_index,
                                   ScalarType section_value, ScalarType max_time, ScalarType dt,
                                   ScalarType* period);

/**
 * @brief 行列の固有値を計算（QR法またはPower法）
 * @tparam ScalarType スカラー型
 * @param matrix 正方行列 (6x6)
 * @return 固有値のベクトル（複素数）
 */
template <typename ScalarType>
std::vector<std::complex<ScalarType>> ComputeEigenvalues(
    const std::array<std::array<ScalarType, 6>, 6>& matrix);

/**
 * @brief 固有ベクトルを計算
 * @tparam ScalarType スカラー型
 * @param matrix 正方行列 (6x6)
 * @param eigenvalue 対応する固有値
 * @return 固有ベクトル (6次元)
 */
template <typename ScalarType>
std::array<ScalarType, 6> ComputeEigenvector(const std::array<std::array<ScalarType, 6>, 6>& matrix,
                                             const std::complex<ScalarType>& eigenvalue);

// ---------------------------------------------------------------------------
// 実装部（ヘッダーオンリーライブラリ）
// ---------------------------------------------------------------------------

/**
 * @brief ポアンカレ写像を1回適用
 */
template <typename ScalarType>
State<ScalarType> ApplyPoincareMap(const State<ScalarType>& state, ScalarType mu, int section_index,
                                   ScalarType section_value, ScalarType max_time, ScalarType dt,
                                   ScalarType* period) {
  return ApplyPoincareMapImpl(state, mu, section_index, section_value, max_time, dt, period);
}

/**
 * @brief Newton-Raphson法による周期軌道の精密化
 */
template <typename ScalarType>
PeriodicOrbit<ScalarType> RefinePeriodicOrbit(const State<ScalarType>& initial_guess, ScalarType mu,
                                              int section_index, ScalarType section_value,
                                              int max_iterations, ScalarType tolerance,
                                              ScalarType max_time, ScalarType dt,
                                              NewtonConvergenceInfo<ScalarType>* convergence_info) {
  return RefinePeriodicOrbitImpl(initial_guess, mu, section_index, section_value, max_iterations,
                                 tolerance, max_time, dt, convergence_info);
}

/**
 * @brief モノドロミー行列の計算
 */
template <typename ScalarType>
std::array<std::array<ScalarType, 6>, 6> ComputeMonodromyMatrix(
    const PeriodicOrbit<ScalarType>& orbit, ScalarType mu, ScalarType dt) {
  return ComputeMonodromyMatrixImpl(orbit, mu, dt);
}

/**
 * @brief 安定性解析
 */
template <typename ScalarType>
void AnalyzeStability(PeriodicOrbit<ScalarType>* orbit, ScalarType eigenvalue_threshold) {
  AnalyzeStabilityImpl(orbit, eigenvalue_threshold);
}

/**
 * @brief 不変多様体の計算（詳細実装は後ほど追加）
 */
template <typename ScalarType>
std::vector<ManifoldTrajectory<ScalarType>> ComputeInvariantManifolds(
    const PeriodicOrbit<ScalarType>& orbit, ScalarType mu, const ManifoldConfig<ScalarType>& config,
    ScalarType dt);  // 実装は periodic_orbit_manifold.hpp に移動

/**
 * @brief 固有値計算（Eigenライブラリを使用）
 */
template <typename ScalarType>
std::vector<std::complex<ScalarType>> ComputeEigenvalues(
    const std::array<std::array<ScalarType, 6>, 6>& matrix);

/**
 * @brief 固有ベクトル計算（Eigenライブラリを使用）
 */
template <typename ScalarType>
std::array<ScalarType, 6> ComputeEigenvector(const std::array<std::array<ScalarType, 6>, 6>& matrix,
                                             const std::complex<ScalarType>& eigenvalue);

/**
 * @brief 固有値と固有ベクトルを同時に計算
 */
template <typename ScalarType>
struct EigenDecomposition {
  std::vector<std::complex<ScalarType>> eigenvalues;
  std::vector<std::array<std::complex<ScalarType>, 6>> eigenvectors;
};

template <typename ScalarType>
EigenDecomposition<ScalarType> ComputeEigenDecomposition(
    const std::array<std::array<ScalarType, 6>, 6>& matrix);

}  // namespace periodic_orbit

// 実装部をインクルード
#include "periodic_orbit_impl.hpp"
#include "periodic_orbit_manifold.hpp"

#endif  // PERIODIC_ORBIT_HPP
