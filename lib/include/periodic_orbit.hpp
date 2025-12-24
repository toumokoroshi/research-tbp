/**
 * @file periodic_orbit.hpp
 * @brief 周期軌道の検出・安定性解析・不変多様体計算のためのライブラリ
 * @details CR3BP（円制限三体問題）における周期軌道の解析ツール
 * @date 2025-12-24
 */

#ifndef PERIODIC_ORBIT_HPP
#define PERIODIC_ORBIT_HPP

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
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
 */
template <typename ScalarType>
struct PeriodicOrbit {
  State<ScalarType> initial_state;   ///< 初期状態 (x, y, z, vx, vy, vz)
  ScalarType period = 0.0;           ///< 周期 (無次元時間)
  ScalarType jacobi_constant = 0.0;  ///< ヤコビ積分の値

  // 安定性解析結果
  bool stability_computed = false;                    ///< 安定性解析が実行されたか
  std::vector<std::complex<ScalarType>> eigenvalues;  ///< モノドロミー行列の固有値
  bool is_stable = false;                             ///< 安定か不安定か (SPO/UPO)
  ScalarType stability_index = 0.0;                   ///< 安定性指数 (最大固有値の絶対値)

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
 */
template <typename ScalarType>
struct ManifoldConfig {
  ScalarType epsilon = 1e-6;          ///< 固有ベクトル方向への初期変位量
  ScalarType forward_time = 100.0;    ///< 正方向積分時間
  ScalarType backward_time = -100.0;  ///< 負方向積分時間
  int num_initial_points = 20;        ///< 各多様体の初期点数
  bool compute_stable = true;         ///< 安定多様体を計算するか
  bool compute_unstable = true;       ///< 不安定多様体を計算するか
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
    ScalarType dt) {
  // TODO: 実装予定
  // 固有ベクトル方向に変位を加えて時間反転積分
  throw std::runtime_error("ComputeInvariantManifolds: Not implemented yet");
}

/**
 * @brief 固有値計算（詳細実装は後ほど追加）
 */
template <typename ScalarType>
std::vector<std::complex<ScalarType>> ComputeEigenvalues(
    const std::array<std::array<ScalarType, 6>, 6>& matrix) {
  // TODO: 実装予定
  // QR法またはEigenライブラリを使用
  throw std::runtime_error("ComputeEigenvalues: Not implemented yet");
}

/**
 * @brief 固有ベクトル計算（詳細実装は後ほど追加）
 */
template <typename ScalarType>
std::array<ScalarType, 6> ComputeEigenvector(const std::array<std::array<ScalarType, 6>, 6>& matrix,
                                             const std::complex<ScalarType>& eigenvalue) {
  // TODO: 実装予定
  throw std::runtime_error("ComputeEigenvector: Not implemented yet");
}

}  // namespace periodic_orbit

// 実装部をインクルード
#include "periodic_orbit_impl.hpp"

#endif  // PERIODIC_ORBIT_HPP
