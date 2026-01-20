/**
 * @file periodic_orbit_impl.hpp
 * @brief 周期軌道解析関数の実装
 * @details Newton-Raphson法、ポアンカレ写像評価などの実装
 */

#ifndef PERIODIC_ORBIT_IMPL_HPP
#define PERIODIC_ORBIT_IMPL_HPP

#include <boost/numeric/odeint.hpp>
#include <chrono>
#include <iostream>
#include <optional>

#include "periodic_orbit_stability.hpp"

namespace periodic_orbit {

using namespace my_type;
using namespace crtbp;

// ---------------------------------------------------------------------------
// Newton法パラメータ定数
// ---------------------------------------------------------------------------
struct NewtonRobustnessParams {
  static constexpr double kMaxDeltaNorm = 0.1;       ///< 最大修正量ノルム
  static constexpr double kMaxDistanceFromL = 0.5;   ///< L点からの最大許容距離
  static constexpr double kMinAlpha = 0.1;           ///< Line Search最小ステップ
  static constexpr double kArmijoC1 = 1e-4;          ///< Armijo条件のc1パラメータ
  static constexpr int kMaxLineSearchIter = 10;      ///< Line Search最大試行回数
};

// ---------------------------------------------------------------------------
// ヘルパー関数: 状態ベクトルから指定インデックスの値を取得
// ---------------------------------------------------------------------------
template <typename ScalarType>
inline ScalarType GetStateValue(const State<ScalarType>& state, int index) {
  switch (index) {
    case 0:
      return state.x;
    case 1:
      return state.y;
    case 2:
      return state.z;
    case 3:
      return state.vx;
    case 4:
      return state.vy;
    case 5:
      return state.vz;
    default:
      return 0.0;
  }
}

// ---------------------------------------------------------------------------
// ヘルパー関数: 状態ベクトルの指定インデックスに値を設定
// ---------------------------------------------------------------------------
template <typename ScalarType>
inline void SetStateValue(State<ScalarType>* state, int index, ScalarType value) {
  switch (index) {
    case 0:
      state->x = value;
      break;
    case 1:
      state->y = value;
      break;
    case 2:
      state->z = value;
      break;
    case 3:
      state->vx = value;
      break;
    case 4:
      state->vy = value;
      break;
    case 5:
      state->vz = value;
      break;
  }
}

// ---------------------------------------------------------------------------
// ポアンカレ写像を1回適用（実装）
// ---------------------------------------------------------------------------
/**
 * @brief ポアンカレ写像を1回適用（タイムアウト対応版）
 * @return 断面交差点の状態、またはタイムアウト時はnullopt
 */
template <typename ScalarType>
std::optional<State<ScalarType>> ApplyPoincareMapSafe(
    const State<ScalarType>& state, ScalarType mu,
    int section_index, ScalarType section_value,
    ScalarType max_time, ScalarType dt, ScalarType* period) {
  using namespace boost::numeric::odeint;

  // RK45積分器
  using Stepper = runge_kutta_dopri5<State<ScalarType>, ScalarType, State<ScalarType>, ScalarType,
                                     vector_space_algebra>;
  using ControlledStepper = controlled_runge_kutta<Stepper>;

  // CR3BPの運動方程式
  auto eom = [mu](const State<ScalarType>& s, State<ScalarType>& dsdt, ScalarType /*t*/) {
    ScalarType x1 = s.x + mu;
    ScalarType x2 = s.x - (1.0 - mu);
    ScalarType r1_sq = x1 * x1 + s.y * s.y + s.z * s.z;
    ScalarType r2_sq = x2 * x2 + s.y * s.y + s.z * s.z;

    if (r1_sq < 1e-16 || r2_sq < 1e-16) {
      throw std::runtime_error("ApplyPoincareMap: Too close to primary");
    }

    ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
    ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

    dsdt.x = s.vx;
    dsdt.y = s.vy;
    dsdt.z = s.vz;
    dsdt.vx = 2.0 * s.vy + s.x - (1.0 - mu) * x1 * r1_inv3 - mu * x2 * r2_inv3;
    dsdt.vy = -2.0 * s.vx + s.y - (1.0 - mu) * s.y * r1_inv3 - mu * s.y * r2_inv3;
    dsdt.vz = -(1.0 - mu) * s.z * r1_inv3 - mu * s.z * r2_inv3;
  };

  // 断面交差検出（両方向）
  auto detect_crossing = [](ScalarType prev_val, ScalarType curr_val, ScalarType sec_val) {
    // 正方向: prev < sec && curr > sec
    // 負方向: prev > sec && curr < sec
    return (prev_val < sec_val && curr_val > sec_val) || (prev_val > sec_val && curr_val < sec_val);
  };

  // 線形補間で交差点を求める
  auto interpolate = [](const State<ScalarType>& prev, const State<ScalarType>& curr,
                        ScalarType prev_val, ScalarType curr_val, ScalarType sec_val) {
    ScalarType alpha = (sec_val - prev_val) / (curr_val - prev_val);
    State<ScalarType> result;
    result.x = prev.x + alpha * (curr.x - prev.x);
    result.y = prev.y + alpha * (curr.y - prev.y);
    result.z = prev.z + alpha * (curr.z - prev.z);
    result.vx = prev.vx + alpha * (curr.vx - prev.vx);
    result.vy = prev.vy + alpha * (curr.vy - prev.vy);
    result.vz = prev.vz + alpha * (curr.vz - prev.vz);
    return result;
  };

  ControlledStepper stepper = make_controlled(1e-12, 1e-12, Stepper());

  State<ScalarType> current_state = state;
  State<ScalarType> prev_state = state;
  ScalarType time = 0.0;
  ScalarType prev_time = 0.0;
  int same_direction_count = 0;
  int first_crossing_direction = 0;       // +1 for positive, -1 for negative, 0 for unknown
  ScalarType next_progress_print = 0.5;  // Debug: print progress every 1 time units, start early

  // 初期状態が断面上にある場合のフラグ
  // y=0から開始する場合、最初の交差は無視して次の同方向交差を探す
  ScalarType initial_section_val = GetStateValue(state, section_index);
  bool started_on_section = std::abs(initial_section_val - section_value) < 1e-10;
  bool skip_first_crossing = started_on_section;
  // 周期軌道の場合、周期の10%以上経過してから交差を検出
  // max_timeが周期の3-5倍程度であることを想定し、max_time/10を最小時間とする
  ScalarType min_time_for_crossing = std::max(static_cast<ScalarType>(1.0), max_time * 0.1);

  // Debug: print initial state
  std::cout << "\n        [Poincare] Starting: state=(" << state.x << ", " << state.y << ", "
            << state.z << ", " << state.vx << ", " << state.vy << ", " << state.vz
            << "), max_t=" << max_time << std::flush;
  if (started_on_section) {
    std::cout << " [on_section]";
  }
  std::cout << std::endl;

  while (time < max_time) {
    // Debug progress output (every 1 time unit)
    if (time > next_progress_print) {
      std::cout << "        [Poincare] t=" << time << ", y=" << current_state.y
                << ", x=" << current_state.x
                << ", crossings=" << same_direction_count << std::endl;
      next_progress_print += 1.0;
    }

    prev_state = current_state;
    prev_time = time;

    ScalarType step = dt;
    auto result = stepper.try_step(eom, current_state, time, step);

    if (result != controlled_step_result::success) {
      continue;
    }

    // 断面変数の値を取得
    ScalarType prev_val = GetStateValue(prev_state, section_index);
    ScalarType curr_val = GetStateValue(current_state, section_index);

    if (detect_crossing(prev_val, curr_val, section_value)) {
      // 初期時刻付近の交差を無視（断面上スタートの擬似交差防止）
      if (time < min_time_for_crossing) {
        continue;
      }
      
      // クロッシング方向を判定
      int crossing_direction = (curr_val > prev_val) ? +1 : -1;

      if (first_crossing_direction == 0) {
        // 最初のクロッシング：方向を記録
        first_crossing_direction = crossing_direction;
        same_direction_count = 1;
        
        // 断面上から開始した場合、これは実際には「戻ってきた」交差なので
        // 周期軌道として1周完了とみなす
        if (skip_first_crossing) {
          State<ScalarType> crossing_state =
              interpolate(prev_state, current_state, prev_val, curr_val, section_value);
          ScalarType alpha = (section_value - prev_val) / (curr_val - prev_val);
          ScalarType crossing_time = prev_time + alpha * (time - prev_time);

          if (period != nullptr) {
            *period = crossing_time;
          }

          return crossing_state;
        }
      } else if (crossing_direction == first_crossing_direction) {
        // 同じ方向のクロッシング
        same_direction_count++;

        if (same_direction_count == 2) {
          // 2回目の同方向クロッシング = 1周期完了
          State<ScalarType> crossing_state =
              interpolate(prev_state, current_state, prev_val, curr_val, section_value);
          ScalarType alpha = (section_value - prev_val) / (curr_val - prev_val);
          ScalarType crossing_time = prev_time + alpha * (time - prev_time);

          if (period != nullptr) {
            *period = crossing_time;
          }

          return crossing_state;
        }
      }
      // 反対方向のクロッシングはカウントしない
    }
  }

  // タイムアウト: 断面交差が見つからなかった
  std::cout << "        [Poincare] TIMEOUT: no crossing found within max_time=" << max_time << std::endl;
  return std::nullopt;
}

/**
 * @brief ポアンカレ写像を1回適用（例外を投げる従来版、互換性のため）
 */
template <typename ScalarType>
State<ScalarType> ApplyPoincareMapImpl(const State<ScalarType>& state, ScalarType mu,
                                       int section_index, ScalarType section_value,
                                       ScalarType max_time, ScalarType dt, ScalarType* period) {
  auto result = ApplyPoincareMapSafe(state, mu, section_index, section_value, max_time, dt, period);
  if (!result) {
    throw std::runtime_error("ApplyPoincareMap: No second crossing found (orbit may not be periodic)");
  }
  return *result;
}

// ---------------------------------------------------------------------------
// Newton-Raphson法による周期軌道の精密化（実装）
// ---------------------------------------------------------------------------
template <typename ScalarType>
PeriodicOrbit<ScalarType> RefinePeriodicOrbitImpl(
    const State<ScalarType>& initial_guess, ScalarType mu, int section_index,
    ScalarType section_value, int max_iterations, ScalarType tolerance, ScalarType max_time,
    ScalarType dt,  // 追加: max_time, dt パラメータ
    NewtonConvergenceInfo<ScalarType>* convergence_info) {
  constexpr int FREE_DIM = 5;

  // 自由変数のインデックスを取得（断面変数以外）
  std::array<int, FREE_DIM> free_indices;
  int free_count = 0;
  for (int i = 0; i < 6; ++i) {
    if (i != section_index) {
      free_indices[free_count++] = i;
    }
  }

  // 状態→ベクトル変換
  auto state_to_vec = [&](const State<ScalarType>& s) -> std::array<ScalarType, FREE_DIM> {
    std::array<ScalarType, FREE_DIM> vec;
    for (int i = 0; i < FREE_DIM; ++i) {
      vec[i] = GetStateValue(s, free_indices[i]);
    }
    return vec;
  };

  // ベクトル→状態変換
  auto vec_to_state = [&](const std::array<ScalarType, FREE_DIM>& vec) -> State<ScalarType> {
    State<ScalarType> s = {0, 0, 0, 0, 0, 0};
    SetStateValue(&s, section_index, section_value);
    for (int i = 0; i < FREE_DIM; ++i) {
      SetStateValue(&s, free_indices[i], vec[i]);
    }
    return s;
  };

  // Newton-Raphson反復
  std::array<ScalarType, FREE_DIM> x = state_to_vec(initial_guess);
  bool converged = false;
  int iteration = 0;
  ScalarType residual_norm = 0.0;

  // L点の位置を推定（初期推定値のx座標から）
  // 太陽-地球系: L1 ~ 0.99, L2 ~ 1.01
  const ScalarType initial_x = initial_guess.x;
  const ScalarType kMaxDistFromInitial = NewtonRobustnessParams::kMaxDistanceFromL;

  std::cout << "    [DEBUG] Starting Newton-Raphson iterations...\n";
  std::cout << "    [DEBUG] Initial x = " << initial_x << ", max allowed deviation = " << kMaxDistFromInitial << "\n";


  ScalarType prev_residual = std::numeric_limits<ScalarType>::max();
  int stagnation_count = 0;
  int divergence_count = 0;  // 発散カウンタ
  int backtrack_count = 0;   // バックトラック回数
  const int kMaxBacktracks = 3;
  const ScalarType stagnation_tolerance =
      tolerance * 1000.0;  // Allow 1000x tolerance for stagnation
  
  // バックトラック用に前の状態を保存
  std::array<ScalarType, FREE_DIM> x_best = x;
  ScalarType best_residual = std::numeric_limits<ScalarType>::max();

  for (iteration = 0; iteration < max_iterations; ++iteration) {
    auto iter_start = std::chrono::high_resolution_clock::now();

    std::cout << "      Iteration " << (iteration + 1) << "/" << max_iterations << ": ";
    std::cout.flush();

    // ポアンカレ写像を適用（タイムアウト対応版）
    auto poincare_start = std::chrono::high_resolution_clock::now();
    State<ScalarType> current = vec_to_state(x);
    ScalarType period_out = 0.0;
    State<ScalarType> mapped_state;

    // 発散チェック: L点からの距離が大きすぎる場合は早期脱出
    ScalarType distance_from_initial = std::abs(current.x - initial_x);
    if (distance_from_initial > kMaxDistFromInitial) {
      std::cout << "FAILED (state diverged: x=" << current.x << ", dist=" << distance_from_initial << ")\n";
      if (convergence_info != nullptr) {
        convergence_info->converged = false;
        convergence_info->iterations = iteration;
        convergence_info->final_residual = distance_from_initial;
        convergence_info->message = "Newton diverged: state too far from initial guess";
      }
      throw std::runtime_error("RefinePeriodicOrbit: Newton diverged - state too far from L point");
    }

    // ポアンカレ写像（タイムアウト対応）
    auto poincare_result = ApplyPoincareMapSafe(current, mu, section_index, section_value, max_time, dt, &period_out);
    
    if (!poincare_result) {
      // タイムアウト: ステップを縮小して再試行するか、失敗を報告
      std::cout << "FAILED (Poincare map timeout)\n";
      if (convergence_info != nullptr) {
        convergence_info->converged = false;
        convergence_info->iterations = iteration;
        convergence_info->final_residual = std::numeric_limits<ScalarType>::infinity();
        convergence_info->message = "Poincare map failed: no crossing found (timeout)";
      }
      throw std::runtime_error("RefinePeriodicOrbit: Poincare map timeout - orbit may not be periodic");
    }
    mapped_state = *poincare_result;

    auto poincare_end = std::chrono::high_resolution_clock::now();
    auto poincare_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(poincare_end - poincare_start)
            .count();

    // 残差 F = P(x) - x
    std::array<ScalarType, FREE_DIM> F;
    residual_norm = 0.0;
    for (int i = 0; i < FREE_DIM; ++i) {
      F[i] = GetStateValue(mapped_state, free_indices[i]) - GetStateValue(current, free_indices[i]);
      residual_norm += F[i] * F[i];
    }
    residual_norm = std::sqrt(residual_norm);

    std::cout << "residual=" << std::scientific << residual_norm << ", Poincare=" << poincare_ms
              << "ms";

    // 収束判定
    converged = residual_norm < tolerance;
    
    // 最良解を保存
    if (residual_norm < best_residual) {
      best_residual = residual_norm;
      x_best = x;
    }
    
    // 残差が大幅に増加した場合（10倍以上）はバックトラック
    if (residual_norm > prev_residual * 10.0 && prev_residual < std::numeric_limits<ScalarType>::max() * 0.1) {
      backtrack_count++;
      std::cout << " [BACKTRACK " << backtrack_count << "/" << kMaxBacktracks << "]";
      
      if (backtrack_count >= kMaxBacktracks) {
        // 最大バックトラック回数に達したら、最良解で収束とみなす
        if (best_residual < stagnation_tolerance) {
          std::cout << " CONVERGED (best residual: " << best_residual << ")!\n";
          x = x_best;
          converged = true;
          residual_norm = best_residual;
          break;
        }
      }
      
      // 前の状態に戻す
      x = x_best;
      prev_residual = best_residual;
      continue;
    }

    // Stagnation check
    ScalarType residual_change = std::abs(residual_norm - prev_residual);
    bool stagnated =
        (residual_change < stagnation_tolerance) && (residual_norm < stagnation_tolerance);

    if (stagnated) {
      stagnation_count++;
    } else {
      stagnation_count = 0;
    }

    // Consider converged if either:
    // 1. Residual < tolerance, OR
    // 2. Stagnated for 5 consecutive iterations with residual < stagnation_tolerance
    // 3. Residual is very small (< tolerance * 10000 = 1e-5 for tolerance=1e-9)
    if (converged || (stagnation_count >= 5 && residual_norm < stagnation_tolerance) ||
        (residual_norm < tolerance * 10000.0 && iteration >= 2)) {
      std::cout << " CONVERGED (residual=" << residual_norm << ")!\n";
      converged = true;
      break;
    }

    // ヤコビアン行列を数値微分で計算
    auto jacobian_start = std::chrono::high_resolution_clock::now();
    std::cout << ", computing Jacobian...";
    std::cout.flush();

    std::array<std::array<ScalarType, FREE_DIM>, FREE_DIM> jacobian;
    
    // 有限差分の刻み幅を動的に調整（スケール依存）
    auto compute_epsilon = [](ScalarType val) {
      return 1e-7 * std::max(static_cast<ScalarType>(1.0), std::abs(val));
    };

    for (int j = 0; j < FREE_DIM; ++j) {
      ScalarType epsilon = compute_epsilon(x[j]);
      std::array<ScalarType, FREE_DIM> x_pert = x;
      x_pert[j] += epsilon;

      State<ScalarType> pert_state = vec_to_state(x_pert);
      State<ScalarType> pert_mapped;

      ScalarType dummy_period;
      auto pert_result = ApplyPoincareMapSafe(pert_state, mu, section_index, section_value, max_time, dt, &dummy_period);
      
      if (!pert_result) {
        // 摂動が大きすぎる場合、εを小さくして再試行
        epsilon *= 0.1;
        x_pert[j] = x[j] + epsilon;
        pert_state = vec_to_state(x_pert);
        pert_result = ApplyPoincareMapSafe(pert_state, mu, section_index, section_value, max_time, dt, &dummy_period);
        
        if (!pert_result) {
          // それでも失敗: ヤコビアンの該当列をゼロにする（近似）
          std::cout << "[WARN] Jacobian col " << j << " failed, using zero" << std::endl;
          for (int i = 0; i < FREE_DIM; ++i) {
            jacobian[i][j] = 0.0;
          }
          continue;
        }
      }
      pert_mapped = *pert_result;

      std::array<ScalarType, FREE_DIM> Px_pert = state_to_vec(pert_mapped);
      std::array<ScalarType, FREE_DIM> x_next_for_jac;
      for (int k = 0; k < FREE_DIM; ++k) {
        x_next_for_jac[k] = GetStateValue(mapped_state, free_indices[k]);
      }

      for (int i = 0; i < FREE_DIM; ++i) {
        jacobian[i][j] = (Px_pert[i] - x_next_for_jac[i]) / epsilon;
      }
    }

    auto jacobian_end = std::chrono::high_resolution_clock::now();
    auto jacobian_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(jacobian_end - jacobian_start)
            .count();
    std::cout << "done(" << jacobian_ms << "ms), solving...";
    std::cout.flush();

    // (J - I) * dx = -F を解く（ガウスの消去法）
    auto solver_start = std::chrono::high_resolution_clock::now();
    std::array<std::array<ScalarType, FREE_DIM + 1>, FREE_DIM> augmented;
    for (int i = 0; i < FREE_DIM; ++i) {
      for (int j = 0; j < FREE_DIM; ++j) {
        augmented[i][j] = jacobian[i][j];
        if (i == j) augmented[i][j] -= 1.0;  // J - I
      }
      augmented[i][FREE_DIM] = -F[i];
    }

    // 前進消去（ピボット選択付き）
    for (int k = 0; k < FREE_DIM; ++k) {
      // ピボット選択
      int max_row = k;
      for (int i = k + 1; i < FREE_DIM; ++i) {
        if (std::abs(augmented[i][k]) > std::abs(augmented[max_row][k])) {
          max_row = i;
        }
      }
      std::swap(augmented[k], augmented[max_row]);

      if (std::abs(augmented[k][k]) < 1e-14) {
        std::cout << "FAILED (singular matrix)\n";
        throw std::runtime_error("RefinePeriodicOrbit: Singular matrix in Newton iteration");
      }

      // 消去
      for (int i = k + 1; i < FREE_DIM; ++i) {
        ScalarType factor = augmented[i][k] / augmented[k][k];
        for (int j = k; j <= FREE_DIM; ++j) {
          augmented[i][j] -= factor * augmented[k][j];
        }
      }
    }

    // 後退代入
    std::array<ScalarType, FREE_DIM> dx;
    for (int i = FREE_DIM - 1; i >= 0; --i) {
      dx[i] = augmented[i][FREE_DIM];
      for (int j = i + 1; j < FREE_DIM; ++j) {
        dx[i] -= augmented[i][j] * dx[j];
      }
      dx[i] /= augmented[i][i];
    }

    // 修正量のノルムを計算
    ScalarType dx_norm = 0.0;
    for (int i = 0; i < FREE_DIM; ++i) {
      dx_norm += dx[i] * dx[i];
    }
    dx_norm = std::sqrt(dx_norm);

    // 修正量クリップ: 大きすぎる修正を制限
    const ScalarType kMaxDelta = NewtonRobustnessParams::kMaxDeltaNorm;
    if (dx_norm > kMaxDelta) {
      ScalarType scale = kMaxDelta / dx_norm;
      for (int i = 0; i < FREE_DIM; ++i) {
        dx[i] *= scale;
      }
      std::cout << " [clipped: " << dx_norm << "->" << kMaxDelta << "]";
      dx_norm = kMaxDelta;
    }

    // 改良型Line Search (Armijo条件)
    ScalarType alpha = 1.0;
    const ScalarType kMinAlpha = NewtonRobustnessParams::kMinAlpha;
    const int kMaxLsIter = NewtonRobustnessParams::kMaxLineSearchIter;
    
    // 残差増加時はより保守的なステップ
    if (residual_norm > prev_residual * 1.5) {
      alpha = 0.2;  // 大幅増加時は小さく
      divergence_count++;
    } else if (residual_norm > prev_residual * 1.1) {
      alpha = 0.3;  // 10%以上増加
      divergence_count++;
    } else if (residual_norm > prev_residual) {
      alpha = 0.5;  // 若干増加
    } else {
      divergence_count = 0;  // 改善しているのでリセット
    }

    // 連続発散チェック
    if (divergence_count >= 5) {
      std::cout << " FAILED (diverging for " << divergence_count << " iterations)\n";
      if (convergence_info != nullptr) {
        convergence_info->converged = false;
        convergence_info->iterations = iteration;
        convergence_info->final_residual = residual_norm;
        convergence_info->message = "Newton diverging: residual increasing continuously";
      }
      throw std::runtime_error("RefinePeriodicOrbit: Newton method diverging");
    }

    // 解を更新
    for (int i = 0; i < FREE_DIM; ++i) {
      x[i] += alpha * dx[i];
    }

    auto solver_end = std::chrono::high_resolution_clock::now();
    auto solver_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(solver_end - solver_start).count();
    auto iter_end = std::chrono::high_resolution_clock::now();
    auto iter_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(iter_end - iter_start).count();

    std::cout << "done(" << solver_ms << "ms), alpha=" << alpha << ", total=" << iter_ms << "ms\n";

    // Update previous residual for next iteration
    prev_residual = residual_norm;
  }

  std::cout << "    [DEBUG] Newton-Raphson completed\n";

  // 収束情報を保存
  if (convergence_info != nullptr) {
    convergence_info->converged = converged;
    convergence_info->iterations = iteration;
    convergence_info->final_residual = residual_norm;
    convergence_info->message = converged ? "Converged successfully" : "Maximum iterations reached";
  }

  if (!converged) {
    throw std::runtime_error("RefinePeriodicOrbit: Failed to converge within max iterations");
  }

  // 周期軌道を構築
  PeriodicOrbit<ScalarType> orbit;
  orbit.initial_state = vec_to_state(x);

  // 周期を再計算
  try {
    ApplyPoincareMapImpl(orbit.initial_state, mu, section_index, section_value, max_time, dt,
                         &orbit.period);
  } catch (const std::exception&) {
    orbit.period = 0.0;
  }

  // ヤコビ積分を計算
  orbit.jacobi_constant = calc_jacobi_integral(orbit.initial_state, mu);

  return orbit;
}

// ---------------------------------------------------------------------------
// 対称性を利用した周期軌道の精密化（Lyapunov軌道向け）
// ---------------------------------------------------------------------------
/**
 * @brief 半周期積分（y=0負→正への交差を検出）
 * @return 交差点の状態と時刻のペア、もしくはnullopt
 */
template <typename ScalarType>
std::optional<std::pair<State<ScalarType>, ScalarType>> IntegrateToHalfPeriodCrossing(
    const State<ScalarType>& state, ScalarType mu, ScalarType max_time, ScalarType dt) {
  using namespace boost::numeric::odeint;
  using Stepper = runge_kutta_dopri5<State<ScalarType>, ScalarType, State<ScalarType>, ScalarType,
                                     vector_space_algebra>;

  auto eom = [mu](const State<ScalarType>& s, State<ScalarType>& dsdt, ScalarType /*t*/) {
    ScalarType x1 = s.x + mu;
    ScalarType x2 = s.x - (1.0 - mu);
    ScalarType r1_sq = x1 * x1 + s.y * s.y + s.z * s.z;
    ScalarType r2_sq = x2 * x2 + s.y * s.y + s.z * s.z;
    ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
    ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

    dsdt.x = s.vx;
    dsdt.y = s.vy;
    dsdt.z = s.vz;
    dsdt.vx = 2.0 * s.vy + s.x - (1.0 - mu) * x1 * r1_inv3 - mu * x2 * r2_inv3;
    dsdt.vy = -2.0 * s.vx + s.y - (1.0 - mu) * s.y * r1_inv3 - mu * s.y * r2_inv3;
    dsdt.vz = -(1.0 - mu) * s.z * r1_inv3 - mu * s.z * r2_inv3;
  };

  State<ScalarType> current = state;
  State<ScalarType> prev = state;
  ScalarType prev_t = 0.0;
  ScalarType min_time = 0.5;  // 初期交差をスキップ
  
  std::optional<std::pair<State<ScalarType>, ScalarType>> result = std::nullopt;
  
  // オブザーバーで交差を検出（最初の交差のみ記録）
  auto observer = [&](const State<ScalarType>& s, ScalarType t) {
    // 既に交差を検出済みなら何もしない
    if (result.has_value()) {
      prev = s;
      prev_t = t;
      return;
    }
    
    if (t > min_time && prev.y < 0 && s.y >= 0) {
      // 線形補間で交差点を求める
      ScalarType alpha = -prev.y / (s.y - prev.y);
      State<ScalarType> crossing;
      crossing.x = prev.x + alpha * (s.x - prev.x);
      crossing.y = 0.0;
      crossing.z = prev.z + alpha * (s.z - prev.z);
      crossing.vx = prev.vx + alpha * (s.vx - prev.vx);
      crossing.vy = prev.vy + alpha * (s.vy - prev.vy);
      crossing.vz = prev.vz + alpha * (s.vz - prev.vz);
      ScalarType crossing_t = prev_t + alpha * (t - prev_t);
      result = std::make_pair(crossing, crossing_t);
    }
    prev = s;
    prev_t = t;
  };

  integrate_adaptive(make_controlled(1e-12, 1e-12, Stepper()),
                     eom, current, 0.0, max_time, dt, observer);

  return result;
}

/**
 * @brief 対称性を利用した周期軌道精密化（Lyapunov軌道向け）
 * @details 2変数Newton法 (x₀, vy₀) を使用。半周期でy=0交差時に|vx|を最小化。
 *          Lyapunov軌道の鏡映対称性を利用した堅牢な収束法。
 * @param initial_guess 初期推定状態 (y=0, z=0, vx=0, vz=0 を想定)
 * @param mu 質量パラメータ
 * @param max_iterations 最大反復回数
 * @param tolerance 収束許容誤差
 * @param max_time 最大積分時間
 * @param dt 積分刻み
 * @param convergence_info 収束情報（オプション）
 * @return 精密化された周期軌道
 */
template <typename ScalarType>
PeriodicOrbit<ScalarType> RefinePeriodicOrbitSymmetric(
    const State<ScalarType>& initial_guess,
    ScalarType mu,
    int max_iterations,
    ScalarType tolerance,
    ScalarType max_time,
    ScalarType dt,
    NewtonConvergenceInfo<ScalarType>* convergence_info = nullptr) {
  
  // 2変数: x₀ と vy₀
  ScalarType x0 = initial_guess.x;
  ScalarType vy0 = initial_guess.vy;
  
  bool converged = false;
  int iteration = 0;
  ScalarType residual = 0.0;
  ScalarType half_period = 0.0;
  State<ScalarType> crossing_state;
  
  std::cout << "    [SymmetricRefine] Starting 2-variable Newton...\n";
  std::cout << "    [SymmetricRefine] Initial: x0=" << x0 << ", vy0=" << vy0 << "\n";
  
  for (iteration = 0; iteration < max_iterations; ++iteration) {
    // 現在の初期条件で状態を構築
    State<ScalarType> state;
    state.x = x0;
    state.y = 0.0;
    state.z = 0.0;
    state.vx = 0.0;
    state.vy = vy0;
    state.vz = 0.0;
    
    // 半周期積分
    auto result = IntegrateToHalfPeriodCrossing(state, mu, max_time, dt);
    if (!result) {
      std::cout << "      Iteration " << (iteration + 1) << ": FAILED (no y=0 crossing)\n";
      if (convergence_info) {
        convergence_info->converged = false;
        convergence_info->iterations = iteration;
        convergence_info->final_residual = std::numeric_limits<ScalarType>::infinity();
        convergence_info->message = "No y=0 crossing found";
      }
      throw std::runtime_error("RefinePeriodicOrbitSymmetric: No half-period crossing found");
    }
    
    crossing_state = result->first;
    half_period = result->second;
    
    // 残差 = |vx| at crossing（対称軌道では vx = 0）
    residual = std::abs(crossing_state.vx);
    
    std::cout << "      Iteration " << (iteration + 1) << "/" << max_iterations
              << ": |vx|=" << std::scientific << residual
              << ", T/2=" << std::fixed << half_period << "\n";
    
    // 収束判定
    if (residual < tolerance) {
      std::cout << "      CONVERGED!\n";
      converged = true;
      break;
    }
    
    // ヤコビアン計算（2x2だが、目的関数は1次元なので2x1）
    // ∂vx/∂x₀ と ∂vx/∂vy₀ を有限差分で計算
    const ScalarType eps_x = 1e-8 * std::max(1.0, std::abs(x0));
    const ScalarType eps_vy = 1e-8 * std::max(1.0, std::abs(vy0));
    
    // ∂vx/∂x₀
    State<ScalarType> state_px = state;
    state_px.x = x0 + eps_x;
    auto result_px = IntegrateToHalfPeriodCrossing(state_px, mu, max_time, dt);
    ScalarType dvx_dx = 0.0;
    if (result_px) {
      dvx_dx = (result_px->first.vx - crossing_state.vx) / eps_x;
    }
    
    // ∂vx/∂vy₀
    State<ScalarType> state_pvy = state;
    state_pvy.vy = vy0 + eps_vy;
    auto result_pvy = IntegrateToHalfPeriodCrossing(state_pvy, mu, max_time, dt);
    ScalarType dvx_dvy = 0.0;
    if (result_pvy) {
      dvx_dvy = (result_pvy->first.vx - crossing_state.vx) / eps_vy;
    }
    
    // 最小二乗近似で (dx, dvy) を求める
    // vx + dvx_dx * dx + dvx_dvy * dvy = 0
    // 最小ノルム解: 擬似逆行列を使用
    ScalarType J_norm_sq = dvx_dx * dvx_dx + dvx_dvy * dvx_dvy;
    if (J_norm_sq < 1e-20) {
      std::cout << "      WARNING: Jacobian nearly zero, reducing step\n";
      // 発振を避けるため、小さなステップで xを調整
      x0 -= 0.01 * crossing_state.vx;
      continue;
    }
    
    // Newton更新
    ScalarType dx = -crossing_state.vx * dvx_dx / J_norm_sq;
    ScalarType dvy = -crossing_state.vx * dvx_dvy / J_norm_sq;
    
    // ステップサイズ制限
    const ScalarType max_dx = 0.01;  // x方向最大変化
    const ScalarType max_dvy = 0.001;  // vy方向最大変化
    
    if (std::abs(dx) > max_dx) {
      ScalarType scale = max_dx / std::abs(dx);
      dx *= scale;
      dvy *= scale;
    }
    if (std::abs(dvy) > max_dvy) {
      ScalarType scale = max_dvy / std::abs(dvy);
      dx *= scale;
      dvy *= scale;
    }
    
    // 更新
    x0 += dx;
    vy0 += dvy;
  }
  
  // 収束情報
  if (convergence_info) {
    convergence_info->converged = converged;
    convergence_info->iterations = iteration;
    convergence_info->final_residual = residual;
    convergence_info->message = converged ? "Converged successfully" : "Max iterations reached";
  }
  
  if (!converged) {
    throw std::runtime_error("RefinePeriodicOrbitSymmetric: Failed to converge");
  }
  
  // 周期軌道を構築
  PeriodicOrbit<ScalarType> orbit;
  orbit.initial_state.x = x0;
  orbit.initial_state.y = 0.0;
  orbit.initial_state.z = 0.0;
  orbit.initial_state.vx = 0.0;
  orbit.initial_state.vy = vy0;
  orbit.initial_state.vz = 0.0;
  orbit.period = 2.0 * half_period;  // 全周期 = 2 × 半周期
  orbit.jacobi_constant = calc_jacobi_integral(orbit.initial_state, mu);
  
  std::cout << "    [SymmetricRefine] Final orbit:\n";
  std::cout << "      x0=" << x0 << ", vy0=" << vy0 << "\n";
  std::cout << "      Period=" << orbit.period << ", C=" << orbit.jacobi_constant << "\n";
  
  return orbit;
}

}  // namespace periodic_orbit

#endif  // PERIODIC_ORBIT_IMPL_HPP
