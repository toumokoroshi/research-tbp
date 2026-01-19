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

#include "periodic_orbit_stability.hpp"

namespace periodic_orbit {

using namespace my_type;
using namespace crtbp;

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
template <typename ScalarType>
State<ScalarType> ApplyPoincareMapImpl(const State<ScalarType>& state, ScalarType mu,
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
  ScalarType next_progress_print = 10.0;  // Debug: print progress every 10 time units

  // Debug: print initial state
  std::cout << "\n        [Poincare] Starting: state=(" << state.x << ", " << state.y << ", "
            << state.z << ", " << state.vx << ", " << state.vy << ", " << state.vz
            << "), max_t=" << max_time << std::flush;
  std::cout << std::endl;

  while (time < max_time) {
    // Debug progress output
    if (time > next_progress_print) {
      std::cout << "        [Poincare] t=" << time << ", y=" << current_state.y
                << ", crossings=" << same_direction_count << std::endl;
      next_progress_print += 10.0;
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
      // クロッシング方向を判定
      int crossing_direction = (curr_val > prev_val) ? +1 : -1;

      if (first_crossing_direction == 0) {
        // 最初のクロッシング：方向を記録
        first_crossing_direction = crossing_direction;
        same_direction_count = 1;
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

  throw std::runtime_error(
      "ApplyPoincareMap: No second crossing found (orbit may not be periodic)");
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

  std::cout << "    [DEBUG] Starting Newton-Raphson iterations...\n";

  ScalarType prev_residual = std::numeric_limits<ScalarType>::max();
  int stagnation_count = 0;
  const ScalarType stagnation_tolerance =
      tolerance * 1000.0;  // Allow 1000x tolerance for stagnation (1e-9)

  for (iteration = 0; iteration < max_iterations; ++iteration) {
    auto iter_start = std::chrono::high_resolution_clock::now();

    std::cout << "      Iteration " << (iteration + 1) << "/" << max_iterations << ": ";
    std::cout.flush();

    // ポアンカレ写像を適用
    auto poincare_start = std::chrono::high_resolution_clock::now();
    State<ScalarType> current = vec_to_state(x);
    ScalarType period_out = 0.0;
    State<ScalarType> mapped_state;

    try {
      mapped_state = ApplyPoincareMapImpl(current, mu, section_index, section_value, max_time, dt,
                                          &period_out);
    } catch (const std::exception& e) {
      std::cout << "FAILED (Poincare map error)\n";
      if (convergence_info != nullptr) {
        convergence_info->converged = false;
        convergence_info->iterations = iteration;
        convergence_info->final_residual = std::numeric_limits<ScalarType>::infinity();
        convergence_info->message = std::string("Poincare map failed: ") + e.what();
      }
      throw;
    }

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
    if (converged || (stagnation_count >= 5 && residual_norm < stagnation_tolerance)) {
      std::cout << " CONVERGED!\n";
      converged = true;
      break;
    }

    // ヤコビアン行列を数値微分で計算
    auto jacobian_start = std::chrono::high_resolution_clock::now();
    std::cout << ", computing Jacobian...";
    std::cout.flush();

    std::array<std::array<ScalarType, FREE_DIM>, FREE_DIM> jacobian;
    ScalarType epsilon = 1e-8;

    for (int j = 0; j < FREE_DIM; ++j) {
      std::array<ScalarType, FREE_DIM> x_pert = x;
      x_pert[j] += epsilon;

      State<ScalarType> pert_state = vec_to_state(x_pert);
      State<ScalarType> pert_mapped;

      try {
        ScalarType dummy_period;
        pert_mapped = ApplyPoincareMapImpl(pert_state, mu, section_index, section_value, max_time,
                                           dt, &dummy_period);
      } catch (const std::exception&) {
        // 摂動が大きすぎる場合、εを小さくして再試行
        epsilon *= 0.1;
        x_pert[j] = x[j] + epsilon;
        pert_state = vec_to_state(x_pert);
        ScalarType dummy_period;
        pert_mapped = ApplyPoincareMapImpl(pert_state, mu, section_index, section_value, max_time,
                                           dt, &dummy_period);
      }

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

    // 適応的ダンピング: 残差が増加した場合はステップを縮小
    // 発散を防ぎながら収束時には高速化
    ScalarType alpha = 1.0;
    if (residual_norm > prev_residual * 1.1) {
      // 残差が10%以上増加した場合、ステップを半減
      alpha = 0.3;
    } else if (residual_norm > prev_residual) {
      // 若干増加した場合
      alpha = 0.5;
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

}  // namespace periodic_orbit

#endif  // PERIODIC_ORBIT_IMPL_HPP
