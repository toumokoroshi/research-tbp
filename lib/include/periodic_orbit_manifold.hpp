/**
 * @file periodic_orbit_manifold.hpp
 * @brief 不変多様体計算の実装
 * @details 固有値・固有ベクトル計算と安定/不安定多様体の計算
 * @date 2026-01-05
 */

#ifndef PERIODIC_ORBIT_MANIFOLD_HPP
#define PERIODIC_ORBIT_MANIFOLD_HPP

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/numeric/odeint.hpp>
#include <cmath>

#include "periodic_orbit_stability.hpp"

namespace periodic_orbit {

using namespace my_type;
using namespace crtbp;

// ---------------------------------------------------------------------------
// 固有値・固有ベクトル計算（Eigenライブラリを使用）
// ---------------------------------------------------------------------------

/**
 * @brief 固有値計算の実装
 */
template <typename ScalarType>
std::vector<std::complex<ScalarType>> ComputeEigenvalues(
    const std::array<std::array<ScalarType, 6>, 6>& matrix) {
  // std::arrayからEigen行列に変換
  Eigen::Matrix<ScalarType, 6, 6> M;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      M(i, j) = matrix[i][j];
    }
  }

  // 固有値計算
  Eigen::EigenSolver<Eigen::Matrix<ScalarType, 6, 6>> solver(M);
  auto eigenvalues = solver.eigenvalues();

  // 結果をstd::vectorに変換
  std::vector<std::complex<ScalarType>> result(6);
  for (int i = 0; i < 6; ++i) {
    result[i] = eigenvalues(i);
  }
  return result;
}

/**
 * @brief 固有値と固有ベクトルを同時に計算
 */
template <typename ScalarType>
EigenDecomposition<ScalarType> ComputeEigenDecomposition(
    const std::array<std::array<ScalarType, 6>, 6>& matrix) {
  // std::arrayからEigen行列に変換
  Eigen::Matrix<ScalarType, 6, 6> M;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      M(i, j) = matrix[i][j];
    }
  }

  // 固有値・固有ベクトル計算
  Eigen::EigenSolver<Eigen::Matrix<ScalarType, 6, 6>> solver(M);
  auto eigenvalues = solver.eigenvalues();
  auto eigenvectors = solver.eigenvectors();

  EigenDecomposition<ScalarType> result;
  result.eigenvalues.resize(6);
  result.eigenvectors.resize(6);

  for (int i = 0; i < 6; ++i) {
    result.eigenvalues[i] = eigenvalues(i);
    for (int j = 0; j < 6; ++j) {
      result.eigenvectors[i][j] = eigenvectors(j, i);
    }
  }
  return result;
}

/**
 * @brief 固有ベクトル計算の実装
 */
template <typename ScalarType>
std::array<ScalarType, 6> ComputeEigenvector(const std::array<std::array<ScalarType, 6>, 6>& matrix,
                                             const std::complex<ScalarType>& eigenvalue) {
  auto decomp = ComputeEigenDecomposition(matrix);

  // 指定された固有値に最も近い固有ベクトルを探す
  int best_idx = 0;
  ScalarType min_diff = std::abs(decomp.eigenvalues[0] - eigenvalue);
  for (int i = 1; i < 6; ++i) {
    ScalarType diff = std::abs(decomp.eigenvalues[i] - eigenvalue);
    if (diff < min_diff) {
      min_diff = diff;
      best_idx = i;
    }
  }

  // 実数部のみを抽出（実固有値の場合）
  std::array<ScalarType, 6> result;
  for (int i = 0; i < 6; ++i) {
    result[i] = decomp.eigenvectors[best_idx][i].real();
  }
  return result;
}

// ---------------------------------------------------------------------------
// Jacobi surface projection
// ---------------------------------------------------------------------------

template <typename ScalarType>
inline void ProjectToJacobiSurface(State<ScalarType>* state, ScalarType mu, ScalarType target_cj) {
  const ScalarType v_sq = state->vx * state->vx + state->vy * state->vy + state->vz * state->vz;
  if (v_sq <= static_cast<ScalarType>(0)) {
    throw std::runtime_error("ProjectToJacobiSurface: zero velocity magnitude");
  }

  const ScalarType U_star = calc_potential_U(state->x, state->y, state->z, mu);
  const ScalarType v_desired_sq = 2.0 * U_star - target_cj;
  if (v_desired_sq < static_cast<ScalarType>(0)) {
    throw std::runtime_error("ProjectToJacobiSurface: negative v^2 for target Cj");
  }

  const ScalarType scale = std::sqrt(v_desired_sq / v_sq);
  state->vx *= scale;
  state->vy *= scale;
  state->vz *= scale;
}

// ---------------------------------------------------------------------------
// 安定多様体計算
// ---------------------------------------------------------------------------

/**
 * @brief 安定・不安定多様体を計算
 */
template <typename ScalarType>
std::vector<ManifoldTrajectory<ScalarType>> ComputeInvariantManifolds(
    const PeriodicOrbit<ScalarType>& orbit, ScalarType mu, const ManifoldConfig<ScalarType>& config,
    ScalarType dt) {
  std::vector<ManifoldTrajectory<ScalarType>> result;

  // モノドロミー行列が計算済みか確認
  if (!orbit.stability_computed) {
    throw std::runtime_error(
        "ComputeInvariantManifolds: Stability analysis must be performed first");
  }

  // 固有値・固有ベクトルを計算
  auto decomp = ComputeEigenDecomposition(orbit.monodromy_matrix);

  // 安定固有値（|λ| < 1）と不安定固有値（|λ| > 1）のインデックスを特定
  std::vector<int> stable_indices;
  std::vector<int> unstable_indices;
  const ScalarType tolerance = 1e-6;

  for (int i = 0; i < 6; ++i) {
    ScalarType abs_eigenvalue = std::abs(decomp.eigenvalues[i]);
    if (abs_eigenvalue < 1.0 - tolerance) {
      stable_indices.push_back(i);
    } else if (abs_eigenvalue > 1.0 + tolerance) {
      unstable_indices.push_back(i);
    }
    // |λ| ≈ 1 の固有値は周期軌道自体の方向なので無視
  }

  if (stable_indices.empty() && unstable_indices.empty()) {
    std::cout << "[Info] No stable/unstable directions found (orbit may be stable)\n";
    return result;
  }

  // 周期軌道上の点を取得するために軌道を積分
  std::vector<State<ScalarType>> orbit_points;
  std::vector<std::array<std::array<ScalarType, 6>, 6>> orbit_stms;
  std::vector<ScalarType> orbit_times;

  ScalarType t = 0.0;
  ScalarType phase_step = orbit.period / config.num_initial_points;

  VariationalEquation<ScalarType> var_eq(mu);
  STMState<ScalarType> state;
  state.orbit = orbit.initial_state;

  using Stepper =
      boost::numeric::odeint::runge_kutta_dopri5<STMState<ScalarType>, ScalarType,
                                                 STMState<ScalarType>, ScalarType,
                                                 boost::numeric::odeint::vector_space_algebra>;
  auto controlled_stepper = boost::numeric::odeint::make_controlled(
      static_cast<ScalarType>(1e-12), static_cast<ScalarType>(1e-12), Stepper());

  std::cout << "    Computing orbit points..." << std::flush;
  for (int i = 0; i < config.num_initial_points; ++i) {
    ScalarType target_time = i * phase_step;
    boost::numeric::odeint::integrate_adaptive(controlled_stepper, var_eq, state, t, target_time,
                                               dt);
    t = target_time;

    State<ScalarType> point = state.orbit;
    orbit_points.push_back(point);
    orbit_stms.push_back(state.stm);
    orbit_times.push_back(t);
  }

  std::cout << " done (" << orbit_points.size() << " points)\n";

  // CR3BP equations of motion for manifold propagation
  auto cr3bp_eom = [mu](const std::array<ScalarType, 6>& s, std::array<ScalarType, 6>& dsdt,
                        ScalarType /* t */) {
    ScalarType x1 = s[0] + mu;
    ScalarType x2 = s[0] - (1.0 - mu);
    ScalarType r1_sq = x1 * x1 + s[1] * s[1] + s[2] * s[2];
    ScalarType r2_sq = x2 * x2 + s[1] * s[1] + s[2] * s[2];

    if (r1_sq < 1e-16 || r2_sq < 1e-16) {
      throw std::runtime_error("ComputeInvariantManifolds: Too close to primary");
    }

    ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
    ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

    dsdt[0] = s[3];                                                                    // dx/dt = vx
    dsdt[1] = s[4];                                                                    // dy/dt = vy
    dsdt[2] = s[5];                                                                    // dz/dt = vz
    dsdt[3] = 2.0 * s[4] + s[0] - (1.0 - mu) * x1 * r1_inv3 - mu * x2 * r2_inv3;       // dvx/dt
    dsdt[4] = -2.0 * s[3] + s[1] - (1.0 - mu) * s[1] * r1_inv3 - mu * s[1] * r2_inv3;  // dvy/dt
    dsdt[5] = -(1.0 - mu) * s[2] * r1_inv3 - mu * s[2] * r2_inv3;                      // dvz/dt
  };

  auto AddState = [](const std::array<ScalarType, 6>& s,
                     std::vector<State<ScalarType>>* states,
                     std::vector<ScalarType>* times, ScalarType t) {
    State<ScalarType> state = {s[0], s[1], s[2], s[3], s[4], s[5]};
    states->push_back(state);
    times->push_back(t);
  };

  auto IntegrateTrajectory = [&](std::array<ScalarType, 6>& manifold_state, ScalarType start_time,
                                 ScalarType end_time, ScalarType dt_step,
                                 std::vector<State<ScalarType>>* states,
                                 std::vector<ScalarType>* times) {
    const ScalarType dt_signed = (end_time >= start_time) ? std::abs(dt_step) : -std::abs(dt_step);

    if (config.integrator_type == "rk45") {
      using Stepper = boost::numeric::odeint::runge_kutta_dopri5<std::array<ScalarType, 6>>;
      auto controlled_stepper = boost::numeric::odeint::make_controlled(
          static_cast<ScalarType>(config.abs_tolerance),
          static_cast<ScalarType>(config.rel_tolerance), Stepper());

      auto observer = [&](const std::array<ScalarType, 6>& s, ScalarType t) {
        if (!times->empty()) {
          const ScalarType dt_last = std::abs(t - times->back());
          if (dt_last <= static_cast<ScalarType>(1e-15)) {
            return;
          }
        }
        AddState(s, states, times, t);
      };

      AddState(manifold_state, states, times, start_time);
      boost::numeric::odeint::integrate_adaptive(controlled_stepper, cr3bp_eom, manifold_state,
                                                 start_time, end_time, dt_signed, observer);
    } else {
      const ScalarType total_time = std::abs(end_time - start_time);
      const int num_steps = static_cast<int>(total_time / std::abs(dt_signed));
      ScalarType t = start_time;

      AddState(manifold_state, states, times, t);
      for (int step = 0; step < num_steps; ++step) {
        boost::numeric::odeint::integrate_const(
            boost::numeric::odeint::runge_kutta4<std::array<ScalarType, 6>>(), cr3bp_eom,
            manifold_state, t, t + dt_signed, dt_signed / 10.0);
        t += dt_signed;
        AddState(manifold_state, states, times, t);
      }
    }
  };

  // 安定多様体�計算（逆時間�分）
  if (config.compute_stable && !stable_indices.empty()) {
    int stable_idx = stable_indices[0];  // 最も強�安定方�
    auto& eigenvec = decomp.eigenvectors[stable_idx];

    // 固有ベクトル�実数部（t=0）
    std::array<ScalarType, 6> eigenvec0{};
    for (int i = 0; i < 6; ++i) {
      eigenvec0[i] = eigenvec[i].real();
    }

    std::cout << "    Computing stable manifolds: " << std::flush;
    for (size_t pt_idx = 0; pt_idx < orbit_points.size(); ++pt_idx) {
      if (pt_idx % 10 == 0) {
        std::cout << pt_idx << "/" << orbit_points.size() << "..." << std::flush;
      }
      const auto& base_point = orbit_points[pt_idx];
      const auto& stm = orbit_stms[pt_idx];

      // Transport eigenvector using STM at this phase
      std::array<ScalarType, 6> transported{};
      for (int i = 0; i < 6; ++i) {
        ScalarType sum = 0.0;
        for (int j = 0; j < 6; ++j) {
          sum += stm[i][j] * eigenvec0[j];
        }
        transported[i] = sum;
      }

      ScalarType norm = 0.0;
      for (int i = 0; i < 6; ++i) {
        norm += transported[i] * transported[i];
      }
      norm = std::sqrt(norm);

      // +方向と-方向の両方に摂動
      for (int dir = -1; dir <= 1; dir += 2) {
        ManifoldTrajectory<ScalarType> traj;
        traj.type = ManifoldTrajectory<ScalarType>::Type::STABLE;

        // 固有ベクトル方向に微小摂動を加える
        State<ScalarType> perturbed;
        perturbed.x = base_point.x + dir * config.epsilon * transported[0] / norm;
        perturbed.y = base_point.y + dir * config.epsilon * transported[1] / norm;
        perturbed.z = base_point.z + dir * config.epsilon * transported[2] / norm;
        perturbed.vx = base_point.vx + dir * config.epsilon * transported[3] / norm;
        perturbed.vy = base_point.vy + dir * config.epsilon * transported[4] / norm;
        perturbed.vz = base_point.vz + dir * config.epsilon * transported[5] / norm;

        ProjectToJacobiSurface(&perturbed, mu, orbit.jacobi_constant);

        traj.initial_displacement = perturbed;

        // 逆時間積分
        std::array<ScalarType, 6> manifold_state = {perturbed.x,  perturbed.y,  perturbed.z,
                                                    perturbed.vx, perturbed.vy, perturbed.vz};
        IntegrateTrajectory(manifold_state, static_cast<ScalarType>(0.0),
                           static_cast<ScalarType>(config.backward_time), dt,
                           &traj.trajectory, &traj.times);

        result.push_back(traj);
      }
    }
    std::cout << "done\n";
  }

  // 不安定多様体の計算（正時間積分）
  if (config.compute_unstable && !unstable_indices.empty()) {
    int unstable_idx = unstable_indices[0];  // 最も強い不安定方向
    auto& eigenvec = decomp.eigenvectors[unstable_idx];

    // Eigenvector at t=0
    std::array<ScalarType, 6> eigenvec0{};
    for (int i = 0; i < 6; ++i) {
      eigenvec0[i] = eigenvec[i].real();
    }

    std::cout << "    Computing unstable manifolds: " << std::flush;
    for (size_t pt_idx = 0; pt_idx < orbit_points.size(); ++pt_idx) {
      if (pt_idx % 10 == 0) {
        std::cout << pt_idx << "/" << orbit_points.size() << "..." << std::flush;
      }
      const auto& base_point = orbit_points[pt_idx];

      const auto& stm = orbit_stms[pt_idx];

      // Transport eigenvector using STM at this phase
      std::array<ScalarType, 6> transported{};
      for (int i = 0; i < 6; ++i) {
        ScalarType sum = 0.0;
        for (int j = 0; j < 6; ++j) {
          sum += stm[i][j] * eigenvec0[j];
        }
        transported[i] = sum;
      }

      ScalarType norm = 0.0;
      for (int i = 0; i < 6; ++i) {
        norm += transported[i] * transported[i];
      }
      norm = std::sqrt(norm);

      for (int dir = -1; dir <= 1; dir += 2) {
        ManifoldTrajectory<ScalarType> traj;
        traj.type = ManifoldTrajectory<ScalarType>::Type::UNSTABLE;

        State<ScalarType> perturbed;
        perturbed.x = base_point.x + dir * config.epsilon * transported[0] / norm;
        perturbed.y = base_point.y + dir * config.epsilon * transported[1] / norm;
        perturbed.z = base_point.z + dir * config.epsilon * transported[2] / norm;
        perturbed.vx = base_point.vx + dir * config.epsilon * transported[3] / norm;
        perturbed.vy = base_point.vy + dir * config.epsilon * transported[4] / norm;
        perturbed.vz = base_point.vz + dir * config.epsilon * transported[5] / norm;

        ProjectToJacobiSurface(&perturbed, mu, orbit.jacobi_constant);

        traj.initial_displacement = perturbed;

        // 正時間積分
        std::array<ScalarType, 6> manifold_state = {perturbed.x,  perturbed.y,  perturbed.z,
                                                    perturbed.vx, perturbed.vy, perturbed.vz};
        IntegrateTrajectory(manifold_state, static_cast<ScalarType>(0.0),
                           static_cast<ScalarType>(config.forward_time), dt,
                           &traj.trajectory, &traj.times);

        result.push_back(traj);
      }
    }
    std::cout << "done\n";
  }

  std::cout << "[Info] Computed " << result.size() << " manifold trajectories\n";
  return result;
}

}  // namespace periodic_orbit

#endif  // PERIODIC_ORBIT_MANIFOLD_HPP
