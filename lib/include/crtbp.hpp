////////////////////////////////////////////////////////////////////////////////
/// @file           crtbp.hpp
/// @brief          円制限三体問題(CRTBP)を扱うクラス
/// @author         tabata
/// @date           2024/8/1
/// @par            edittor/ date/ version/ description
///                 tabata/ 2024/8/1/ 1.0/ 初版作成
///                 tabata/ 2024/8/1/ 2.0/ 初版作成
///                 tabata/ 2025/11/06 2.1 テンプレート化
///
////////////////////////////////////////////////////////////////////////////////

#ifndef CRTBP_H_
#define CRTBP_H_

#include <array>
#include <chrono>
#include <cmath>
#include <thread>
#include <vector>

#include "vector3d.hpp"

template <typename T>
class CRTBP {
 private:
  std::array<T, 3> init_posi_;
  std::array<T, 3> posi_;
  std::array<T, 3> init_v_;
  std::array<T, 3> v_;
  std::array<T, 3> init_q_;
  std::array<T, 3> init_p_;
  std::array<T, 3> q_;
  std::array<T, 3> p_;
  //! ratio of M_earth and M_earth+M_sun
  T mu_;
  //! time step for symplectic integration
  T timestep_;
  // coeficient c in the forth order symplectic method
  std::array<T, 4> c_ = {1. / (2. * (2. - std::pow(2., 1.0 / 3.0))),
                         (1. - std::pow(2., 1.0 / 3.0)) / (2. * (2. - std::pow(2., 1.0 / 3.0))),
                         (1. - std::pow(2., 1.0 / 3.0)) / (2. * (2. - std::pow(2., 1.0 / 3.0))),
                         1. / (2. * (2. - std::pow(2, 1.0 / 3.0)))};
  // coeficient d in the forth order symplectic method
  std::array<T, 4> d_ = {1. / (2. - std::pow(2., 1.0 / 3.0)),
                         -std::pow(2., 1.0 / 3.0) / (2. - std::pow(2., 1.0 / 3.0)),
                         1. / (2. - std::pow(2., 1.0 / 3.0)), 0.0};

  T calc_Jacobi_integral() const;

 public:
  /**
   * @brief Constructor
   * @param x Initial x-coordinate
   * @param y Initial y-coordinate
   * @param z Initial z-coordinate
   * @param vx Initial x-velocity
   * @param vy Initial y-velocity
   * @param vz Initial z-velocity
   */

  CRTBP(T init_x, T init_y, T init_z, T init_vx, T init_vy, T init_vz, T mu, T timestep = 0.001);

  CRTBP(std::array<Vector3d<T>, 2> init_state, T mu, T timestep = 0.001);
  CRTBP(std::array<T, 6> init_state, T mu, T timestep = 0.001);

  /// @brief 太陽とNEOの距離を計算する
  /// @return 太陽とNEOの距離
  T calc_r1() const;

  /// @brief 地球とNEOの距離を計算する
  /// @return 地球とNEOの距離
  T calc_r2() const;

  std::array<T, 6> init_state() const;
  std::array<T, 6> current_state() const;

  void RK4_step_noncanonical();
  void RK4_step_canonical();

  std::vector<std::array<T, 8>> trajectory_propagation(T endtime);
};

#include <array>
#include <cmath>
#include <iostream>
#include <memory>

#include "crtbp.hpp"
#include "vector3d.hpp"

//clang-format off
template <typename T>
CRTBP<T>::CRTBP(T init_x, T init_y, T init_z, T init_vx, T init_vy, T init_vz, T mu, T timestep)
    : init_posi_{init_x, init_y, init_z},
      posi_{init_x, init_y, init_z},
      init_v_{init_vx, init_vy, init_vz},
      v_{init_vx, init_vy, init_vz},
      init_q_{init_x, init_y, init_z},
      init_p_{init_vx - init_y, init_vy + init_x, init_vz},
      q_{init_x, init_y, init_z},
      mu_(mu),
      timestep_(timestep) {}

template <typename T>
CRTBP<T>::CRTBP(std::array<Vector3d<T>, 2> init_state, T mu, T timestep)
    : init_posi_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      posi_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      init_v_{init_state[1].x(), init_state[1].y(), init_state[1].z()},
      v_{init_state[1].x(), init_state[1].y(), init_state[1].z()},
      init_q_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      init_p_{init_state[1].x() - init_state[0].y(), init_state[1].y() + init_state[0].x(),
              init_state[1].z()},
      q_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      p_{init_state[1].x() - init_state[0].y(), init_state[1].y() + init_state[0].x(),
         init_state[1].z()},
      mu_(mu),
      timestep_(timestep) {}

template <typename T>
CRTBP<T>::CRTBP(std::array<T, 6> init_state, T mu, T timestep)
    : init_posi_{init_state[0], init_state[1], init_state[2]},
      posi_{init_state[0], init_state[1], init_state[2]},
      init_v_{init_state[3], init_state[4], init_state[5]},
      v_{init_state[3], init_state[4], init_state[5]},
      init_q_{init_state[0], init_state[1], init_state[2]},
      init_p_{init_state[3] - init_state[0], init_state[4] + init_state[1], init_state[5]},
      q_{init_state[0], init_state[1], init_state[2]},
      p_{init_state[3] - init_state[0], init_state[4] + init_state[1], init_state[5]},
      mu_(mu),
      timestep_(timestep) {}

template <typename T>
T CRTBP<T>::calc_Jacobi_integral() const {
  T r1 = calc_r1();
  T r2 = calc_r2();
  return posi_[0] * posi_[0] + posi_[1] * posi_[1] -
         (v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]) + 2. * (1. - mu_) / r1 + 2. * mu_ / r2 +
         mu_ * (1. - mu_);
}

template <typename T>
T CRTBP<T>::calc_r1() const {
  return std::sqrt((posi_[0] + mu_) * (posi_[0] + mu_) + posi_[1] * posi_[1] + posi_[2] * posi_[2]);
}

template <typename T>
T CRTBP<T>::calc_r2() const {
  return std::sqrt((posi_[0] - 1. + mu_) * (posi_[0] - 1. + mu_) + posi_[1] * posi_[1] +
                   posi_[2] * posi_[2]);
}

template <typename T>
std::array<T, 6> CRTBP<T>::init_state() const {
  std::array<T, 6> state;
  for (int i = 0; i < 3; i++) {
    state[i] = init_posi_[i];
  }
  state[3] = init_v_[0];
  state[4] = init_v_[1];
  state[5] = init_v_[2];
  return state;
}

template <typename T>
std::array<T, 6> CRTBP<T>::current_state() const {
  std::array<T, 6> state;
  for (int i = 0; i < 3; i++) {
    state[i] = posi_[i];
  }
  state[3] = v_[0];
  state[4] = v_[1];
  state[5] = v_[2];
  return state;
}

// 位置と速度をRK4で更新する
template <typename T>
void CRTBP<T>::RK4_step_noncanonical() {
  std::array<T, 6> k1, k2, k3, k4;
  std::array<T, 3> y1, y2, y3;
  // std::cout << std::setprecision(15) << "posi_: " << posi_[0] << ", " << posi_[1] << ", "
  //           << posi_[2] << ", " << v_[0] << ", " << v_[1] << ", " << v_[2] << std::endl;
  // k1を計算
  k1[0] = v_[0];
  k1[1] = v_[1];
  k1[2] = v_[2];

  T r1_3 = std::pow(
      ((posi_[0] + mu_) * (posi_[0] + mu_) + posi_[1] * posi_[1] + posi_[2] * posi_[2]), 3. / 2.);
  T r2_3 = std::pow(
      ((posi_[0] - 1. + mu_) * (posi_[0] - 1. + mu_) + posi_[1] * posi_[1] + posi_[2] * posi_[2]),
      3. / 2.);
  // (修正) r1_sq (2乗) を計算
  T r1_sq = (posi_[0] + mu_) * (posi_[0] + mu_) + posi_[1] * posi_[1] + posi_[2] * posi_[2];
  T r2_sq =
      (posi_[0] - 1. + mu_) * (posi_[0] - 1. + mu_) + posi_[1] * posi_[1] + posi_[2] * posi_[2];
  // (修正) r1_inv3 (逆3乗) を計算
  T r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
  T r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));
  //clang-format off
  k1[3] = 2 * v_[1] + posi_[0] - (1 - mu_) * (posi_[0] + mu_) / r1_3 -
          mu_ * (posi_[0] - 1 + mu_) / r2_3;
  k1[4] = -2 * v_[0] + posi_[1] - (1 - mu_) * posi_[1] / r1_3 - mu_ * posi_[1] / r2_3;
  k1[5] = -(1 - mu_) * posi_[2] / r1_3 - mu_ * posi_[2] / r2_3;
  // std::cout << std::setprecision(15)
  //           << "2 * vy, x, -(1-mu)*(x+mu)/r1^3, -mu*(x-1+mu)/r2^3: " << 2 * v_[1] << ", "
  //           << posi_[0] << ", " << -(mu_) << ", " << -mu_ * (posi_[0] - 1 + mu_) * r2_inv3
  //           << std::endl;
  // std::cout << std::setprecision(15) << "k1: " << k1[0] << ", " << k1[1] << ", " << k1[2] << ", "
  //           << k1[3] << ", " << k1[4] << ", " << k1[5] << std::endl;
  // y1を計算
  for (int i = 0; i < 3; i++) {
    y1[i] = posi_[i] + k1[i] * timestep_ / 2;
  }

  // k2を計算
  k2[0] = v_[0] + k1[3] * timestep_ / 2;
  k2[1] = v_[1] + k1[4] * timestep_ / 2;
  k2[2] = v_[2] + k1[5] * timestep_ / 2;

  r1_3 = std::pow(((y1[0] + mu_) * (y1[0] + mu_) + y1[1] * y1[1] + y1[2] * y1[2]), 3. / 2);
  r2_3 = std::pow(((y1[0] - 1 + mu_) * (y1[0] - 1 + mu_) + y1[1] * y1[1] + y1[2] * y1[2]), 3. / 2);

  k2[3] = 2 * k2[1] + y1[0] - (1 - mu_) * (y1[0] + mu_) / r1_3 - mu_ * (y1[0] - 1 + mu_) / r2_3;
  k2[4] = -2 * k2[0] + y1[1] - (1 - mu_) * y1[1] / r1_3 - mu_ * y1[1] / r2_3;
  k2[5] = -(1 - mu_) * y1[2] / r1_3 - mu_ * y1[2] / r2_3;

  // std::cout << std::setprecision(15) << "k2: " << k2[0] << ", " << k2[1] << ", " << k2[2] << ", "
  // << k2[3] << ", " << k2[4] << ", " << k2[5] << std::endl;

  // y2を計算
  for (int i = 0; i < 3; i++) {
    y2[i] = posi_[i] + k2[i] * timestep_ / 2;
  }

  // k3を計算
  k3[0] = v_[0] + k2[3] * timestep_ / 2;
  k3[1] = v_[1] + k2[4] * timestep_ / 2;
  k3[2] = v_[2] + k2[5] * timestep_ / 2;

  r1_3 = std::pow(((y2[0] + mu_) * (y2[0] + mu_) + y2[1] * y2[1] + y2[2] * y2[2]), 3. / 2);
  r2_3 = std::pow(((y2[0] - 1 + mu_) * (y2[0] - 1 + mu_) + y2[1] * y2[1] + y2[2] * y2[2]), 3. / 2);

  k3[3] = 2 * k3[1] + y2[0] - (1 - mu_) * (y2[0] + mu_) / r1_3 - mu_ * (y2[0] - 1 + mu_) / r2_3;
  k3[4] = -2 * k3[0] + y2[1] - (1 - mu_) * y2[1] / r1_3 - mu_ * y2[1] / r2_3;
  k3[5] = -(1 - mu_) * y2[2] / r1_3 - mu_ * y2[2] / r2_3;

  // std::cout << std::setprecision(15) << "k3: " << k3[0] << ", " << k3[1] << ", " << k3[2] << ", "
  // << k3[3] << ", " << k3[4] << ", " << k3[5] << std::endl;

  // y3を計算
  for (int i = 0; i < 3; i++) {
    y3[i] = posi_[i] + k3[i] * timestep_;
  }

  // k4を計算
  k4[0] = v_[0] + k3[3] * timestep_;
  k4[1] = v_[1] + k3[4] * timestep_;
  k4[2] = v_[2] + k3[5] * timestep_;

  r1_3 = std::pow(((y3[0] + mu_) * (y3[0] + mu_) + y3[1] * y3[1] + y3[2] * y3[2]), 3. / 2);
  r2_3 = std::pow(((y3[0] - 1 + mu_) * (y3[0] - 1 + mu_) + y3[1] * y3[1] + y3[2] * y3[2]), 3. / 2);

  k4[3] = 2 * k4[1] + y3[0] - (1 - mu_) * (y3[0] + mu_) / r1_3 - mu_ * (y3[0] - 1 + mu_) / r2_3;
  k4[4] = -2 * k4[0] + y3[1] - (1 - mu_) * y3[1] / r1_3 - mu_ * y3[1] / r2_3;
  k4[5] = -(1 - mu_) * y3[2] / r1_3 - mu_ * y3[2] / r2_3;

  // std::cout << std::setprecision(15) << "k4: " << k4[0] << ", " << k4[1] << ", " << k4[2] << ", "
  //           << k4[3] << ", " << k4[4] << ", " << k4[5] << std::endl;
  // 位置と速度を更新
  for (int i = 0; i < 3; i++) {
    posi_[i] += timestep_ / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    v_[i] += timestep_ / 6 * (k1[i + 3] + 2 * k2[i + 3] + 2 * k3[i + 3] + k4[i + 3]);
  }
}
template <typename T>
void CRTBP<T>::RK4_step_canonical() {
  std::array<T, 3> q;
  std::array<T, 3> p;
  std::array<T, 3> q_buf;
  std::array<T, 3> p_buf;

  auto calc_r1_ = [&](const std::array<T, 3>& q_temp) {
    return std::sqrt((q_temp[0] + mu_) * (q_temp[0] + mu_) + q_temp[1] * q_temp[1] +
                     q_temp[2] * q_temp[2]);
  };

  auto calc_r2_ = [&](const std::array<T, 3>& q_temp_) {
    return std::sqrt((q_temp_[0] - 1. + mu_) * (q_temp_[0] - 1. + mu_) + q_temp_[1] * q_temp_[1] +
                     q_temp_[2] * q_temp_[2]);
  };

  // copy for for loop
  for (int i = 0; i < 3; i++) {
    q[i] = posi_[i];
  }
  p[0] = v_[0] - posi_[1];
  p[1] = v_[1] + posi_[0];
  p[2] = v_[2];

  std::array<T, 3> k1, k2, k3, k4;
  std::array<T, 3> l1, l2, l3, l4;

  // k1を計算
  k1[0] = p[0] + q[1];
  k1[1] = p[1] - q[0];
  k1[2] = p[2];

  // l1を計算
  T r1_3 = std::pow(calc_r1_(q), 3);
  T r2_3 = std::pow(calc_r2_(q), 3);

  l1[0] = -(-p[1] - (1 - mu_) * q[0] / r1_3 - mu_ * q[0] / r2_3);
  l1[1] = -(p[0] - (1 - mu_) * q[1] / r1_3 - mu_ * q[1] / r2_3);
  l1[2] = -(-(1 - mu_) * q[2] / r1_3 - mu_ * q[2] / r2_3);

  for (int i = 0; i < 3; i++) {
    q_buf[i] = q[i] + k1[i] * timestep_ / 2;
    p_buf[i] = p[i] + l1[i] * timestep_ / 2;
  }

  // k2を計算
  k2[0] = p_buf[0] + q_buf[1];
  k2[1] = p_buf[1] - q_buf[0];
  k2[2] = p_buf[2];

  // l2を計算
  r1_3 = std::pow(calc_r1_(q_buf), 3);
  r2_3 = std::pow(calc_r2_(q_buf), 3);

  l2[0] = -(-p_buf[1] - (1 - mu_) * q_buf[0] / r1_3 - mu_ * q_buf[0] / r2_3);
  l2[1] = -(p_buf[0] - (1 - mu_) * q_buf[1] / r1_3 - mu_ * q_buf[1] / r2_3);
  l2[2] = -(-(1 - mu_) * q_buf[2] / r1_3 - mu_ * q_buf[2] / r2_3);

  for (int i = 0; i < 3; i++) {
    q_buf[i] = q[i] + k2[i] * timestep_ / 2;
    p_buf[i] = p[i] + l2[i] * timestep_ / 2;
  }

  // k3を計算
  k3[0] = p_buf[0] + q_buf[1];
  k3[1] = p_buf[1] - q_buf[0];
  k3[2] = p_buf[2];

  // l3を計算

  r1_3 = std::pow(calc_r1_(q_buf), 3);
  r2_3 = std::pow(calc_r2_(q_buf), 3);

  l3[0] = -(-p_buf[1] - (1 - mu_) * q_buf[0] / r1_3 - mu_ * q_buf[0] / r2_3);
  l3[1] = -(p_buf[0] - (1 - mu_) * q_buf[1] / r1_3 - mu_ * q_buf[1] / r2_3);
  l3[2] = -(-(1 - mu_) * q_buf[2] / r1_3 - mu_ * q_buf[2] / r2_3);

  for (int i = 0; i < 3; i++) {
    q_buf[i] = q[i] + k3[i] * timestep_;
    p_buf[i] = p[i] + l3[i] * timestep_;
  }

  // k4を計算
  k4[0] = p_buf[0] + q_buf[1];
  k4[1] = p_buf[1] - q_buf[0];
  k4[2] = p_buf[2];

  // l4を計算
  r1_3 = std::pow(calc_r1_(q_buf), 3);
  r2_3 = std::pow(calc_r2_(q_buf), 3);

  l4[0] = -(-p_buf[1] - (1 - mu_) * q_buf[0] / r1_3 - mu_ * q_buf[0] / r2_3);
  l4[1] = -(p_buf[0] - (1 - mu_) * q_buf[1] / r1_3 - mu_ * q_buf[1] / r2_3);
  l4[2] = -(-(1 - mu_) * q_buf[2] / r1_3 - mu_ * q_buf[2] / r2_3);

  // qとpを更新
  for (int i = 0; i < 3; i++) {
    q[i] += timestep_ / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    p[i] += timestep_ / 6 * (l1[i] + 2 * l2[i] + 2 * l3[i] + l4[i]);
  }

  // 位置と速度を更新
  for (int i = 0; i < 3; i++) {
    posi_[i] = q[i];
  }
  v_[0] = p[0] - q[1];
  v_[1] = p[1] + q[0];
  v_[2] = p[2];
}
template <typename T>
std::vector<std::array<T, 8>> CRTBP<T>::trajectory_propagation(T endtime) {
  T init_jacobi_integral = calc_Jacobi_integral();
  std::vector<std::array<T, 8>> trajectory;
  trajectory.reserve((endtime + 1) / timestep_);
  T time = 0.0;
  std::array<T, 8> state1;
  state1[0] = time;
  state1[1] = posi_[0];
  state1[2] = posi_[1];
  state1[3] = posi_[2];
  state1[4] = v_[0];
  state1[5] = v_[1];
  state1[6] = v_[2];
  state1[7] = calc_Jacobi_integral() - init_jacobi_integral;
  trajectory.push_back(state1);

  while (time < endtime) {
    RK4_step_noncanonical();

    time += timestep_;
    std::array<T, 8> state;
    state[0] = time;
    state[1] = posi_[0];
    state[2] = posi_[1];
    state[3] = posi_[2];
    state[4] = v_[0];
    state[5] = v_[1];
    state[6] = v_[2];
    state[7] = calc_Jacobi_integral() - init_jacobi_integral;

    trajectory.push_back(state);
  }
  return trajectory;
}

#endif
