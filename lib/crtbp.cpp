#include "crtbp.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>

#include "vector3d.hpp"

//clang-format off
CRTBP::CRTBP(double init_x,
             double init_y,
             double init_z,
             double init_vx,
             double init_vy,
             double init_vz,
             double mu,
             double timestep)
    : init_posi_{init_x, init_y, init_z},
      posi_{init_x, init_y, init_z},
      init_v_{init_vx, init_vy, init_vz},
      v_{init_vx, init_vy, init_vz},
      init_q_{init_x, init_y, init_z},
      init_p_{init_vx - init_y, init_vy + init_x, init_vz},
      q_{init_x, init_y, init_z},
      mu_(3.003e-6),
      timestep_(timestep) {}

CRTBP::CRTBP(std::array<Vector3d, 2> init_state, double mu, double timestep)
    : init_posi_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      posi_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      init_v_{init_state[1].x(), init_state[1].y(), init_state[1].z()},
      v_{init_state[1].x(), init_state[1].y(), init_state[1].z()},
      init_q_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      init_p_{init_state[1].x() - init_state[0].y(),
              init_state[1].y() + init_state[0].x(),
              init_state[1].z()},
      q_{init_state[0].x(), init_state[0].y(), init_state[0].z()},
      p_{init_state[1].x() - init_state[0].y(),
         init_state[1].y() + init_state[0].x(),
         init_state[1].z()},
      mu_(mu),
      timestep_(timestep) {}

CRTBP::CRTBP(std::array<double, 6> init_state, double mu, double timestep)
    : init_posi_{init_state[0], init_state[1], init_state[2]},
      posi_{init_state[0], init_state[1], init_state[2]},
      init_v_{init_state[3], init_state[4], init_state[5]},
      v_{init_state[3], init_state[4], init_state[5]},
      init_q_{init_state[0], init_state[1], init_state[2]},
      init_p_{init_state[3] - init_state[0], init_state[4] + init_state[1], init_state[5]},
      q_{init_state[0], init_state[1], init_state[2]},
      p_{init_state[3] - init_state[0], init_state[4] + init_state[1], init_state[5]},
      mu_(3.003e-6),
      timestep_(timestep) {}
//clang-format on
double CRTBP::calc_Jacobi_integral() const {
    double r1 = calc_r1();
    double r2 = calc_r2();
    return posi_[0] * posi_[0] + posi_[1] * posi_[1] - (v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]) +
           2. * (1. - mu_) / r1 + 2. * mu_ / r2 + mu_ * (1. - mu_);
}

double CRTBP::calc_r1() const {
    return std::sqrt((posi_[0] + mu_) * (posi_[0] + mu_) + posi_[1] * posi_[1] +
                     posi_[2] * posi_[2]);
}

double CRTBP::calc_r2() const {
    return std::sqrt((posi_[0] - 1 + mu_) * (posi_[0] - 1 + mu_) + posi_[1] * posi_[1] +
                     posi_[2] * posi_[2]);
}

std::array<double, 6> CRTBP::init_state() const {
    std::array<double, 6> state;
    for (int i = 0; i < 3; i++) {
        state[i] = init_posi_[i];
    }
    state[3] = init_v_[0];
    state[4] = init_v_[1];
    state[5] = init_v_[2];
    return state;
}

std::array<double, 6> CRTBP::current_state() const {
    std::array<double, 6> state;
    for (int i = 0; i < 3; i++) {
        state[i] = posi_[i];
    }
    state[3] = v_[0];
    state[4] = v_[1];
    state[5] = v_[2];
    return state;
}

// 位置と速度をRK4で更新する
void CRTBP::RK4_step_noncanonical() {
    std::array<double, 6> k1, k2, k3, k4;
    std::array<double, 3> y1, y2, y3;

    // k1を計算
    k1[0] = v_[0];
    k1[1] = v_[1];
    k1[2] = v_[2];

    double r1_3 = std::pow(
        ((posi_[0] + mu_) * (posi_[0] + mu_) + posi_[1] * posi_[1] + posi_[2] * posi_[2]), 3. / 2.);
    double r2_3 = std::pow(
        ((posi_[0] - 1. + mu_) * (posi_[0] - 1. + mu_) + posi_[1] * posi_[1] + posi_[2] * posi_[2]),
        3. / 2.);
    //clang-format off
    k1[3] = 2*v_[1]+posi_[0]-(1-mu_)*(posi_[0]+mu_)/r1_3-mu_*(posi_[0]-1+mu_)/r2_3;
    k1[4] = -2*v_[0]+posi_[1]-(1-mu_)*posi_[1]/r1_3-mu_*posi_[1]/r2_3;
    k1[5] = -(1 - mu_) * posi_[2] / r1_3 - mu_ * posi_[2] / r2_3;
    //clang-format on
    // y1を計算
    for (int i = 0; i < 3; i++) {
        y1[i] = posi_[i] + k1[i] * timestep_ / 2;
    }

    // k2を計算
    k2[0] = v_[0] + k1[3] * timestep_ / 2;
    k2[1] = v_[1] + k1[4] * timestep_ / 2;
    k2[2] = v_[2] + k1[5] * timestep_ / 2;

    r1_3 = std::pow(((y1[0] + mu_) * (y1[0] + mu_) + y1[1] * y1[1] + y1[2] * y1[2]), 3. / 2);
    r2_3 =
        std::pow(((y1[0] - 1 + mu_) * (y1[0] - 1 + mu_) + y1[1] * y1[1] + y1[2] * y1[2]), 3. / 2);

    k2[3] = 2 * k2[1] + y1[0] - (1 - mu_) * (y1[0] + mu_) / r1_3 - mu_ * (y1[0] - 1 + mu_) / r2_3;
    k2[4] = -2 * k2[0] + y1[1] - (1 - mu_) * y1[1] / r1_3 - mu_ * y1[1] / r2_3;
    k2[5] = -(1 - mu_) * y1[2] / r1_3 - mu_ * y1[2] / r2_3;

    // y2を計算
    for (int i = 0; i < 3; i++) {
        y2[i] = posi_[i] + k2[i] * timestep_ / 2;
    }

    // k3を計算
    k3[0] = v_[0] + k2[3] * timestep_ / 2;
    k3[1] = v_[1] + k2[4] * timestep_ / 2;
    k3[2] = v_[2] + k2[5] * timestep_ / 2;

    r1_3 = std::pow(((y2[0] + mu_) * (y2[0] + mu_) + y2[1] * y2[1] + y2[2] * y2[2]), 3. / 2);
    r2_3 =
        std::pow(((y2[0] - 1 + mu_) * (y2[0] - 1 + mu_) + y2[1] * y2[1] + y2[2] * y2[2]), 3. / 2);

    k3[3] = 2 * k3[1] + y2[0] - (1 - mu_) * (y2[0] + mu_) / r1_3 - mu_ * (y2[0] - 1 + mu_) / r2_3;
    k3[4] = -2 * k3[0] + y2[1] - (1 - mu_) * y2[1] / r1_3 - mu_ * y2[1] / r2_3;
    k3[5] = -(1 - mu_) * y2[2] / r1_3 - mu_ * y2[2] / r2_3;

    // y3を計算
    for (int i = 0; i < 3; i++) {
        y3[i] = posi_[i] + k3[i] * timestep_;
    }

    // k4を計算
    k4[0] = v_[0] + k3[3] * timestep_;
    k4[1] = v_[1] + k3[4] * timestep_;
    k4[2] = v_[2] + k3[5] * timestep_;

    r1_3 = std::pow(((y3[0] + mu_) * (y3[0] + mu_) + y3[1] * y3[1] + y3[2] * y3[2]), 3. / 2);
    r2_3 =
        std::pow(((y3[0] - 1 + mu_) * (y3[0] - 1 + mu_) + y3[1] * y3[1] + y3[2] * y3[2]), 3. / 2);

    k4[3] = 2 * k4[1] + y3[0] - (1 - mu_) * (y3[0] + mu_) / r1_3 - mu_ * (y3[0] - 1 + mu_) / r2_3;
    k4[4] = -2 * k4[0] + y3[1] - (1 - mu_) * y3[1] / r1_3 - mu_ * y3[1] / r2_3;
    k4[5] = -(1 - mu_) * y3[2] / r1_3 - mu_ * y3[2] / r2_3;

    // 位置と速度を更新
    for (int i = 0; i < 3; i++) {
        posi_[i] += timestep_ / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        v_[i] += timestep_ / 6 * (k1[i + 3] + 2 * k2[i + 3] + 2 * k3[i + 3] + k4[i + 3]);
    }
}
void CRTBP::RK4_step_canonical() {
    std::array<double, 3> q;
    std::array<double, 3> p;
    std::array<double, 3> q_buf;
    std::array<double, 3> p_buf;

    auto calc_r1_ = [&](const std::array<double, 3>& q_temp) {
        return std::sqrt((q_temp[0] + mu_) * (q_temp[0] + mu_) + q_temp[1] * q_temp[1] +
                         q_temp[2] * q_temp[2]);
    };

    auto calc_r2_ = [&](const std::array<double, 3>& q_temp_) {
        return std::sqrt((q_temp_[0] - 1. + mu_) * (q_temp_[0] - 1. + mu_) +
                         q_temp_[1] * q_temp_[1] + q_temp_[2] * q_temp_[2]);
    };

    // copy for for loop
    for (int i = 0; i < 3; i++) {
        q[i] = posi_[i];
    }
    p[0] = v_[0] - posi_[1];
    p[1] = v_[1] + posi_[0];
    p[2] = v_[2];

    std::array<double, 3> k1, k2, k3, k4;
    std::array<double, 3> l1, l2, l3, l4;

    // k1を計算
    k1[0] = p[0] + q[1];
    k1[1] = p[1] - q[0];
    k1[2] = p[2];

    // l1を計算
    double r1_3 = std::pow(calc_r1_(q), 3);
    double r2_3 = std::pow(calc_r2_(q), 3);

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

std::vector<std::array<double, 8>> CRTBP::trajectory_propagation(double endtime) {

    double init_jacobi_integral = calc_Jacobi_integral();
    std::vector<std::array<double, 8>> trajectory;
    trajectory.reserve((endtime + 1) / timestep_);
    double time = 0.0;
    std::array<double, 8> state1;
    state1[0] = time;
    state1[1] = posi_[0];
    state1[2] = posi_[1];
    state1[3] = posi_[2];
    state1[4] = v_[0];
    state1[5] = v_[1]; 
    state1[6] = v_[2];
    state1[7] = calc_Jacobi_integral()-init_jacobi_integral;
    trajectory.push_back(state1);

    while (time < endtime) {

        RK4_step_noncanonical();

        time += timestep_;
        std::array<double, 8> state;
        state[0] = time;
        state[1] = posi_[0];
        state[2] = posi_[1];
        state[3] = posi_[2];
        state[4] = v_[0];
        state[5] = v_[1];
        state[6] = v_[2];
        state[7] = calc_Jacobi_integral()-init_jacobi_integral;
        
        trajectory.push_back(state);
    }
    return trajectory;
}
