////////////////////////////////////////////////////////////////////////////////
/// @file           crtbp.hpp
/// @brief          円制限三体問題(CRTBP)を扱うクラス
/// @author         tabata
/// @date           2024/8/1
/// @par            edittor/ date/ version/ description
///                 tabata/ 2024/8/1/ 1.0/ 初版作成
///                 tabata/ 2024/8/1/ 2.0/ 初版作成
///
////////////////////////////////////////////////////////////////////////////////

#ifndef CRTBP_H_
#define CRTBP_H_

#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "vector3d.hpp"

class CRTBP {
   private:
    std::array<double, 3> init_posi_;
    std::array<double, 3> posi_;
    std::array<double, 3> init_v_;
    std::array<double, 3> v_;
    std::array<double, 3> init_q_;
    std::array<double, 3> init_p_;
    std::array<double, 3> q_;
    std::array<double, 3> p_;
    //! ratio of M_earth and M_earth+M_sun
    double mu_;
    //! time step for symplectic integration
    double timestep_;
    // coeficient c in the forth order symplectic method
    std::array<double, 4> c_ = {
        1. / (2. * (2. - std::pow(2., 1.0 / 3.0))),
        (1. - std::pow(2., 1.0 / 3.0)) / (2. * (2. - std::pow(2., 1.0 / 3.0))),
        (1. - std::pow(2., 1.0 / 3.0)) / (2. * (2. - std::pow(2., 1.0 / 3.0))),
        1. / (2. * (2. - std::pow(2, 1.0 / 3.0)))};
    // coeficient d in the forth order symplectic method
    std::array<double, 4> d_ = {1. / (2. - std::pow(2., 1.0 / 3.0)),
                                -std::pow(2., 1.0 / 3.0) / (2. - std::pow(2., 1.0 / 3.0)),
                                1. / (2. - std::pow(2., 1.0 / 3.0)),
                                0.0};

    double calc_Jacobi_integral() const;

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

    CRTBP(double init_x,
          double init_y,
          double init_z,
          double init_vx,
          double init_vy,
          double init_vz,
          double mu,
          double timestep = 0.001);

    CRTBP(std::array<Vector3d, 2> init_state, double mu, double timestep = 0.001);
    CRTBP(std::array<double, 6> init_state, double mu, double timestep = 0.001);

    /// @brief 太陽とNEOの距離を計算する
    /// @return 太陽とNEOの距離
    double calc_r1() const;

    /// @brief 地球とNEOの距離を計算する
    /// @return 地球とNEOの距離
    double calc_r2() const;

    std::array<double, 6> init_state() const;
    std::array<double, 6> current_state() const;

    void RK4_step_noncanonical();
    void RK4_step_canonical();

    std::vector<std::array<double, 8>> trajectory_propagation(double endtime);
};

#endif