/**
 * @file rtbp.hpp
 * @author tabata
 * @brief 制限三体問題の関数をまとめたライブラリ
 * @version 1.0
 * @date 2025-01-24
 * @note C++17
 * @note odeint
 * @note 楕円制限三体問題の関数が必要になったら名前空間と共に追加する
 * @par history
 * - 2025-01-24 tabata version 1.0
 */
#ifndef RTBP_HPP
#define RTBP_HPP

#include <array>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <iostream>
#include <utils.hpp>
// #include <vector3d.hpp>
#include <vector>

#include "utils.hpp"
#include "vector3d.hpp"

namespace rtbp {
inline static constexpr double G = 6.67430e-11;
inline static constexpr double AU = 1.49597870700e11;
inline static constexpr double M_SUN = 1.989400e30;
inline static constexpr double DEFAULT_MU = 3.003e-6;

// 3D座標を表す型のエイリアス
using Coord3D = std::array<double, 3>;
// 6D状態を表す型のエイリアス
using state = std::array<double, 6>;

/**
 * @brief Create a frame conversion matrix object
 * @param unit_n 規格化されたか回転軸
 * @param theta
 * @return std::array<std::array<double, 3>, 3>
 */
static std::array<std::array<double, 3>, 3>
create_frame_conversion_matrix(const Vector3d &unit_n, double theta) {
  double half_theta = theta / 2.0;
  double q0 = std::cos(half_theta);
  double sin_half_theta = std::sin(half_theta);
  double q1 = unit_n.x() * sin_half_theta;
  double q2 = unit_n.y() * sin_half_theta;
  double q3 = unit_n.z() * sin_half_theta;

  double q0q0 = q0 * q0;
  double q1q1 = q1 * q1;
  double q2q2 = q2 * q2;
  double q3q3 = q3 * q3;
  double q0q1 = q0 * q1;
  double q0q2 = q0 * q2;
  double q0q3 = q0 * q3;
  double q1q2 = q1 * q2;
  double q1q3 = q1 * q3;
  double q2q3 = q2 * q3;

  // 回転行列
  std::array<std::array<double, 3>, 3> convert_matrix = {
      {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 + q0q3), 2.0 * (q1q3 - q0q2)},
       {2.0 * (q1q2 - q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 + q0q1)},
       {2.0 * (q1q3 + q0q2), 2.0 * (q2q3 - q0q1), q0q0 - q1q1 - q2q2 + q3q3}}};
  return convert_matrix;
};

/**
 * @brief Create a coordinate rotation matrix object
 *
 * @param unit_n 規格化されたか回転軸
 * @param theta
 * @return std::array<std::array<double, 3>, 3>
 */
static std::array<std::array<double, 3>, 3>
create_coordinate_rotation_matrix(const Vector3d &unit_n, double theta) {
  double half_theta = theta / 2.0;
  double q0 = std::cos(half_theta);
  double sin_half_theta = std::sin(half_theta);
  double q1 = unit_n.x() * sin_half_theta;
  double q2 = unit_n.y() * sin_half_theta;
  double q3 = unit_n.z() * sin_half_theta;

  double q0q0 = q0 * q0;
  double q1q1 = q1 * q1;
  double q2q2 = q2 * q2;
  double q3q3 = q3 * q3;
  double q0q1 = q0 * q1;
  double q0q2 = q0 * q2;
  double q0q3 = q0 * q3;
  double q1q2 = q1 * q2;
  double q1q3 = q1 * q3;
  double q2q3 = q2 * q3;

  std::array<std::array<double, 3>, 3> rot_matrix = {
      {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 - q0q3), 2.0 * (q1q3 + q0q2)},
       {2.0 * (q1q2 + q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 - q0q1)},
       {2.0 * (q1q3 - q0q2), 2.0 * (q2q3 + q0q1), q0q0 - q1q1 - q2q2 + q3q3}}};
  return rot_matrix;
};

/**
 * @brief Create a Rodrigues rotation matrix object
 *
 * @param unit_n 規格化された回転軸
 * @param theta
 * @return std::array<std::array<double, 3>, 3>
 */
static std::array<std::array<double, 3>, 3>
create_rodrigues_matrix(const Vector3d &unit_n, double theta) {
  double c = std::cos(theta);
  double s = std::sin(theta);
  double t = 1.0 - c;
  double x = unit_n.x();
  double y = unit_n.y();
  double z = unit_n.z();

  std::array<std::array<double, 3>, 3> rot_matrix = {
      {{t * x * x + c, t * x * y - s * z, t * x * z + s * y},
       {t * x * y + s * z, t * y * y + c, t * y * z - s * x},
       {t * x * z - s * y, t * y * z + s * x, t * z * z + c}}};
  return rot_matrix;
};

/**
 * @brief Apply a rotation matrix to a vector *
 * @param rot_matrix matrix to apply
 * @param v vector to rotate
 * @return Vector3d
 */
static Vector3d apply_matrix(std::array<std::array<double, 3>, 3> rot_matrix,
                             const Vector3d &v) {
  return {rot_matrix[0][0] * v.x() + rot_matrix[0][1] * v.y() +
              rot_matrix[0][2] * v.z(),
          rot_matrix[1][0] * v.x() + rot_matrix[1][1] * v.y() +
              rot_matrix[1][2] * v.z(),
          rot_matrix[2][0] * v.x() + rot_matrix[2][1] * v.y() +
              rot_matrix[2][2] * v.z()};
};

/**
 * @brief 円制限三体問題の関数をまとめた名前空間
 *
 */
namespace crtbp {

double calc_SALI(const std::array<double, 6> &ref_state,
                 const std::array<double, 6> &perturbed_state1,
                 const std::array<double, 6> &perturbed_state2) {
  Vector3d deviation_vec1{perturbed_state1[0] - ref_state[0],
                          perturbed_state1[1] - ref_state[1],
                          perturbed_state1[2] - ref_state[2]};
  Vector3d deviation_vec2{perturbed_state2[0] - ref_state[0],
                          perturbed_state2[1] - ref_state[1],
                          perturbed_state2[2] - ref_state[2]};

  Vector3d normalized_dev_vec1 = deviation_vec1.normalise();
  Vector3d normalized_dev_vec2 = deviation_vec2.normalise();

  Vector3d sa = normalized_dev_vec1 - normalized_dev_vec2;
  Vector3d wa = normalized_dev_vec1 + normalized_dev_vec2;

  double sa_norm = sa.magnitude();
  double wa_norm = wa.magnitude();

  // SALIの計算
  return std::min(sa_norm, wa_norm);
}

// std::array<double ,3>かVector3dのみ受け付ける
template <typename T> double get_x(const T &r) {
  if constexpr (std::is_same_v<T, Vector3d>)
    return r.x();
  else if constexpr (std::is_same_v<T, Coord3D>)
    return r[0];
  else
    throw std::invalid_argument("Unsupported coordinate type");
}

// std::array<double ,3>かVector3dのみ受け付ける
template <typename T> double get_y(const T &r) {
  if constexpr (std::is_same_v<T, Vector3d>)
    return r.y();
  else if constexpr (std::is_same_v<T, Coord3D>)
    return r[1];
  else
    throw std::invalid_argument("Unsupported coordinate type");
}

// std::array<double ,3>かVector3dのみ受け付ける
template <typename T> double get_z(const T &r) {
  if constexpr (std::is_same_v<T, Vector3d>)
    return r.z();
  else if constexpr (std::is_same_v<T, Coord3D>)
    return r[2];
  else
    throw std::invalid_argument("Unsupported coordinate type");
}

// r1の計算,
template <typename T>
double calc_r1(const T &r, const AstroConstants &astro_params) {
  double x = get_x(r);
  double y = get_y(r);
  double z = get_z(r);
  return std::sqrt((x + astro_params.mu) * (x + astro_params.mu) + y * y +
                   z * z);
}

// r2の計算,
template <typename T> double calc_r2(const T &r, AstroConstants &astro_params) {
  double x = get_x(r);
  double y = get_y(r);
  double z = get_z(r);
  return std::sqrt((x - 1 + astro_params.mu) * (x - 1 + astro_params.mu) +
                   y * y + z * z);
}

// ポテンシャルエネルギーの計算
// template <typename T>
// double calc_potential_energy(const T &r, const AstroConstants &astro_params)
// {
//   double x = get_x(r);
//   double y = get_y(r);
//   double r1 = calc_r1(r, astro_params);
//   double r2 = calc_r2(r, astro_params);

//   if (r1 <= 0 || r2 <= 0) {
//     throw std::domain_error("Invalid distance calculated");
//   }

//   return 0.5 * (x * x + y * y) + (1 - astro_params.mu) / r1 +
//          astro_params.mu / r2;
// }

// template <typename T, std::size_t N>
// double calc_jacobi_integral(const std::array<T, N> &x,
//                             const AstroConstants &astro_params) {
//   double x_, y_, z_, vx_, vy_, vz_;
//   if constexpr (std::is_same_v<T, double> && N == 6) {
//     x_ = x[0];
//     y_ = x[1];
//     z_ = x[2];
//     vx_ = x[3];
//     vy_ = x[4];
//     vz_ = x[5];
//   } else if constexpr (std::is_same_v<T, Vector3d> && N == 2) {
//     x_ = x[0].x();
//     y_ = x[0].y();
//     z_ = x[0].z();
//     vx_ = x[1].x();
//     vy_ = x[1].y();
//     vz_ = x[1].z();
//   }

//   double r1 = std::sqrt((x_ + astro_params.mu) * (x_ + astro_params.mu) +
//                         y_ * y_ + z_ * z_);
//   double r2 =
//       std::sqrt((x_ - 1 + astro_params.mu) * (x_ - 1 + astro_params.mu) +
//                 y_ * y_ + z_ * z_);

//   return -(vx_ * vx_ + vy_ * vy_ + vz_ * vz_) +
//          2. * calc_potential_energy(x, astro_params);
// }

template <typename T>
double calc_v_abs(const T &r, const double JACOBI_INTEGRAL,
                  const AstroConstants &astro_params) {
  double r1 = calc_r1(r, astro_params.mu);
  double r2 = calc_r2(r, astro_params.mu);
  double x = get_x(r);
  double y = get_y(r);
  double z = get_z(r);
  return std::sqrt(x * x + y * y + 2. * (1. - astro_params.mu) / r1 +
                   2. * astro_params.mu / r2 +
                   astro_params.mu * (1. - astro_params.mu) - JACOBI_INTEGRAL);
}

template <typename T>
Vector3d calc_velocity(const T &point, const double v_abs,
                       const AstroConstants &astro_params,
                       const double inclination, const double OMEGA,
                       const double theta) {
  // 回転軸と回転角からクオータニオン経由で回転行列を生成
  auto create_rot_matrix =
      [](const Vector3d &unit_n,
         double theta_) -> std::array<std::array<double, 3>, 3> {
    double half_theta = theta_ / 2.0;
    double q0 = std::cos(half_theta);
    double sin_half_theta = std::sin(half_theta);
    double q1 = unit_n.x() * sin_half_theta;
    double q2 = unit_n.y() * sin_half_theta;
    double q3 = unit_n.z() * sin_half_theta;

    double q0q0 = q0 * q0;
    double q1q1 = q1 * q1;
    double q2q2 = q2 * q2;
    double q3q3 = q3 * q3;
    double q0q1 = q0 * q1;
    double q0q2 = q0 * q2;
    double q0q3 = q0 * q3;
    double q1q2 = q1 * q2;
    double q1q3 = q1 * q3;
    double q2q3 = q2 * q3;

    std::array<std::array<double, 3>, 3> rot_matrix = {
        {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 - q0q3), 2.0 * (q1q3 + q0q2)},
         {2.0 * (q1q2 + q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 - q0q1)},
         {2.0 * (q1q3 - q0q2), 2.0 * (q2q3 + q0q1),
          q0q0 - q1q1 - q2q2 + q3q3}}};
    return rot_matrix;
  };

  // 入力されたベクトルを入力された行列で変換
  auto convert = [](std::array<std::array<double, 3>, 3> convert_matrix,
                    const Vector3d &v) -> Vector3d {
    return {convert_matrix[0][0] * v.x() + convert_matrix[0][1] * v.y() +
                convert_matrix[0][2] * v.z(),
            convert_matrix[1][0] * v.x() + convert_matrix[1][1] * v.y() +
                convert_matrix[1][2] * v.z(),
            convert_matrix[2][0] * v.x() + convert_matrix[2][1] * v.y() +
                convert_matrix[2][2] * v.z()};
  };
  double vx_, vy_, vz_ = 0;

  // inclinationとOMEGAを用いて軌道面の法線ベクトルを計算
  Vector3d normal_vector{std::sin(inclination) * std::cos(OMEGA),
                         std::sin(inclination) * std::sin(OMEGA),
                         std::cos(inclination)};
  // 法線ベクトルと位置ベクトルの外積を計算
  double x, y, z;
  if constexpr (std::is_same_v<T, Vector3d>) {
    x = point.x();
    y = point.y();
    z = point.z();
  } else if constexpr (std::is_same_v<T, Coord3D>) {
    x = point[0];
    y = point[1];
    z = point[2];
  }
  Vector3d r2_vector{x - 1. + astro_params.mu, y, z};
  Vector3d h_vector = r2_vector.gaiseki(normal_vector);
  Vector3d normalized_h_vector = h_vector.normalise();

  // 速度ベクトルを計算
  vx_ = v_abs * normalized_h_vector.x();
  vy_ = v_abs * normalized_h_vector.y();
  vz_ = v_abs * normalized_h_vector.z();

  if (theta == 0.0)
    return Vector3d(vx_, vy_, vz_);
  else {
    // クオータニオンを用いて速度ベクトルをnormal_vector周りにthetaだけ回転させる
    std::array<std::array<double, 3>, 3> rot_matrix =
        create_rot_matrix(normal_vector, theta);
    Vector3d velocity{vx_, vy_, vz_};
    Vector3d rotated_velocity = convert(rot_matrix, velocity);

    return rotated_velocity;
  }
}

template <typename T, std::size_t N>
const std::array<T, N>
ConvertInertial2Rotating(const std::array<T, N> &ast_state,
                         const std::array<T, N> &p2_state,
                         const AstroConstants &astro_params) {
  // std::array<double, 6>かstd::array<Vector3d, 2>のみ受け付ける
  static_assert((std::is_same_v<T, double> && N == 6) ||
                    (std::is_same_v<T, Vector3d> && N == 2),
                "Array must be either std::array<double, 6> or "
                "std::array<Vector3d, 2>");
  std::array<Vector3d, 2> init_ast_state_G;
  std::array<Vector3d, 2> init_p2_state_G;
  if constexpr (std::is_same_v<T, double> && N == 6) {
    init_ast_state_G = {Vector3d(ast_state[0], ast_state[1], ast_state[2]),
                        Vector3d(ast_state[3], ast_state[4], ast_state[5])};
    init_p2_state_G = {Vector3d(p2_state[0], p2_state[1], p2_state[2]),
                       Vector3d(p2_state[3], p2_state[4], p2_state[5])};
  } else if constexpr (std::is_same_v<T, Vector3d> && N == 2) {
    init_ast_state_G = {ast_state[0], ast_state[1]};
    init_p2_state_G = {p2_state[0], p2_state[1]};
  }

  // 慣性座標系での基底ベクトル
  Vector3d G_x{1.0, 0.0, 0.0};
  Vector3d G_y{0.0, 1.0, 0.0};
  Vector3d G_z{0.0, 0.0, 1.0};

  /*x軸を一致させる中間座標系rot1と移す  */
  // 慣性系から見た地球ベクトルと慣性系のx軸がなす角度を計算
  double theta1 = std::atan2(init_p2_state_G[0].gaiseki(G_x).magnitude(),
                             init_p2_state_G[0].naiseki(G_x));
  std::cout << std::setprecision(15);

  // 回転軸を決定
  Vector3d rotax1 = G_x.gaiseki(init_p2_state_G[0]).normalise();
  // 変換行列を生成
  std::array<std::array<double, 3>, 3> convert_G_to_R1 =
      create_rodrigues_matrix(rotax1, -theta1);
  std::array<std::array<double, 3>, 3> rotate1 =
      create_rodrigues_matrix(rotax1, theta1);

  /* 地球の軌道面の法線ベクトルと中間座標系rot1のz軸を一致させ目的の回転座標系rot2へ移す
   */
  // 地球の軌道面の法線ベクトル
  Vector3d n_plane = init_p2_state_G[0].gaiseki(init_p2_state_G[1]).normalise();

  // 慣性系から見た中間座標系rot1のz軸
  Vector3d rot1_z = apply_matrix(rotate1, G_z);

  //    慣性座標系から見た地球軌道面の法線ベクトルと慣性系から見た中間座標系rot1のz軸がなす角度を計算
  double theta2 =
      std::atan2(n_plane.gaiseki(rot1_z).magnitude(), n_plane.naiseki(rot1_z));
  Vector3d rotax2 = {1.0, 0.0, 0.0};

  // 変換行列を生成
  std::array<std::array<double, 3>, 3> convert_R1_to_R2 =
      create_rodrigues_matrix(rotax2, -theta2);
  std::array<std::array<double, 3>, 3> rotate2 =
      create_rodrigues_matrix(rotax2, theta2);
  // 小惑星の位置を回転座標系に変換
  Vector3d ast_pos_R1 = apply_matrix(convert_G_to_R1, init_ast_state_G[0]);
  Vector3d epos_R1 = apply_matrix(convert_G_to_R1, init_p2_state_G[0]);
  Vector3d ast_pos_R2 = apply_matrix(convert_R1_to_R2, ast_pos_R1);
  Vector3d epos_R2 = apply_matrix(convert_R1_to_R2, epos_R1);
  Vector3d ast_pos_R{ast_pos_R2.x() + 1.0 - epos_R2.x() - astro_params.mu,
                     ast_pos_R2.y(), ast_pos_R2.z()};

  /* calc the velocity of rotating frame */
  //  変換行列convert_G_to_R1の微分を計算
  // theta1の微分

  double ND_time_ref = std::sqrt(astro_params.au * astro_params.au /
                                 astro_params.gm_sun); // Non-Dimensional time
  double mean_motion = 1 / ND_time_ref;
  std::cout << "mean motion = " << mean_motion << std::endl;
  // conbert AU/day to Non-Dimensional
  Vector3d ND_init_e_velo_G = n_plane.gaiseki(init_p2_state_G[0]).normalise() *
                              init_p2_state_G[1].magnitude() * ND_time_ref /
                              86400.;

  double e_G_x = init_p2_state_G[0].x();
  double e_G_y = init_p2_state_G[0].y();
  double e_G_z = init_p2_state_G[0].z();
  double ND_e_G_vx = ND_init_e_velo_G.x();
  double ND_e_G_vy = ND_init_e_velo_G.y();
  double ND_e_G_vz = ND_init_e_velo_G.z();

  double poo = init_p2_state_G[0].magnitude();
  double sin_theta1 = std::sqrt(1. - std::cos(theta1) * std::cos(theta1));

  // double dot_theta1 = -(ND_e_G_vx * poo * poo -
  //                       e_G_x * (e_G_x * ND_e_G_vx + e_G_y * ND_e_G_vy +
  //                       e_G_z * ND_e_G_vz)) /
  //                       <-
  //                       円制限三体問題だからちきゅうの速度と位置の内積は0
  //                     poo / poo / poo / sin_theta1;
  double dot_theta1 = -ND_e_G_vx / poo / sin_theta1;

  // theta2の微分
  double sin_theta2 = std::sqrt(1. - std::cos(theta2) * std::cos(theta2));
  if (n_plane.gaiseki(apply_matrix(rotate1, G_x)).naiseki(rot1_z) < 0) {
    sin_theta2 = -sin_theta2;
  }

  Vector3d k{n_plane.x(), n_plane.y(), n_plane.z()};

  double n1 = std::sqrt(e_G_y * e_G_y + e_G_z * e_G_z);
  double a = e_G_y * ND_e_G_vy + e_G_z * ND_e_G_vz;
  double b = ND_e_G_vy * e_G_z + e_G_y * ND_e_G_vz;

  double dot_theta2 =
      (k.x() * (ND_e_G_vz * n1 * n1 - e_G_z * a) * sin_theta1 / n1 / n1 / n1 -
       k.x() * e_G_z * std ::cos(theta1) * dot_theta1 / n1 +
       k.y() * (b * n1 * n1 - 2. * e_G_z * e_G_y * a) *
           (1. - std::cos(theta1)) / n1 / n1 / n1 / n1 +
       k.y() * e_G_z * ND_e_G_vz * sin_theta1 * dot_theta1 / n1 / n1 -
       2. * k.z() * (e_G_y * ND_e_G_vy * n1 * n1 - e_G_y * a) *
           (1. - std::cos(theta1)) / n1 / n1 / n1 / n1 -
       k.z() * e_G_y * e_G_y * sin_theta1 * dot_theta1 / n1 / n1 +
       k.z() * sin_theta1 * dot_theta1) /
      sin_theta2;

  double c1 = G_x.naiseki(init_p2_state_G[0]) / init_p2_state_G[0].magnitude();
  double s1 = sin_theta1;
  double c2 = std::cos(theta2);
  double s2 = sin_theta2;
  std::array<std::array<double, 3>, 3> dot_convert_G_to_R1{
      {{-s1 * dot_theta1,
        (ND_e_G_vy * n1 * n1 - e_G_y * a) * s1 / n1 / n1 / n1 +
            e_G_y * c1 * dot_theta1 / n1,
        (ND_e_G_vz * n1 * n1 - e_G_z * a) * s1 / n1 / n1 / n1 +
            e_G_z * c1 * dot_theta1 / n1},
       {-(ND_e_G_vy * n1 * n1 - e_G_y * a) * s1 / n1 / n1 / n1 -
            e_G_y * c1 * dot_theta1 / n1,
        (2. * e_G_z * ND_e_G_vz * n1 * n1 - 2. * e_G_z * e_G_z * a) *
                (1.0 - c1) / n1 / n1 / n1 / n1 +
            e_G_z * e_G_z * s1 * dot_theta1 / n1 / n1 - s1 * dot_theta1,
        -(b * n1 * n1 - 2. * e_G_z * e_G_y * a) * (1. - c1) / n1 / n1 / n1 /
                n1 -
            e_G_z * e_G_y * s1 * dot_theta1 / n1 / n1},
       {-(ND_e_G_vz * n1 * n1 - e_G_z * a) * s1 / n1 / n1 / n1 -
            e_G_z * c1 * dot_theta1 / n1,
        -(b * n1 * n1 - 2. * e_G_z * e_G_y * a) * (1.0 - c1) / n1 / n1 / n1 /
                n1 -
            e_G_z * e_G_y * s1 * dot_theta1 / n1 / n1,
        (2. * e_G_y * ND_e_G_vy * n1 * n1 - e_G_y * e_G_y * a) * (1. - c1) /
                n1 / n1 / n1 / n1 +
            e_G_y * e_G_y * s1 * dot_theta1 / n1 / n1 - s1 * dot_theta1}}};

  std::array<std::array<double, 3>, 3> dot_convert_R1_to_R2{
      {{0.0, 0.0, 0.0},
       {0.0, -s2 * dot_theta2, c2 * dot_theta2},
       {0.0, -c2 * dot_theta2, -s2 * dot_theta2}}};
  // 速度を座標変換
  Vector3d ast_vel_R_1 = apply_matrix(
      dot_convert_R1_to_R2, apply_matrix(convert_G_to_R1, init_ast_state_G[0]));
  Vector3d ast_vel_R_2 = apply_matrix(
      convert_R1_to_R2, apply_matrix(dot_convert_G_to_R1, init_ast_state_G[0]));
  Vector3d ast_vel_R_3 = apply_matrix(
      convert_R1_to_R2, apply_matrix(convert_G_to_R1, init_ast_state_G[1]));

  Vector3d ast_vel_R = ast_vel_R_1 + ast_vel_R_2 + ast_vel_R_3;
  // check
  Vector3d e_vel_R_1 = {0, 0, 0};
  // apply_matrix(dot_convert_R1_to_R2, apply_matrix(convert_G_to_R1,
  // init_p2_state_G[0]));
  Vector3d e_vel_R_2 = {0, 0, 0};
  // apply_matrix(convert_R1_to_R2, apply_matrix(dot_convert_G_to_R1,
  // init_p2_state_G[0]));
  Vector3d e_vel_R_3 = apply_matrix(
      convert_R1_to_R2, apply_matrix(convert_G_to_R1, ND_init_e_velo_G));
  Vector3d e_vel_R = e_vel_R_1 + e_vel_R_2 + e_vel_R_3;

#ifdef _debug

  // 中間座標系にから見た地球の速度
  // OK!!!!!

#endif

  // check
  std::cout << "converted earth velocity = " << e_vel_R.x() << " "
            << e_vel_R.y() << " " << e_vel_R.z() << std::endl;
  std::cout << "<>               Asteroid data converted to rotating frame"
            << std::endl;
  if constexpr (std::is_same_v<T, double> && N == 6) {
    return {ast_pos_R.x(), ast_pos_R.y(), ast_pos_R.z(),
            ast_vel_R.x(), ast_vel_R.y(), ast_vel_R.z()};
  } else if constexpr (std::is_same_v<T, Vector3d> && N == 2) {
    return {ast_pos_R, ast_vel_R};
  }
}

// CRTBPの運動方程式
// EquationOfMotionを継承
// class EquationOfMotion : public rtbp::EquationOfMotion {
// private:
//   AstroConstants astro_params_;

// public:
//   explicit EquationOfMotion(const AstroConstants &astro_params)
//       : astro_params_(astro_params) {}

//   void operator()(const state &x, state &dxdt, const double /* t */) const {
//     const double x_ = x[0], y_ = x[1], z_ = x[2];
//     const double vx_ = x[3], vy_ = x[4], vz_ = x[5];

//     const double r1 = calc_r1(x, astro_params_);
//     const double r2 = calc_r2(x, astro_params_);

//     dxdt[0] = vx_;
//     dxdt[1] = vy_;
//     dxdt[2] = vz_;
//     dxdt[3] =
//         2. * vy_ + x_ -
//         (1. - astro_params_.mu) * (x_ + astro_params_.mu) / (r1 * r1 * r1) -
//         astro_params_.mu * (x_ - 1. + astro_params_.mu) / (r2 * r2 * r2);
//     dxdt[4] = -2. * vx_ + y_ - (1. - astro_params_.mu) * y_ / (r1 * r1 * r1)
//     -
//               astro_params_.mu * y_ / (r2 * r2 * r2);
//     dxdt[5] = -(1. - astro_params_.mu) * z_ / (r1 * r1 * r1) -
//               astro_params_.mu * z_ / (r2 * r2 * r2);
//   }
// };
// class observer {
// private:
//   double init_jacobi_;
//   AstroConstants astro_params_;

// public:
//   std::vector<std::array<double, 8>> &m_out;
//   explicit observer(double init_jacobi, const AstroConstants &astro_params,
//                     std::vector<std::array<double, 8>> &out)
//       : init_jacobi_(init_jacobi), astro_params_(astro_params), m_out(out) {}

//   void operator()(const state &x, double t) const {
//     double jacobi_integral = calc_jacobi_integral(x, astro_params_);
//     m_out.push_back({t, x[0], x[1], x[2], x[3], x[4], x[5],
//                      init_jacobi_ - jacobi_integral});
//   }
// };

// class Trajectory {
// private:
//   double jacobi_integral_;
//   AstroConstants astro_params_;

// public:
//   std::vector<std::array<double, 8>> trajectory_;

//   Trajectory(double jacobi_integral, const AstroConstants &astro_params)
//       : jacobi_integral_(jacobi_integral), astro_params_(astro_params) {}

//   std::vector<std::array<double, 8>> integrate(const state &x0, double t0,
//                                                double tf, double dt) {
//     state init_state = x0;
//     trajectory_.clear();
//     trajectory_.reserve(static_cast<size_t>((tf - t0) / dt) + 1);

//     auto eom = EquationOfMotion(astro_params_);
//     auto obs = observer(jacobi_integral_, astro_params_, trajectory_);
//     boost::numeric::odeint::controlled_runge_kutta<
//         boost::numeric::odeint::runge_kutta_fehlberg78<state>>
//         stepper;

//     boost::numeric::odeint::integrate_const(stepper, eom, init_state, t0, tf,
//                                             dt, obs);

//     return trajectory_;
//   }
// };

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>

// (ユーザー提供の Vector3d クラス定義がここにあると仮定)
// (ユーザー提供の AstroConstants 構造体定義がここにあると仮定)
// AstroConstants には .au, .gm_sun, .gm_earth が [m] と [m^3/s^2] で
// 正しく設定されている必要がある。

namespace {

// 3x3 Matrix type alias
using Matrix3x3 = std::array<std::array<double, 3>, 3>;

// State vector type aliases
using StateVector6 = std::array<double, 6>;
using StateVectorPair = std::array<Vector3d, 2>;

// 物理定数
constexpr double kSecondsPerDay = 86400.0;

/**
 * @brief 3x3回転行列を3Dベクトルに適用します。
 */
Vector3d ApplyMatrix(const Matrix3x3 &matrix, const Vector3d &vector) {
  return Vector3d(matrix[0][0] * vector.x() + matrix[0][1] * vector.y() +
                      matrix[0][2] * vector.z(),
                  matrix[1][0] * vector.x() + matrix[1][1] * vector.y() +
                      matrix[1][2] * vector.z(),
                  matrix[2][0] * vector.x() + matrix[2][1] * vector.y() +
                      matrix[2][2] * vector.z());
}

/**
 * @brief 6要素の状態配列 [x, y, z, vx, vy, vz] を
 * [pos, vel] の Vector3d ペアに変換します。
 */
StateVectorPair ConvertStateToPair(const StateVector6 &state) {
  return {Vector3d(state[0], state[1], state[2]),
          Vector3d(state[3], state[4], state[5])};
}

/**
 * @brief [pos, vel] の Vector3d ペアを 6要素の状態配列
 * [x, y, z, vx, vy, vz] に変換します。
 */
StateVector6 ConvertPairToState(const StateVectorPair &state_pair) {
  return {state_pair[0].x(), state_pair[0].y(), state_pair[0].z(),
          state_pair[1].x(), state_pair[1].y(), state_pair[1].z()};
}

} // namespace

/**
 * @brief J2000慣性系（太陽中心, AU, AU/day）の状態ベクトルを、
 * 太陽-P2回転座標系（共通重心中心, 無次元）に変換します。
 *
 * この関数は、P2（例：地球）の一般的な（傾斜した）軌道を扱います。
 * 回転座標系 (R) は以下のように定義されます：
 * - 原点: 太陽-P2 共通重心
 * - x軸: 太陽からP2を指す方向
 * - z軸: P2の軌道角運動量ベクトル (h = r_p2 x v_p2) に平行
 * - y軸: 右手系を完成させる (y = z x x)
 *
 * 出力は、長さ単位(LU)=r_p2、時間単位(TU)=1/n で無次元化されます。
 *
 * @tparam T std::arrayの要素型 (double または Vector3d)。
 * @tparam N std::arrayのサイズ (6 または 2)。
 * @param ast_state J2000慣性系（太陽中心）における小惑星の状態ベクトル
 * [pos(AU), vel(AU/day)]。
 * @param p2_state J2000慣性系（太陽中心）におけるP2（地球）の状態ベクトル
 * [pos(AU), vel(AU/day)]。
 * @param astro_params 天文定数（au(m), gm_sun(m^3/s^2), gm_earth(m^3/s^2)）。
 * @return 回転座標系（無次元）における小惑星の状態ベクトル [pos_nd, vel_nd]。
 */
template <typename T, std::size_t N>
const std::array<T, N> FrameTransformation(const std::array<T, N> &ast_state,
                                           const std::array<T, N> &p2_state,
                                           const AstroConstants &astro_params) {
  // 正しいテンプレート引数か静的に検証
  static_assert(
      (std::is_same_v<T, double> && N == 6) ||
          (std::is_same_v<T, Vector3d> && N == 2),
      "Input must be std::array<double, 6> or std::array<Vector3d, 2>");

  // 0. 単位換算 (m^3/s^2 -> AU^3/day^2)
  if (astro_params.au == 0) {
    throw std::runtime_error("AstroConstants::au cannot be zero.");
  }
  const double kM3s2_to_AU3day2 =
      (kSecondsPerDay * kSecondsPerDay) /
      (astro_params.au * astro_params.au * astro_params.au);
  const double gm_sun_AD = astro_params.gm_sun * kM3s2_to_AU3day2;
  const double gm_earth_AD = astro_params.gm_earth * kM3s2_to_AU3day2;
  const double gm_total_AD = gm_sun_AD + gm_earth_AD;

  // 1. 入力を {pos, vel} の Vector3d ペアに統一
  //    (G = J2000 Inertial Frame, 単位: AU, AU/day)
  StateVectorPair ast_state_G;
  StateVectorPair p2_state_G;
  if constexpr (std::is_same_v<T, double> && N == 6) {
    ast_state_G = ConvertStateToPair(ast_state);
    p2_state_G = ConvertStateToPair(p2_state);
  } else if constexpr (std::is_same_v<T, Vector3d> && N == 2) {
    ast_state_G = {ast_state[0], ast_state[1]};
    p2_state_G = {p2_state[0], p2_state[1]};
  }

  const Vector3d r_ast_G = ast_state_G[0];
  const Vector3d v_ast_G = ast_state_G[1];
  const Vector3d r_p2_G = p2_state_G[0];
  const Vector3d v_p2_G = p2_state_G[1];

  // 2. 質量比と共通重心 (Barycenter, B) の状態を計算
  // 質量比は単位系に依存しない
  const double mu =
      astro_params.gm_earth / (astro_params.gm_sun + astro_params.gm_earth);

  // 太陽（G系の原点）から見た共通重心の状態 (AU, AU/day)
  const Vector3d r_B_G = r_p2_G * mu;
  const Vector3d v_B_G = v_p2_G * mu;

  // 3. G系における回転座標系 (R) の基底ベクトルを定義
  const double r_p2_mag = r_p2_G.magnitude(); // これが LU (AU)
  if (r_p2_mag == 0) {
    throw std::runtime_error("P2 position magnitude is zero.");
  }
  const Vector3d i_R_G = r_p2_G / r_p2_mag; // x軸 (太陽 -> P2)

  const Vector3d h_G = r_p2_G.gaiseki(v_p2_G); // P2の軌道角運動量 (AU^2/day)
  const double h_mag = h_G.magnitude();
  if (h_mag == 0) {
    throw std::runtime_error(
        "P2 angular momentum magnitude is zero (e.g., linear motion).");
  }
  const Vector3d k_R_G = h_G / h_mag; // z軸 (h の方向)

  const Vector3d j_R_G = k_R_G.gaiseki(i_R_G); // y軸 (右手系を完成)

  // 4. 回転行列 C (G系 -> R系) を作成
  const Matrix3x3 C_G2R = {{{i_R_G.x(), i_R_G.y(), i_R_G.z()},
                            {j_R_G.x(), j_R_G.y(), j_R_G.z()},
                            {k_R_G.x(), k_R_G.y(), k_R_G.z()}}};

  // 5. 位置ベクトルの変換 (単位: AU)
  const Vector3d r_ast_rel_B_G = r_ast_G - r_B_G;
  const Vector3d ast_pos_R = ApplyMatrix(C_G2R, r_ast_rel_B_G); // (AU)

  // 6. 回転系の角速度ベクトル (omega_R) を計算 (単位: rad/day)
  // P2の加速度 (AU/day^2)
  const Vector3d a_p2_G =
      r_p2_G * (-gm_sun_AD / (r_p2_mag * r_p2_mag * r_p2_mag));

  const double dot_r_p2_mag = r_p2_G.naiseki(v_p2_G) / r_p2_mag; // (AU/day)
  const Vector3d dot_i_R_G =
      (v_p2_G / r_p2_mag) -
      (r_p2_G * (dot_r_p2_mag / (r_p2_mag * r_p2_mag))); // (1/day)

  const Vector3d dot_h_G = r_p2_G.gaiseki(a_p2_G);       // (AU^2/day^2)
  const double dot_h_mag = h_G.naiseki(dot_h_G) / h_mag; // (AU^2/day^2)
  const Vector3d dot_k_R_G =
      (dot_h_G / h_mag) - (h_G * (dot_h_mag / (h_mag * h_mag))); // (1/day)

  const Vector3d dot_j_R_G =
      dot_k_R_G.gaiseki(i_R_G) + k_R_G.gaiseki(dot_i_R_G); // (1/day)

  const double omega_x_R = dot_j_R_G.naiseki(k_R_G); // (rad/day)
  const double omega_y_R = dot_k_R_G.naiseki(i_R_G); // (rad/day)
  const double omega_z_R = dot_i_R_G.naiseki(j_R_G); // (rad/day)

  const Vector3d omega_R(omega_x_R, omega_y_R, omega_z_R); // (rad/day)

  // 7. 速度ベクトルの変換 (単位: AU/day)
  const Vector3d v_ast_rel_B_G = v_ast_G - v_B_G;                   // (AU/day)
  const Vector3d v_R_transport = ApplyMatrix(C_G2R, v_ast_rel_B_G); // (AU/day)
  const Vector3d v_R_coriolis = omega_R.gaiseki(ast_pos_R); // (rad/day * AU)
  const Vector3d ast_vel_R = v_R_transport - v_R_coriolis;  // (AU/day)

  // 8. 無次元化
  // 長さの単位 (LU)
  const double lu = r_p2_mag; // (AU)

  // 時間の単位 (TU)
  // n = sqrt( (GM_sun + GM_earth) / LU^3 )  (rad/day)
  const double mean_motion_AD = std::sqrt(gm_total_AD / (lu * lu * lu));
  if (mean_motion_AD == 0) {
    throw std::runtime_error("Mean motion is zero, cannot non-dimensionalize.");
  }
  const double tu = 1.0 / mean_motion_AD; // (day)

  // 速度の単位 (VU)
  const double vu = lu / tu; // (AU/day)

  // 無次元化の実行
  const Vector3d ast_pos_R_nd = ast_pos_R / lu;
  const Vector3d ast_vel_R_nd = ast_vel_R / vu;

  // 9. 出力のフォーマットを整える
  if constexpr (std::is_same_v<T, double> && N == 6) {
    return ConvertPairToState({ast_pos_R_nd, ast_vel_R_nd});
  } else if constexpr (std::is_same_v<T, Vector3d> && N == 2) {
    return {ast_pos_R_nd, ast_vel_R_nd};
  }
}

} // namespace crtbp
} // namespace rtbp

#endif // RTBP_HPP