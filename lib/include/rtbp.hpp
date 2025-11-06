/**
 * @file rtbp.hpp
 * @author tabata
 * @brief 制限三体問題の関数をまとめたライブラリ
 * @version 1.1
 * @date 2025-11-06
 * @note C++17
 * @note odeint
 * @note 楕円制限三体問題の関数が必要になったら名前空間と共に追加する
 * @par history
 * - 2025-01-24 tabata version 1.0
 * - 2025-11-06 refactored template parameter names
 */
#ifndef RTBP_HPP
#define RTBP_HPP

#include <array>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <iostream>
#include <vector3d.hpp>

namespace my_type {
// 3D座標を表す型のエイリアス
template <typename ScalarType>
using Coord3D = Vector3d<ScalarType>;
// 6D状態を表す型のエイリアス
template <typename ScalarType>
using StateVector6 = std::array<ScalarType, 6>;
// 3x3 Matrix type alias
template <typename ScalarType>
using Matrix3x3 = std::array<std::array<ScalarType, 3>, 3>;

template <typename ScalarType>
using StateVectorPair = std::array<Vector3d<ScalarType>, 2>;
/**
 * @brief　普通の座標を保存
 */
template <typename ScalarType>
struct State {
  Vector3d<ScalarType> r;
  Vector3d<ScalarType> v;

  ScalarType& x() { return r.x; }
  const ScalarType& x() const { return r.x; }

  ScalarType& y() { return r.y; }
  const ScalarType& y() const { return r.y; }

  // z の読み書き
  ScalarType& z() { return r.z; }
  const ScalarType& z() const { return r.z; }

  ScalarType& vx() { return v.x; }
  const ScalarType& vx() const { return v.x; }

  ScalarType& vy() { return v.y; }
  const ScalarType& vy() const { return v.y; }

  ScalarType& vz() { return v.z; }
  const ScalarType& vz() const { return v.z; }
};
};  // namespace my_type
/**
 * @brief スタティックな関数
 */
namespace {
using namespace my_type;
// 物理定数
constexpr double kSecondsPerDay = 86400.0;  // P -> double

template <typename ScalarType, typename CoordType>
ScalarType get_x(const CoordType& r) {
  if constexpr (std::is_same_v<CoordType, Vector3d<ScalarType>>)
    return r.x();
  else if constexpr (std::is_same_v<CoordType, Coord3D<ScalarType>>)
    return r[0];
  else if constexpr (std::is_same_v<CoordType, StateVector6<ScalarType>>)
    return r[0];
  else
    throw std::invalid_argument("Unsupported coordinate type for get_x");
}

template <typename ScalarType, typename CoordType>
ScalarType get_y(const CoordType& r) {
  if constexpr (std::is_same_v<CoordType, Vector3d<ScalarType>>)
    return r.y();
  else if constexpr (std::is_same_v<CoordType, Coord3D<ScalarType>>)
    return r[1];
  else if constexpr (std::is_same_v<CoordType, StateVector6<ScalarType>>)
    return r[1];
  else
    throw std::invalid_argument("Unsupported coordinate type for get_y");
}

template <typename ScalarType, typename CoordType>
ScalarType get_z(const CoordType& r) {
  if constexpr (std::is_same_v<CoordType, Vector3d<ScalarType>>)
    return r.z();
  else if constexpr (std::is_same_v<CoordType, Coord3D<ScalarType>>)
    return r[2];
  else if constexpr (std::is_same_v<CoordType, StateVector6<ScalarType>>)
    return r[2];
  else
    throw std::invalid_argument("Unsupported coordinate type for get_z");
}

/**
 * @brief Create a frame conversion matrix object
 * @param unit_n 規格化された回転軸
 * @param theta 回転角
 * @return Matrix3x3<ScalarType>
 */
template <typename ScalarType>
constexpr Matrix3x3<ScalarType> create_frame_conversion_matrix(const Vector3d<ScalarType>& unit_n,
                                                               ScalarType theta) {
  ScalarType half_theta = theta / 2.0;
  ScalarType q0 = std::cos(half_theta);
  ScalarType sin_half_theta = std::sin(half_theta);
  ScalarType q1 = unit_n.x() * sin_half_theta;
  ScalarType q2 = unit_n.y() * sin_half_theta;
  ScalarType q3 = unit_n.z() * sin_half_theta;

  ScalarType q0q0 = q0 * q0;
  ScalarType q1q1 = q1 * q1;
  ScalarType q2q2 = q2 * q2;
  ScalarType q3q3 = q3 * q3;
  ScalarType q0q1 = q0 * q1;
  ScalarType q0q2 = q0 * q2;
  ScalarType q0q3 = q0 * q3;
  ScalarType q1q2 = q1 * q2;
  ScalarType q1q3 = q1 * q3;
  ScalarType q2q3 = q2 * q3;

  // 回転行列
  std::array<std::array<ScalarType, 3>, 3> convert_matrix = {
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
template <typename ScalarType>
constexpr Matrix3x3<ScalarType> create_coordinate_rotation_matrix(
    const Vector3d<ScalarType>& unit_n, ScalarType theta) {
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

  Matrix3x3<ScalarType> rot_matrix = {
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
template <typename ScalarType>
constexpr Matrix3x3<ScalarType> create_rodrigues_matrix(const Vector3d<ScalarType>& unit_n,
                                                        ScalarType theta) {
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
template <typename ScalarType>
Vector3d<ScalarType> ApplyMatrix(const Matrix3x3<ScalarType> rot_matrix,
                                 const Vector3d<ScalarType>& v) {
  return {rot_matrix[0][0] * v.x() + rot_matrix[0][1] * v.y() + rot_matrix[0][2] * v.z(),
          rot_matrix[1][0] * v.x() + rot_matrix[1][1] * v.y() + rot_matrix[1][2] * v.z(),
          rot_matrix[2][0] * v.x() + rot_matrix[2][1] * v.y() + rot_matrix[2][2] * v.z()};
};

/**
 * @brief 6要素の状態配列 [x, y, z, vx, vy, vz] を
 * [pos, vel] の Vector3d ペアに変換
 */
template <typename ScalarType>
StateVectorPair<ScalarType> ConvertStateToPair(const StateVector6<ScalarType>& state) {
  return {Vector3d(state[0], state[1], state[2]), Vector3d(state[3], state[4], state[5])};
}

/**
 * @brief [pos, vel] の Vector3d ペアを 6要素の状態配列
 * [x, y, z, vx, vy, vz] に変換
 */
template <typename ScalarType>
StateVector6<ScalarType> ConvertPairToState(const StateVectorPair<ScalarType>& state_pair) {
  return {state_pair[0].x(), state_pair[0].y(), state_pair[0].z(),
          state_pair[1].x(), state_pair[1].y(), state_pair[1].z()};
}
}  // namespace

namespace param {
/**
 * @brief 構造体の概要説明
 *
 * 詳細説明
 */
template <typename ScalarType>
struct AstroConstants {
  ScalarType au;        ///< 天文単位 (m)
  ScalarType gm_sun;    ///< 太陽の重力定数 (m^3/s^2)
  ScalarType gm_earth;  ///< 地球の重力定数 (m^3/s^2)
  ScalarType mu;        ///< 質量比
  ScalarType G;         ///< 万有引力定数 (kg^-1 m^3 s^-2)
};
}  // namespace param

/**
 * @brief 円制限三体問題の関数をまとめた名前空間
 *
 */
namespace crtbp {
using namespace param;
using namespace my_type;
template <typename ElementType, typename ScalarType, std::size_t Dim>
const ScalarType CalcSALI3d(const std::array<ElementType, Dim>& ref_state,
                            const std::array<ElementType, Dim>& perturbed_state1,
                            const std::array<ElementType, Dim>& perturbed_state2) {
  static_assert((Dim == 6) || (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2),
                "Array must be either std::array<ScalarType, 6> or "
                "std::array<Vector3d<ScalarType>, 2>");
  Vector3d<ScalarType> deviation_vec1;
  Vector3d<ScalarType> deviation_vec2;
  if constexpr (Dim == 6) {
    deviation_vec1 = {perturbed_state1[0] - ref_state[0], perturbed_state1[1] - ref_state[1],
                      perturbed_state1[2] - ref_state[2]};
    deviation_vec2 = {perturbed_state2[0] - ref_state[0], perturbed_state2[1] - ref_state[1],
                      perturbed_state2[2] - ref_state[2]};
  } else if constexpr (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2) {
    deviation_vec1 = {perturbed_state1[0].x() - ref_state[0].x(),
                      perturbed_state1[0].y() - ref_state[0].y(),
                      perturbed_state1[0].z() - ref_state[0].z()};
    deviation_vec2 = {perturbed_state2[0].x() - ref_state[0].x(),
                      perturbed_state2[0].y() - ref_state[0].y(),
                      perturbed_state2[0].z() - ref_state[0].z()};
  }

  Vector3d<ScalarType> normalized_dev_vec1 = deviation_vec1.normalise();
  Vector3d<ScalarType> normalized_dev_vec2 = deviation_vec2.normalise();

  Vector3d<ScalarType> sa = normalized_dev_vec1 - normalized_dev_vec2;
  Vector3d<ScalarType> wa = normalized_dev_vec1 + normalized_dev_vec2;

  ScalarType sa_norm = sa.magnitude();
  ScalarType wa_norm = wa.magnitude();

  // SALIの計算
  return std::min(sa_norm, wa_norm);
}

template <typename ElementType, typename ScalarType, std::size_t Dim>
const ScalarType CalcSALI6d(const std::array<ElementType, Dim>& ref_state,
                            const std::array<ElementType, Dim>& perturbed_state1,
                            const std::array<ElementType, Dim>& perturbed_state2) {
  static_assert((Dim == 6) || (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2),
                "Array must be either std::array<ScalarType, 6> or "
                "std::array<Vector3d<ScalarType>, 2>");
  std::array<ScalarType, 6> deviation_vec1;
  std::array<ScalarType, 6> deviation_vec2;
  if constexpr (Dim == 6) {
    deviation_vec1 = {perturbed_state1[0] - ref_state[0], perturbed_state1[1] - ref_state[1],
                      perturbed_state1[2] - ref_state[2], perturbed_state1[3] - ref_state[3],
                      perturbed_state1[4] - ref_state[4], perturbed_state1[5] - ref_state[5]};
    deviation_vec2 = {perturbed_state2[0] - ref_state[0], perturbed_state2[1] - ref_state[1],
                      perturbed_state2[2] - ref_state[2], perturbed_state2[3] - ref_state[3],
                      perturbed_state2[4] - ref_state[4], perturbed_state2[5] - ref_state[5]};
  } else if constexpr (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2) {
    deviation_vec1 = {
        perturbed_state1[0].x() - ref_state[0].x(), perturbed_state1[0].y() - ref_state[0].y(),
        perturbed_state1[0].z() - ref_state[0].z(), perturbed_state1[1].x() - ref_state[1].x(),
        perturbed_state1[1].y() - ref_state[1].y(), perturbed_state1[1].z() - ref_state[1].z()};
    deviation_vec2 = {
        perturbed_state2[0].x() - ref_state[0].x(), perturbed_state2[0].y() - ref_state[0].y(),
        perturbed_state2[0].z() - ref_state[0].z(), perturbed_state2[1].x() - ref_state[1].x(),
        perturbed_state2[1].y() - ref_state[1].y(), perturbed_state2[1].z() - ref_state[1].z()};
  }
  ScalarType perturbed1_norm =
      std::sqrt(deviation_vec1[0] * deviation_vec1[0] + deviation_vec1[1] * deviation_vec1[1] +
                deviation_vec1[2] * deviation_vec1[2] + deviation_vec1[3] * deviation_vec1[3] +
                deviation_vec1[4] * deviation_vec1[4] + deviation_vec1[5] * deviation_vec1[5]);
  ScalarType perturbed2_norm =
      std::sqrt(deviation_vec2[0] * deviation_vec2[0] + deviation_vec2[1] * deviation_vec2[1] +
                deviation_vec2[2] * deviation_vec2[2] + deviation_vec2[3] * deviation_vec2[3] +
                deviation_vec2[4] * deviation_vec2[4] + deviation_vec2[5] * deviation_vec2[5]);
  std::array<ScalarType, 6> normalized_dev_vec1 = {
      deviation_vec1[0] / perturbed1_norm, deviation_vec1[1] / perturbed1_norm,
      deviation_vec1[2] / perturbed1_norm, deviation_vec1[3] / perturbed1_norm,
      deviation_vec1[4] / perturbed1_norm, deviation_vec1[5] / perturbed1_norm};
  std::array<ScalarType, 6> normalized_dev_vec2 = {
      deviation_vec2[0] / perturbed2_norm, deviation_vec2[1] / perturbed2_norm,
      deviation_vec2[2] / perturbed2_norm, deviation_vec2[3] / perturbed2_norm,
      deviation_vec2[4] / perturbed2_norm, deviation_vec2[5] / perturbed2_norm};

  std::array<ScalarType, 6> sa_vector = {normalized_dev_vec1[0] - normalized_dev_vec2[0],
                                         normalized_dev_vec1[1] - normalized_dev_vec2[1],
                                         normalized_dev_vec1[2] - normalized_dev_vec2[2],
                                         normalized_dev_vec1[3] - normalized_dev_vec2[3],
                                         normalized_dev_vec1[4] - normalized_dev_vec2[4],
                                         normalized_dev_vec1[5] - normalized_dev_vec2[5]};
  std::array<ScalarType, 6> wa_vector = {normalized_dev_vec1[0] + normalized_dev_vec2[0],
                                         normalized_dev_vec1[1] + normalized_dev_vec2[1],
                                         normalized_dev_vec1[2] + normalized_dev_vec2[2],
                                         normalized_dev_vec1[3] + normalized_dev_vec2[3],
                                         normalized_dev_vec1[4] + normalized_dev_vec2[4],
                                         normalized_dev_vec1[5] + normalized_dev_vec2[5]};
  ScalarType sa_norm = std::sqrt(sa_vector[0] * sa_vector[0] + sa_vector[1] * sa_vector[1] +
                                 sa_vector[2] * sa_vector[2] + sa_vector[3] * sa_vector[3] +
                                 sa_vector[4] * sa_vector[4] + sa_vector[5] * sa_vector[5]);
  ScalarType wa_norm = std::sqrt(wa_vector[0] * wa_vector[0] + wa_vector[1] * wa_vector[1] +
                                 wa_vector[2] * wa_vector[2] + wa_vector[3] * wa_vector[3] +
                                 wa_vector[4] * wa_vector[4] + wa_vector[5] * wa_vector[5]);
  // SALIの計算
  ScalarType SALI = std::min(sa_norm, wa_norm);
  return SALI;
}

// r1の計算,
template <typename ScalarType>
ScalarType calc_r1(const Coord3D<ScalarType>& r, const AstroConstants<ScalarType>& astro_params) {
  ScalarType x = r.x();
  ScalarType y = r.y();
  ScalarType z = r.z();
  return std::sqrt((x + astro_params.mu) * (x + astro_params.mu) + y * y + z * z);
}

// r2の計算,
template <typename ScalarType>
ScalarType calc_r2(const Coord3D<ScalarType>& r, const AstroConstants<ScalarType>& astro_params) {
  ScalarType x = r.x();
  ScalarType y = r.y();
  ScalarType z = r.z();
  return std::sqrt((x - 1 + astro_params.mu) * (x - 1 + astro_params.mu) + y * y + z * z);
}

// ポテンシャルエネルギーの計算
template <typename ScalarType>
ScalarType calc_potential_energy(const Coord3D<ScalarType>& r,
                                 const AstroConstants<ScalarType>& astro_params) {
  ScalarType x = get_x<ScalarType>(r);
  ScalarType y = get_y<ScalarType>(r);
  ScalarType r1 = calc_r1<ScalarType>(r, astro_params);
  ScalarType r2 = calc_r2<ScalarType>(r, astro_params);

  if (r1 <= 0 || r2 <= 0) {
    throw std::domain_error("Invalid distance calculated");
  }

  return 0.5 * (x * x + y * y) + (1 - astro_params.mu) / r1 + astro_params.mu / r2;
}

template <typename ScalarType>
ScalarType calc_jacobi_integral(const State<ScalarType>& state,
                                const AstroConstants<ScalarType>& astro_params) {
  ScalarType x_, y_, z_, vx_, vy_, vz_;
  x_ = state.x();
  y_ = state.y();
  z_ = state.z();
  vx_ = state.vx();
  vy_ = state.vy();
  vz_ = state.z();

  const State<ScalarType> pos = {{x_, y_, z_}, {0, 0, 0}};

  ScalarType r1 = calc_r1<ScalarType>(pos, astro_params);
  ScalarType r2 = calc_r2<ScalarType>(pos, astro_params);

  return -(vx_ * vx_ + vy_ * vy_ + vz_ * vz_) +
         2. * calc_potential_energy<ScalarType>(pos, astro_params);
}

template <typename ScalarType>
ScalarType calc_v_abs(const Coord3D<ScalarType>& r, const ScalarType JACOBI_INTEGRAL,
                      const AstroConstants<ScalarType>& astro_params) {
  ScalarType r1 = calc_r1<ScalarType>(r, astro_params.mu);
  ScalarType r2 = calc_r2<ScalarType>(r, astro_params.mu);
  ScalarType x = r.x();
  ScalarType y = r.y();
  ScalarType z = r.z();
  return std::sqrt(x * x + y * y + 2. * (1. - astro_params.mu) / r1 + 2. * astro_params.mu / r2 +
                   astro_params.mu * (1. - astro_params.mu) - JACOBI_INTEGRAL);
}

template <typename ScalarType>
Vector3d<ScalarType> calc_velocity(const Coord3D<ScalarType>& point, const ScalarType v_abs,
                                   const AstroConstants<ScalarType>& astro_params,
                                   const ScalarType inclination, const ScalarType OMEGA,
                                   const ScalarType theta) {
  // 回転軸と回転角からクオータニオン経由で回転行列を生成
  auto create_rot_matrix = [](const Vector3d<ScalarType>& unit_n,
                              ScalarType theta_) -> std::array<std::array<ScalarType, 3>, 3> {
    ScalarType half_theta = theta_ / 2.0;
    ScalarType q0 = std::cos(half_theta);
    ScalarType sin_half_theta = std::sin(half_theta);
    ScalarType q1 = unit_n.x() * sin_half_theta;
    ScalarType q2 = unit_n.y() * sin_half_theta;
    ScalarType q3 = unit_n.z() * sin_half_theta;

    ScalarType q0q0 = q0 * q0;
    ScalarType q1q1 = q1 * q1;
    ScalarType q2q2 = q2 * q2;
    ScalarType q3q3 = q3 * q3;
    ScalarType q0q1 = q0 * q1;
    ScalarType q0q2 = q0 * q2;
    ScalarType q0q3 = q0 * q3;
    ScalarType q1q2 = q1 * q2;
    ScalarType q1q3 = q1 * q3;
    ScalarType q2q3 = q2 * q3;

    std::array<std::array<ScalarType, 3>, 3> rot_matrix = {
        {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 - q0q3), 2.0 * (q1q3 + q0q2)},
         {2.0 * (q1q2 + q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 - q0q1)},
         {2.0 * (q1q3 - q0q2), 2.0 * (q2q3 + q0q1), q0q0 - q1q1 - q2q2 + q3q3}}};
    return rot_matrix;
  };

  // 入力されたベクトルを入力された行列で変換
  auto convert = [](std::array<std::array<ScalarType, 3>, 3> convert_matrix,
                    const Vector3d<ScalarType>& v) -> Vector3d<ScalarType> {
    return {
        convert_matrix[0][0] * v.x() + convert_matrix[0][1] * v.y() + convert_matrix[0][2] * v.z(),
        convert_matrix[1][0] * v.x() + convert_matrix[1][1] * v.y() + convert_matrix[1][2] * v.z(),
        convert_matrix[2][0] * v.x() + convert_matrix[2][1] * v.y() + convert_matrix[2][2] * v.z()};
  };
  ScalarType vx_ = 0.0, vy_ = 0.0, vz_ = 0.0;
  // inclinationとOMEGAを用いて軌道面の法線ベクトルを計算
  Vector3d<ScalarType> normal_vector{std::sin(inclination) * std::cos(OMEGA),
                                     std::sin(inclination) * std::sin(OMEGA),
                                     std::cos(inclination)};
  // 法線ベクトルと位置ベクトルの外積を計算
  ScalarType x, y, z;

  x = point.x();
  y = point.y();
  z = point.z();

  Vector3d<ScalarType> r2_vector{x - 1. + astro_params.mu, y, z};
  Vector3d<ScalarType> h_vector = r2_vector.gaiseki(normal_vector);
  Vector3d<ScalarType> normalized_h_vector = h_vector.normalise();

  // 速度ベクトルを計算
  vx_ = v_abs * normalized_h_vector.x();
  vy_ = v_abs * normalized_h_vector.y();
  vz_ = v_abs * normalized_h_vector.z();

  if (theta == 0.0)
    return Vector3d<ScalarType>(vx_, vy_, vz_);
  else {
    // クオータニオンを用いて速度ベクトルをnormal_vector周りにthetaだけ回転させる
    std::array<std::array<ScalarType, 3>, 3> rot_matrix = create_rot_matrix(normal_vector, theta);
    Vector3d<ScalarType> velocity{vx_, vy_, vz_};
    Vector3d<ScalarType> rotated_velocity = convert(rot_matrix, velocity);

    return rotated_velocity;
  }
}

template <typename ElementType, typename ScalarType, std::size_t Dim>
const std::array<ElementType, Dim> ConvertInertial2Rotating(
    const std::array<ElementType, Dim>& ast_state, const std::array<ElementType, Dim>& p2_state,
    const AstroConstants<ScalarType>& astro_params) {
  // std::array<ScalarType, 6>かstd::array<Vector3d<ScalarType>, 2>のみ受け付ける
  static_assert((std::is_same_v<ElementType, ScalarType> && Dim == 6) ||
                    (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2),
                "Array must be either std::array<ScalarType, 6> or "
                "std::array<Vector3d<ScalarType>, 2>");
  std::array<Vector3d<ScalarType>, 2> init_ast_state_G;
  std::array<Vector3d<ScalarType>, 2> init_p2_state_G;
  if constexpr (std::is_same_v<ElementType, ScalarType> && Dim == 6) {
    init_ast_state_G = {Vector3d<ScalarType>(ast_state[0], ast_state[1], ast_state[2]),
                        Vector3d<ScalarType>(ast_state[3], ast_state[4], ast_state[5])};
    init_p2_state_G = {Vector3d<ScalarType>(p2_state[0], p2_state[1], p2_state[2]),
                       Vector3d<ScalarType>(p2_state[3], p2_state[4], p2_state[5])};
  } else if constexpr (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2) {
    init_ast_state_G = {ast_state[0], ast_state[1]};
    init_p2_state_G = {p2_state[0], p2_state[1]};
  }

  // 慣性座標系での基底ベクトル
  Vector3d<ScalarType> G_x{1.0, 0.0, 0.0};
  Vector3d<ScalarType> G_y{0.0, 1.0, 0.0};
  Vector3d<ScalarType> G_z{0.0, 0.0, 1.0};

  /*x軸を一致させる中間座標系rot1と移す  */
  // 慣性系から見た地球ベクトルと慣性系のx軸がなす角度を計算
  ScalarType theta1 =
      std::atan2(init_p2_state_G[0].gaiseki(G_x).magnitude(), init_p2_state_G[0].naiseki(G_x));
  std::cout << std::setprecision(15);

  // 回転軸を決定
  Vector3d<ScalarType> rotax1 = G_x.gaiseki(init_p2_state_G[0]).normalise();
  // 変換行列を生成
  std::array<std::array<ScalarType, 3>, 3> convert_G_to_R1 =
      create_rodrigues_matrix(rotax1, -theta1);
  std::array<std::array<ScalarType, 3>, 3> rotate1 = create_rodrigues_matrix(rotax1, theta1);

  /* 地球の軌道面の法線ベクトルと中間座標系rot1のz軸を一致させ目的の回転座標系rot2へ移す
   */
  // 地球の軌道面の法線ベクトル
  Vector3d<ScalarType> n_plane = init_p2_state_G[0].gaiseki(init_p2_state_G[1]).normalise();

  // 慣性系から見た中間座標系rot1のz軸
  Vector3d<ScalarType> rot1_z = ApplyMatrix(rotate1, G_z);

  //   慣性座標系から見た地球軌道面の法線ベクトルと慣性系から見た中間座標系rot1のz軸がなす角度を計算
  ScalarType theta2 = std::atan2(n_plane.gaiseki(rot1_z).magnitude(), n_plane.naiseki(rot1_z));
  Vector3d<ScalarType> rotax2 = {1.0, 0.0, 0.0};

  // 変換行列を生成
  std::array<std::array<ScalarType, 3>, 3> convert_R1_to_R2 =
      create_rodrigues_matrix(rotax2, -theta2);
  std::array<std::array<ScalarType, 3>, 3> rotate2 = create_rodrigues_matrix(rotax2, theta2);
  // 小惑星の位置を回転座標系に変換
  Vector3d<ScalarType> ast_pos_R1 = ApplyMatrix(convert_G_to_R1, init_ast_state_G[0]);
  Vector3d<ScalarType> epos_R1 = ApplyMatrix(convert_G_to_R1, init_p2_state_G[0]);
  Vector3d<ScalarType> ast_pos_R2 = ApplyMatrix(convert_R1_to_R2, ast_pos_R1);
  Vector3d<ScalarType> epos_R2 = ApplyMatrix(convert_R1_to_R2, epos_R1);
  Vector3d<ScalarType> ast_pos_R{ast_pos_R2.x() + 1.0 - epos_R2.x() - astro_params.mu,
                                 ast_pos_R2.y(), ast_pos_R2.z()};

  /* calc the velocity of rotating frame */
  //  変換行列convert_G_to_R1の微分を計算
  // theta1の微分

  ScalarType ND_time_ref =
      std::sqrt(astro_params.au * astro_params.au / astro_params.gm_sun);  // Non-Dimensional time
  ScalarType mean_motion = 1 / ND_time_ref;
  std::cout << "mean motion = " << mean_motion << std::endl;
  // conbert AU/day to Non-Dimensional
  Vector3d<ScalarType> ND_init_e_velo_G = n_plane.gaiseki(init_p2_state_G[0]).normalise() *
                                          init_p2_state_G[1].magnitude() * ND_time_ref / 86400.;

  ScalarType e_G_x = init_p2_state_G[0].x();
  ScalarType e_G_y = init_p2_state_G[0].y();
  ScalarType e_G_z = init_p2_state_G[0].z();
  ScalarType ND_e_G_vx = ND_init_e_velo_G.x();
  ScalarType ND_e_G_vy = ND_init_e_velo_G.y();
  ScalarType ND_e_G_vz = ND_init_e_velo_G.z();

  ScalarType poo = init_p2_state_G[0].magnitude();
  ScalarType sin_theta1 = std::sqrt(1. - std::cos(theta1) * std::cos(theta1));

  // ScalarType dot_theta1 = -(ND_e_G_vx * poo * poo -
  //                      e_G_x * (e_G_x * ND_e_G_vx + e_G_y * ND_e_G_vy +
  //                      e_G_z * ND_e_G_vz)) /
  //                      <-
  //                      円制限三体問題だからちきゅうの速度と位置の内積は0
  //                     poo / poo / poo / sin_theta1;
  ScalarType dot_theta1 = -ND_e_G_vx / poo / sin_theta1;

  // theta2の微分
  ScalarType sin_theta2 = std::sqrt(1. - std::cos(theta2) * std::cos(theta2));
  if (n_plane.gaiseki(ApplyMatrix(rotate1, G_x)).naiseki(rot1_z) < 0) {
    sin_theta2 = -sin_theta2;
  }

  Vector3d<ScalarType> k{n_plane.x(), n_plane.y(), n_plane.z()};

  ScalarType n1 = std::sqrt(e_G_y * e_G_y + e_G_z * e_G_z);
  ScalarType a = e_G_y * ND_e_G_vy + e_G_z * ND_e_G_vz;
  ScalarType b = ND_e_G_vy * e_G_z + e_G_y * ND_e_G_vz;

  ScalarType dot_theta2 = (k.x() * (ND_e_G_vz * n1 * n1 - e_G_z * a) * sin_theta1 / n1 / n1 / n1 -
                           k.x() * e_G_z * std ::cos(theta1) * dot_theta1 / n1 +
                           k.y() * (b * n1 * n1 - 2. * e_G_z * e_G_y * a) *
                               (1. - std::cos(theta1)) / n1 / n1 / n1 / n1 +
                           k.y() * e_G_z * ND_e_G_vz * sin_theta1 * dot_theta1 / n1 / n1 -
                           2. * k.z() * (e_G_y * ND_e_G_vy * n1 * n1 - e_G_y * a) *
                               (1. - std::cos(theta1)) / n1 / n1 / n1 / n1 -
                           k.z() * e_G_y * e_G_y * sin_theta1 * dot_theta1 / n1 / n1 +
                           k.z() * sin_theta1 * dot_theta1) /
                          sin_theta2;

  ScalarType c1 = G_x.naiseki(init_p2_state_G[0]) / init_p2_state_G[0].magnitude();
  ScalarType s1 = sin_theta1;
  ScalarType c2 = std::cos(theta2);
  ScalarType s2 = sin_theta2;
  std::array<std::array<ScalarType, 3>, 3> dot_convert_G_to_R1{
      {{-s1 * dot_theta1,
        (ND_e_G_vy * n1 * n1 - e_G_y * a) * s1 / n1 / n1 / n1 + e_G_y * c1 * dot_theta1 / n1,
        (ND_e_G_vz * n1 * n1 - e_G_z * a) * s1 / n1 / n1 / n1 + e_G_z * c1 * dot_theta1 / n1},
       {-(ND_e_G_vy * n1 * n1 - e_G_y * a) * s1 / n1 / n1 / n1 - e_G_y * c1 * dot_theta1 / n1,
        (2. * e_G_z * ND_e_G_vz * n1 * n1 - 2. * e_G_z * e_G_z * a) * (1.0 - c1) / n1 / n1 / n1 /
                n1 +
            e_G_z * e_G_z * s1 * dot_theta1 / n1 / n1 - s1 * dot_theta1,
        -(b * n1 * n1 - 2. * e_G_z * e_G_y * a) * (1. - c1) / n1 / n1 / n1 / n1 -
            e_G_z * e_G_y * s1 * dot_theta1 / n1 / n1},
       {-(ND_e_G_vz * n1 * n1 - e_G_z * a) * s1 / n1 / n1 / n1 - e_G_z * c1 * dot_theta1 / n1,
        -(b * n1 * n1 - 2. * e_G_z * e_G_y * a) * (1.0 - c1) / n1 / n1 / n1 / n1 -
            e_G_z * e_G_y * s1 * dot_theta1 / n1 / n1,
        (2. * e_G_y * ND_e_G_vy * n1 * n1 - e_G_y * e_G_y * a) * (1. - c1) / n1 / n1 / n1 / n1 +
            e_G_y * e_G_y * s1 * dot_theta1 / n1 / n1 - s1 * dot_theta1}}};

  std::array<std::array<ScalarType, 3>, 3> dot_convert_R1_to_R2{
      {{0.0, 0.0, 0.0},
       {0.0, -s2 * dot_theta2, c2 * dot_theta2},
       {0.0, -c2 * dot_theta2, -s2 * dot_theta2}}};
  // 速度を座標変換
  Vector3d<ScalarType> ast_vel_R_1 =
      ApplyMatrix(dot_convert_R1_to_R2, ApplyMatrix(convert_G_to_R1, init_ast_state_G[0]));
  Vector3d<ScalarType> ast_vel_R_2 =
      ApplyMatrix(convert_R1_to_R2, ApplyMatrix(dot_convert_G_to_R1, init_ast_state_G[0]));
  Vector3d<ScalarType> ast_vel_R_3 =
      ApplyMatrix(convert_R1_to_R2, ApplyMatrix(convert_G_to_R1, init_ast_state_G[1]));

  Vector3d<ScalarType> ast_vel_R = ast_vel_R_1 + ast_vel_R_2 + ast_vel_R_3;
  // check
  Vector3d<ScalarType> e_vel_R_1 = {0, 0, 0};
  // ApplyMatrix(dot_convert_R1_to_R2, ApplyMatrix(convert_G_to_R1,
  // init_p2_state_G[0]));
  Vector3d<ScalarType> e_vel_R_2 = {0, 0, 0};
  // ApplyMatrix(convert_R1_to_R2, ApplyMatrix(dot_convert_G_to_R1,
  // init_p2_state_G[0]));
  Vector3d<ScalarType> e_vel_R_3 =
      ApplyMatrix(convert_R1_to_R2, ApplyMatrix(convert_G_to_R1, ND_init_e_velo_G));
  Vector3d<ScalarType> e_vel_R = e_vel_R_1 + e_vel_R_2 + e_vel_R_3;
  // check
  std::cout << "converted earth velocity = " << e_vel_R.x() << " " << e_vel_R.y() << " "
            << e_vel_R.z() << std::endl;
  std::cout << "<>                 Asteroid data converted to rotating frame" << std::endl;
  if constexpr (std::is_same_v<ElementType, ScalarType> && Dim == 6) {
    return {ast_pos_R.x(), ast_pos_R.y(), ast_pos_R.z(),
            ast_vel_R.x(), ast_vel_R.y(), ast_vel_R.z()};
  } else if constexpr (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2) {
    return {ast_pos_R, ast_vel_R};
  }
}

/**
 * @brief J2000慣性系（太陽中心, AU, AU/day）の状態ベクトルを、
 * 太陽-P2回転座標系（共通重心中心, 無次元）に変換
 *
 * この関数は、P2（例：地球）の一般的な（傾斜した）軌道を扱う
 * 回転座標系 (R) は以下のように定義：
 * - 原点: 太陽-P2 共通重心
 * - x軸: 太陽からP2を指す方向
 * - z軸: P2の軌道角運動量ベクトル (h = r_p2 x v_p2) に平行
 * - y軸: 右手系を完成させる (y = z x x)
 *
 * 出力は、長さ単位(LU)=r_p2、時間単位(TU)=1/n で無次元化。
 *
 * @tparam ElementType std::arrayの要素型 (ScalarType または Vector3d<ScalarType>)。
 * @tparam ScalarType 計算に使用するスカラー型 (e.g., double)。
 * @tparam Dim std::arrayのサイズ (6 または 2)。
 * @param ast_state J2000慣性系（太陽中心）における小惑星の状態ベクトル
 * [pos(AU), vel(AU/day)]。
 * @param p2_state J2000慣性系（太陽中心）におけるP2（地球）の状態ベクトル
 * [pos(AU), vel(AU/day)]。
 * @param astro_params 天文定数（au(m), gm_sun(m^3/s^2), gm_earth(m^3/s^2)）。
 * @return 回転座標系（無次元）における小惑星の状態ベクトル [pos_nd, vel_nd]。
 */
template <typename ElementType, typename ScalarType, std::size_t Dim>
const std::array<ElementType, Dim> ConvertInertial2RotatingV2(
    const std::array<ElementType, Dim>& ast_state, const std::array<ElementType, Dim>& p2_state,
    const AstroConstants<ScalarType>& astro_params) {
  // 正しいテンプレート引数か静的に検証
  static_assert((std::is_same_v<ElementType, ScalarType> && Dim == 6) ||
                    (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2),
                "Input must be std::array<ScalarType, 6> or std::array<Vector3d<ScalarType>, 2>");

  // 0. 単位換算 (m^3/s^2 -> AU^3/day^2)
  if (astro_params.au == 0) {
    throw std::runtime_error("AstroConstants::au cannot be zero.");
  }
  const ScalarType kM3s2_to_AU3day2 =
      (kSecondsPerDay * kSecondsPerDay) / (astro_params.au * astro_params.au * astro_params.au);
  const ScalarType gm_sun_AD = astro_params.gm_sun * kM3s2_to_AU3day2;
  const ScalarType gm_earth_AD = astro_params.gm_earth * kM3s2_to_AU3day2;
  const ScalarType gm_total_AD = gm_sun_AD + gm_earth_AD;

  // 1. 入力を {pos, vel} の Vector3d<ScalarType> ペアに統一
  //   (G = J2000 Inertial Frame, 単位: AU, AU/day)
  StateVectorPair<ScalarType> ast_state_G;
  StateVectorPair<ScalarType> p2_state_G;
  if constexpr (std::is_same_v<ElementType, ScalarType> && Dim == 6) {
    ast_state_G = ConvertStateToPair(ast_state);
    p2_state_G = ConvertStateToPair(p2_state);
  } else if constexpr (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2) {
    ast_state_G = {ast_state[0], ast_state[1]};
    p2_state_G = {p2_state[0], p2_state[1]};
  }

  const Vector3d<ScalarType> r_ast_G = ast_state_G[0];
  const Vector3d<ScalarType> v_ast_G = ast_state_G[1];
  const Vector3d<ScalarType> r_p2_G = p2_state_G[0];
  const Vector3d<ScalarType> v_p2_G = p2_state_G[1];

  // 2. 質量比と共通重心 (Barycenter, B) の状態を計算
  // 質量比は単位系に依存しない
  const ScalarType mu = astro_params.gm_earth / (astro_params.gm_sun + astro_params.gm_earth);

  // 太陽（G系の原点）から見た共通重心の状態 (AU, AU/day)
  const Vector3d<ScalarType> r_B_G = r_p2_G * mu;
  const Vector3d<ScalarType> v_B_G = v_p2_G * mu;

  // 3. G系における回転座標系 (R) の基底ベクトルを定義
  const ScalarType r_p2_mag = r_p2_G.magnitude();  // これが LU (AU)
  if (r_p2_mag == 0) {
    throw std::runtime_error("P2 position magnitude is zero.");
  }
  const Vector3d<ScalarType> i_R_G = r_p2_G / r_p2_mag;  // x軸 (太陽 -> P2)

  const Vector3d<ScalarType> h_G = r_p2_G.gaiseki(v_p2_G);  // P2の軌道角運動量 (AU^2/day)
  const ScalarType h_mag = h_G.magnitude();
  if (h_mag == 0) {
    throw std::runtime_error("P2 angular momentum magnitude is zero (e.g., linear motion).");
  }
  const Vector3d<ScalarType> k_R_G = h_G / h_mag;  // z軸 (h の方向)

  const Vector3d<ScalarType> j_R_G = k_R_G.gaiseki(i_R_G);  // y軸 (右手系を完成)

  // 4. 回転行列 C (G系 -> R系) を作成
  const Matrix3x3<ScalarType> C_G2R = {{{i_R_G.x(), i_R_G.y(), i_R_G.z()},
                                        {j_R_G.x(), j_R_G.y(), j_R_G.z()},
                                        {k_R_G.x(), k_R_G.y(), k_R_G.z()}}};

  // 5. 位置ベクトルの変換 (単位: AU)
  const Vector3d<ScalarType> r_ast_rel_B_G = r_ast_G - r_B_G;
  const Vector3d<ScalarType> ast_pos_R = ApplyMatrix(C_G2R, r_ast_rel_B_G);  // (AU)

  // 6. 回転系の角速度ベクトル (omega_R) を計算 (単位: rad/day)
  // P2の加速度 (AU/day^2)
  const Vector3d<ScalarType> a_p2_G = r_p2_G * (-gm_sun_AD / (r_p2_mag * r_p2_mag * r_p2_mag));

  const ScalarType dot_r_p2_mag = r_p2_G.naiseki(v_p2_G) / r_p2_mag;  // (AU/day)
  const Vector3d<ScalarType> dot_i_R_G =
      (v_p2_G / r_p2_mag) - (r_p2_G * (dot_r_p2_mag / (r_p2_mag * r_p2_mag)));  // (1/day)

  const Vector3d<ScalarType> dot_h_G = r_p2_G.gaiseki(a_p2_G);  // (AU^2/day^2)
  const ScalarType dot_h_mag = h_G.naiseki(dot_h_G) / h_mag;    // (AU^2/day^2)
  const Vector3d<ScalarType> dot_k_R_G =
      (dot_h_G / h_mag) - (h_G * (dot_h_mag / (h_mag * h_mag)));  // (1/day)

  const Vector3d<ScalarType> dot_j_R_G =
      dot_k_R_G.gaiseki(i_R_G) + k_R_G.gaiseki(dot_i_R_G);  // (1/day)

  const ScalarType omega_x_R = dot_j_R_G.naiseki(k_R_G);  // (rad/day)
  const ScalarType omega_y_R = dot_k_R_G.naiseki(i_R_G);  // (rad/day)
  const ScalarType omega_z_R = dot_i_R_G.naiseki(j_R_G);  // (rad/day)

  const Vector3d<ScalarType> omega_R(omega_x_R, omega_y_R, omega_z_R);  // (rad/day)

  // 7. 速度ベクトルの変換 (単位: AU/day)
  const Vector3d<ScalarType> v_ast_rel_B_G = v_ast_G - v_B_G;                    // (AU/day)
  const Vector3d<ScalarType> v_R_transport = ApplyMatrix(C_G2R, v_ast_rel_B_G);  // (AU/day)
  const Vector3d<ScalarType> v_R_coriolis = omega_R.gaiseki(ast_pos_R);          // (rad/day * AU)
  const Vector3d<ScalarType> ast_vel_R = v_R_transport - v_R_coriolis;           // (AU/day)

  // 8. 無次元化
  // 長さの単位 (LU)
  const ScalarType lu = r_p2_mag;  // (AU)

  // 時間の単位 (TU)
  // n = sqrt( (GM_sun + GM_earth) / LU^3 )  (rad/day)
  const ScalarType mean_motion_AD = std::sqrt(gm_total_AD / (lu * lu * lu));
  if (mean_motion_AD == 0) {
    throw std::runtime_error("Mean motion is zero, cannot non-dimensionalize.");
  }
  const ScalarType tu = 1.0 / mean_motion_AD;  // (day)

  // 速度の単位 (VU)
  const ScalarType vu = lu / tu;  // (AU/day)

  // 無次元化の実行
  const Vector3d<ScalarType> ast_pos_R_nd = ast_pos_R / lu;
  const Vector3d<ScalarType> ast_vel_R_nd = ast_vel_R / vu;

  // 9. 出力のフォーマットを整える
  if constexpr (std::is_same_v<ElementType, ScalarType> && Dim == 6) {
    return ConvertPairToState({ast_pos_R_nd, ast_vel_R_nd});
  } else if constexpr (std::is_same_v<ElementType, Vector3d<ScalarType>> && Dim == 2) {
    return {ast_pos_R_nd, ast_vel_R_nd};
  }
}

// CRTBPの運動方程式
template <typename ScalarType>
class EquationOfMotion {
 private:
  AstroConstants<ScalarType> astro_params_;

 public:
  explicit EquationOfMotion(const AstroConstants<ScalarType>& astro_params)
      : astro_params_(astro_params) {}

  void operator()(const StateVector6<ScalarType>& x, StateVector6<ScalarType>& dxdt,
                  const ScalarType /* t */) const {
    const ScalarType x_ = x[0], y_ = x[1], z_ = x[2];
    const ScalarType vx_ = x[3], vy_ = x[4], vz_ = x[5];

    const ScalarType r1 = calc_r1<ScalarType>(x, astro_params_);
    const ScalarType r2 = calc_r2<ScalarType>(x, astro_params_);

    dxdt[0] = vx_;
    dxdt[1] = vy_;
    dxdt[2] = vz_;
    dxdt[3] = 2. * vy_ + x_ - (1. - astro_params_.mu) * (x_ + astro_params_.mu) / (r1 * r1 * r1) -
              astro_params_.mu * (x_ - 1. + astro_params_.mu) / (r2 * r2 * r2);
    dxdt[4] = -2. * vx_ + y_ - (1. - astro_params_.mu) * y_ / (r1 * r1 * r1) -
              astro_params_.mu * y_ / (r2 * r2 * r2);
    dxdt[5] =
        -(1. - astro_params_.mu) * z_ / (r1 * r1 * r1) - astro_params_.mu * z_ / (r2 * r2 * r2);
  }
};
template <typename ScalarType>
class observer {
 private:
  ScalarType init_jacobi_;
  AstroConstants<ScalarType> astro_params_;

 public:
  std::vector<std::array<ScalarType, 8>>& m_out;
  explicit observer(ScalarType init_jacobi, const AstroConstants<ScalarType>& astro_params,
                    std::vector<std::array<ScalarType, 8>>& out)
      : init_jacobi_(init_jacobi), astro_params_(astro_params), m_out(out) {}

  void operator()(const StateVector6<ScalarType>& x, ScalarType t) const {
    // T_ jacobi_integral = calc_jacobi_integral(x, astro_params_); //
    // コメントアウトされていた関数 m_out.push_back({t, x[0], x[1], x[2], x[3], x[4], x[5],
    // init_jacobi_ - jacobi_integral});

    // calc_jacobi_integral が未定義のため、ダミーの値をプッシュ
    m_out.push_back({t, x[0], x[1], x[2], x[3], x[4], x[5], 0.0});
  }
};

namespace {

template <typename ScalarType>
struct CanonicalState {
  Vector3d<ScalarType> q;
  Vector3d<ScalarType> p;

  ScalarType& x() { return q.x; }
  const ScalarType& x() const { return q.x; }

  ScalarType& y() { return q.y; }
  const ScalarType& y() const { return q.y; }

  // z の読み書き
  ScalarType& z() { return q.z; }
  const ScalarType& z() const { return q.z; }

  ScalarType& px() { return p.x; }
  const ScalarType& px() const { return p.x; }

  ScalarType& py() { return p.y; }
  const ScalarType& py() const { return p.y; }

  ScalarType& pz() { return p.z; }
  const ScalarType& pz() const { return p.z; }
};

/**
 * @brief 物理状態 (q, q_dot) を正準状態 (q, p) に変換します。
 * p_x = vx + y
 * p_y = vy - x
 * p_z = vz
 */
template <typename ScalarType>
CanonicalState<ScalarType> ConvertToCanonical(const State<ScalarType>& state) {
  return CanonicalState{state.x, state.y, state.z, state.vx + state.y, state.vy - state.x,
                        state.vz};
}

/**
 * @brief 正準状態 (q, p) を物理状態 (q, q_dot) に変換します。
 * vx = p_x - y
 * vy = p_y + x
 * vz = p_z
 */
template <typename ScalarType>
State<ScalarType> ConvertToPhysical(const CanonicalState<ScalarType>& canonical_state) {
  return State{canonical_state.x,
               canonical_state.y,
               canonical_state.z,
               canonical_state.px - canonical_state.y,
               canonical_state.py + canonical_state.x,
               canonical_state.pz};
}

/**
 * @brief ハミルトニアン H_B (コリオリ項) による流れを適用します。
 * これは (x, y) 平面と (px, py) 平面における角度 `angle` の
 * 反時計回りの回転として解析的に解けます。
 * @param state 更新対象の正準状態 (in-out)。
 * @param angle 回転角 (ラジアン)。
 */
template <typename ScalarType>
void ApplyRotation(CanonicalState<ScalarType>* state, ScalarType angle) {
  const ScalarType c = std::cos(angle);
  const ScalarType s = std::sin(angle);

  // (x, y) の回転
  const ScalarType x_new = state->x * c - state->y * s;
  const ScalarType y_new = state->x * s + state->y * c;
  state->x = x_new;
  state->y = y_new;

  // (px, py) の回転
  const ScalarType px_new = state->px * c - state->py * s;
  const ScalarType py_new = state->px * s + state->py * c;
  state->px = px_new;
  state->py = py_new;
}

/**
 * @brief ポテンシャル V(q) の勾配 grad_V を計算します。
 * V = (1-mu)/r1 + mu/r2
 * grad_V = (dV/dx, dV/dy, dV/dz)
 * @param params 質量比 mu を含むパラメータ。
 * @param x, y, z 位置座標。
 * @param[out] grad_V 計算された勾配ベクトル (x, y, z 成分)。
 */
template <typename ScalarType>
void CalculateGradientV(const AstroConstants<ScalarType>& params, ScalarType x, ScalarType y,
                        ScalarType z, ScalarType* grad_V) {
  constexpr double kMinDistanceSq = 1e-10;  // 最小距離の二乗

  const ScalarType mu = params.gm_earth / (params.gm_earth + params.gm_sun);
  const ScalarType mu1 = 1.0 - mu;

  const ScalarType x1 = x + mu;
  const ScalarType x2 = x - (1.0 - mu);

  const ScalarType r1_sq = x1 * x1 + y * y + z * z;
  const ScalarType r2_sq = x2 * x2 + y * y + z * z;

  if (r1_sq < kMinDistanceSq) {
    throw std::runtime_error(
        "Position too close to primary 1 (m1): r1^2 = " + std::to_string(r1_sq) + " at (" +
        std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")");
  }
  if (r2_sq < kMinDistanceSq) {
    throw std::runtime_error(
        "Position too close to primary 2 (m2): r2^2 = " + std::to_string(r2_sq) + " at (" +
        std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")");
  }

  // r1^(-3) と r2^(-3)
  const ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
  const ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

  grad_V[0] = -mu1 * x1 * r1_inv3 - mu * x2 * r2_inv3;
  grad_V[1] = -mu1 * y * r1_inv3 - mu * y * r2_inv3;
  grad_V[2] = -mu1 * z * r1_inv3 - mu * z * r2_inv3;
}

/**
 * @brief H_A の Kick (運動量更新) ステップを適用します。
 * p = p - dt * grad_V(q)
 * @param params 質量比 mu。
 * @param state 更新対象の正準状態 (in-out)。
 * @param dt 時間ステップ幅。
 */
template <typename ScalarType>
void ApplyKick(const AstroConstants<ScalarType>& params, CanonicalState<ScalarType>* state,
               ScalarType dt) {
  double grad_V[3];
  CalculateGradientV(params, state->x, state->y, state->z, grad_V);

  state->px -= dt * grad_V[0];
  state->py -= dt * grad_V[1];
  state->pz -= dt * grad_V[2];
}

/**
 * @brief H_A の Drift (位置更新) ステップを適用します。
 * q = q + dt * p  (注意: H_A の流れでは q_dot = p)
 * @param state 更新対象の正準状態 (in-out)。
 * @param dt 時間ステップ幅。
 */
template <typename ScalarType>
void ApplyDrift(CanonicalState<ScalarType>* state, ScalarType dt) {
  state->x += dt * state->px;
  state->y += dt * state->py;
  state->z += dt * state->pz;
}

}  // namespace

namespace cr3bp {

// ( ... 前回のヘルパー関数群 ... )

/**
 * @brief 2次のシンプレクティック積分ステップ
 *
 * S2(h) = Phi_B(h/2) * Psi_A(h) * Phi_B(h/2)
 *
 * @param params 質量比 mu を含むパラメータ。
 * @param state 積分開始時の物理状態。
 * @param h 1ステップの時間幅。
 * @return 積分後の物理状態 (t + h)。
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep(const AstroConstants<ScalarType>& params,
                                 const State<ScalarType>& state, ScalarType h) {
  CanonicalState<ScalarType> canonical_state = ConvertToCanonical(state);

  // 1. B(h/2) - 回転
  ApplyRotation(&canonical_state, h / 2.0);

  // 2. A(h) = Kick(h/2) * Drift(h) * Kick(h/2)
  // (注意: 符号を修正した「+=」バージョン)
  ApplyKick(params, &canonical_state, h / 2.0);
  ApplyDrift(&canonical_state, h);
  ApplyKick(params, &canonical_state, h / 2.0);

  // 3. B(h/2) - 回転
  ApplyRotation(&canonical_state, h / 2.0);

  return ConvertToPhysical(canonical_state);
}

/**
 * @brief 4次のシンプレクティック積分ステップ (吉田法, 1990)
 *
 * 画像の式 (48) に基づくコンポジション法 (トリプルジャンプ)
 * S4(tau) = S2(x1 * tau) * S2(x0 * tau) * S2(x1 * tau)
 *
 * @param params 質量比 mu を含むパラメータ。
 * @param state 積分開始時の物理状態。
 * @param tau 1ステップの時間幅。
 * @return 積分後の物理状態 (t + tau)。
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep4thOrder(const AstroConstants<ScalarType>& params,
                                         const State<ScalarType>& state, ScalarType tau) {
  // 吉田 (1990) による4次積分のための係数
  // x0 = -2^(1/3) / (2 - 2^(1/3))
  // x1 = 1 / (2 - 2^(1/3))
  const ScalarType kX1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
  const ScalarType kX0 = 1.0 - 2.0 * kX1;  // (x0 + 2*x1 = 1 の関係から)

  // S4(tau) = S2(x1 * tau) * S2(x0 * tau) * S2(x1 * tau)
  State<ScalarType> state1 = SymplecticStep(params, state, kX1 * tau);
  State<ScalarType> state2 = SymplecticStep(params, state1, kX0 * tau);
  State<ScalarType> state3 = SymplecticStep(params, state2, kX1 * tau);

  return state3;
}
template <typename ScalarType>
class observer {
 private:
  double init_jacobi_;
  AstroConstants<ScalarType> astro_params_;

 public:
  std::vector<std::array<ScalarType, 8>>& m_out;
  explicit observer(ScalarType init_jacobi, const AstroConstants<ScalarType>& astro_params,
                    std::vector<std::array<ScalarType, 8>>& out)
      : init_jacobi_(init_jacobi), astro_params_(astro_params), m_out(out) {}

  void operator()(const State<ScalarType>& x, ScalarType t) const {
    double jacobi_integral = calc_jacobi_integral(x, astro_params_);
    m_out.push_back({t, x[0], x[1], x[2], x[3], x[4], x[5], init_jacobi_ - jacobi_integral});
  }
};
/**
 * @brief 4次積分器を使用して指定された時間、積分を実行します。
 */
template <typename ScalarType>
std::vector<State<ScalarType>> Integrate4thOrder(const AstroConstants<ScalarType>& params,
                                                 const State<ScalarType>& initial_state,
                                                 ScalarType time_step, int num_steps, void*) {
  std::vector<State<ScalarType>> trajectory;
  trajectory.reserve(num_steps);

  State<ScalarType> current_state = initial_state;
  for (int i = 0; i < num_steps; ++i) {
    // ここで4次ステップ関数を呼び出す
    current_state = SymplecticStep4thOrder(params, current_state, time_step);
    trajectory.push_back(current_state);
  }
  return trajectory;
}

// class Trajectory {
// private:
//   double jacobi_integral_;
//   AstroConstants astro_params_;
//
// public:
//   std::vector<std::array<double, 8>> trajectory_;
//
//   Trajectory(double jacobi_integral, const AstroConstants &astro_params)
//       : jacobi_integral_(jacobi_integral), astro_params_(astro_params) {}
//
//   std::vector<std::array<double, 8>> integrate(const state &x0, double t0,
//                                                double tf, double dt) {
//     state init_state = x0;
//     trajectory_.clear();
//     trajectory_.reserve(static_cast<size_t>((tf - t0) / dt) + 1);
//
//     auto eom = EquationOfMotion(astro_params_);
//     auto obs = observer(jacobi_integral_, astro_params_, trajectory_);
//     boost::numeric::odeint::controlled_runge_kutta<
//         boost::numeric::odeint::runge_kutta_fehlberg78<state>>
//         stepper;
//
//     boost::numeric::odeint::integrate_const(stepper, eom, init_state, t0,
//     tf,
//                                               dt, obs);
//
//     return trajectory_;
//   }
// };

}  // namespace crtbp
#endif  // RTBP_HPP
