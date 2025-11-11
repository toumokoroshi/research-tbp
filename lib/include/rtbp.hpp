/**
 * @file rtbp.hpp
 * @author tabata (modified by Gemini)
 * @brief 制限三体問題の関数をまとめたライブラリ (SALI対応)
 * @version 2.0 (SALI対応シンプレクティック積分器 修正版)
 * @date 2025-11-10
 * @note C++17
 * @par history
 * - 2025-01-24 tabata version 1.0
 * - 2025-11-06 refactored template parameter names
 * - 2025-11-10 Gemini - Logic corrections, struct refactoring
 */
#ifndef RTBP_HPP
#define RTBP_HPP
#include <algorithm>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
// #include <numeric>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector3d.hpp>
#include <vector>
namespace my_type {
// 3D座標を表す型のエイリアス
template <typename ScalarType>
struct State3d {
  ScalarType x = 0.0, y = 0.0, z = 0.0;

  State3d() = default;
  State3d(ScalarType x, ScalarType y, ScalarType z) : x(x), y(y), z(z) {}
};
// 3x3 Matrix type alias
template <typename ScalarType>
using Matrix3x3 = std::array<std::array<ScalarType, 3>, 3>;

/**
 * @brief 物理状態ベクトル (位置, 速度)
 * (x, y, z, vx, vy, vz)
 */
template <typename ScalarType>
struct State {
  ScalarType x = 0.0, y = 0.0, z = 0.0;
  ScalarType vx = 0.0, vy = 0.0, vz = 0.0;

  // 6スカラーでのコンストラクタ
  State() = default;
  State(ScalarType x, ScalarType y, ScalarType z, ScalarType vx, ScalarType vy, ScalarType vz)
      : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {}
  // RK4のための演算子
  State<ScalarType> operator+(const State<ScalarType>& other) const {
    return {x + other.x, y + other.y, z + other.z, vx + other.vx, vy + other.vy, vz + other.vz};
  }

  State<ScalarType> operator*(ScalarType scalar) const {
    return {x * scalar, y * scalar, z * scalar, vx * scalar, vy * scalar, vz * scalar};
  }
};

/**
 * @brief 正準状態ベクトル (位置 q, 正準運動量 p)
 * (qx, qy, qz, px, py, pz)
 * 偏差ベクトル w = (dq, dp) の表現にも使用
 */
template <typename ScalarType>
struct CanonicalState {
  ScalarType qx = 0.0, qy = 0.0, qz = 0.0;
  ScalarType px = 0.0, py = 0.0, pz = 0.0;

  // 6スカラーでのコンストラクタ
  CanonicalState() = default;
  CanonicalState(ScalarType qx, ScalarType qy, ScalarType qz, ScalarType px, ScalarType py,
                 ScalarType pz)
      : qx(qx), qy(qy), qz(qz), px(px), py(py), pz(pz) {}

  // アクセサ (const /非const)
  ScalarType& x() { return qx; }
  const ScalarType& x() const { return qx; }
  ScalarType& y() { return qy; }
  const ScalarType& y() const { return qy; }
  ScalarType& z() { return qz; }
  const ScalarType& z() const { return qz; }

  // 演算子
  CanonicalState<ScalarType> operator+(const CanonicalState<ScalarType>& other) const {
    return {qx + other.qx, qy + other.qy, qz + other.qz,
            px + other.px, py + other.py, pz + other.pz};
  }

  CanonicalState<ScalarType> operator-(const CanonicalState<ScalarType>& other) const {
    return {qx - other.qx, qy - other.qy, qz - other.qz,
            px - other.px, py - other.py, pz - other.pz};
  }

  CanonicalState<ScalarType> operator/(ScalarType scalar) const {
    return {qx / scalar, qy / scalar, qz / scalar, px / scalar, py / scalar, pz / scalar};
  }

  // 6次元位相空間でのノルム (SALI計算用)
  ScalarType Norm() const {
    return std::sqrt(qx * qx + qy * qy + qz * qz + px * px + py * py + pz * pz);
  }

  // 正規化 (偏差ベクトルの長さを1に戻す)
  void Normalize() {
    ScalarType norm = Norm();
    if (norm > 1e-12) {  // ゼロ除算を回避
      qx /= norm;
      qy /= norm;
      qz /= norm;
      px /= norm;
      py /= norm;
      pz /= norm;
    }
  }
};
/**
 * @brief ポテンシャル V(q) のヘッセ行列 (2階微分)
 * H_ij = d^2 V / (d_qi d_qj)
 */
template <typename ScalarType>
struct HessianMatrix {
  // H = [[hxx, hxy, hxz],
  //      [hxy, hyy, hyz],
  //      [hxz, hyz, hzz]]
  ScalarType hxx = 0.0, hyy = 0.0, hzz = 0.0;
  ScalarType hxy = 0.0, hxz = 0.0, hyz = 0.0;
};

/**
 * @brief SALI計算に必要な全ての状態を保持する構造体
 */
template <typename ScalarType>
struct SaliState {
  CanonicalState<ScalarType> state;  // 主軌道の状態 (q, p)
  CanonicalState<ScalarType> w1;     // 偏差ベクトル1 (dq1, dp1)
  CanonicalState<ScalarType> w2;     // 偏差ベクトル2 (dq2, dp2)
};
};  // namespace my_type
/**
 * @brief スタティックな関数
 */
namespace {

constexpr double kMinDistanceSq = 1e-16;
constexpr double kSecondsPerDay = 86400.0;

using namespace my_type;

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
/**
 * @brief 円制限三体問題の系で物理状態 (q, q_dot) を正準状態 (q, p) に変換
 * p_x = vx + y
 * p_y = vy - x
 * p_z = vz
 */
template <typename ScalarType>
CanonicalState<ScalarType> ConvertToCanonical(const State<ScalarType>& state) {
  return CanonicalState<ScalarType>{
      state.x,
      state.y,
      state.z,
      state.vx + state.y,  // px = vx + y
      state.vy - state.x,  // py = vy - x
      state.vz             // pz = vz
  };
}

/**
 * @brief 円制限三体問題の系で正準状態 (q, p) を物理状態 (q, q_dot) に変換
 * vx = p_x - y
 * vy = p_y + x
 * vz = p_z
 */
template <typename ScalarType>
State<ScalarType> ConvertToPhysical(const CanonicalState<ScalarType>& canonical_state) {
  return State<ScalarType>{
      canonical_state.qx,
      canonical_state.qy,
      canonical_state.qz,
      canonical_state.px - canonical_state.qy,  // vx = px - y
      canonical_state.py + canonical_state.qx,  // vy = py + x
      canonical_state.pz                        // vz = pz
  };
}

// 距離r1 (第三体から m1 への距離) を計算
template <typename ScalarType>
inline ScalarType calc_r1(ScalarType x, ScalarType y, ScalarType z, ScalarType mu) {
  const ScalarType x1 = x + mu;
  return std::sqrt(std::pow(x1, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
}

// 距離r2 (第三体から m2 への距離) を計算
template <typename ScalarType>
inline ScalarType calc_r2(ScalarType x, ScalarType y, ScalarType z, ScalarType mu) {
  const ScalarType x2 = x - (1.0 - mu);
  return std::sqrt(std::pow(x - 1. + mu, 2) + y * y + z * z);
}
/**
 * @brief 有効ポテンシャル U* を計算
 * U* = 1/2(x^2+y^2) + (1-mu)/r1 + mu/r2
 */
template <typename ScalarType>
ScalarType calc_potential_U(ScalarType x, ScalarType y, ScalarType z, ScalarType mu) {
  ScalarType r1 = calc_r1(x, y, z, mu);
  ScalarType r2 = calc_r2(x, y, z, mu);
  if (r1 == 0.0 || r2 == 0.0) {
    throw std::runtime_error("Position coincides with a primary in calc_potential_U.");
  }
  return 0.5 * (x * x + y * y) + (1.0 - mu) / r1 + mu / r2 + mu * (1. - mu) * 0.5;
}

/**
 * @brief ヤコビ積分 C_J = 2*U* - v^2 を計算
 */
template <typename ScalarType>
ScalarType calc_jacobi_integral(const State<ScalarType>& state, const ScalarType mu) {
  const ScalarType v_sq = state.vx * state.vx + state.vy * state.vy + state.vz * state.vz;

  const ScalarType U_star = calc_potential_U(state.x, state.y, state.z, mu);
  return 2.0 * U_star - v_sq;
}

template <typename ScalarType>
ScalarType calc_v_abs(const State3d<ScalarType>& r, const ScalarType JACOBI_INTEGRAL,
                      const ScalarType mu) {
  ScalarType r1 = calc_r1<ScalarType>(r.x, r.y, r.z, mu);
  ScalarType r2 = calc_r2<ScalarType>(r.x, r.y, r.z, mu);

  ScalarType x = r.x;
  ScalarType y = r.y;
  ScalarType z = r.z;
  return std::sqrt(x * x + y * y + 2. * (1. - mu) / r1 + 2. * mu / r2 + mu * (1. - mu) -
                   JACOBI_INTEGRAL);
}

template <typename ScalarType>
State3d<ScalarType> calc_velocity(const State3d<ScalarType>& point, const ScalarType v_abs,
                                  const ScalarType mu, const ScalarType inclination,
                                  const ScalarType OMEGA, const ScalarType theta) {
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
  auto convert = [](Matrix3x3<ScalarType> convert_matrix,
                    const Vector3d<ScalarType>& v) -> State3d<ScalarType> {
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

  x = point.x;
  y = point.y;
  z = point.z;

  Vector3d<ScalarType> r2_vector{x - 1. + mu, y, z};
  Vector3d<ScalarType> h_vector = r2_vector.gaiseki(normal_vector);
  Vector3d<ScalarType> normalized_h_vector = h_vector.normalise();

  // 速度ベクトルを計算
  vx_ = v_abs * normalized_h_vector.x();
  vy_ = v_abs * normalized_h_vector.y();
  vz_ = v_abs * normalized_h_vector.z();

  if (theta == 0.0)
    return State3d<ScalarType>(vx_, vy_, vz_);
  else {
    // クオータニオンを用いて速度ベクトルをnormal_vector周りにthetaだけ回転させる
    std::array<std::array<ScalarType, 3>, 3> rot_matrix = create_rot_matrix(normal_vector, theta);
    Vector3d<ScalarType> velocity{vx_, vy_, vz_};
    State3d<ScalarType> rotated_velocity = convert(rot_matrix, velocity);

    return rotated_velocity;
  }
}

template <typename ElementType, typename ScalarType, std::size_t Dim>
const std::array<ElementType, Dim> ConvertInertial2Rotating____(
    const std::array<ElementType, Dim>& ast_state, const std::array<ElementType, Dim>& p2_state,
    const AstroConstants<ScalarType>& astro_params) {
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
template <typename ScalarType>
const State<ScalarType> ConvertInertial2RotatingV2(const State<ScalarType>& ast_state,
                                                   const State<ScalarType>& p2_state,
                                                   const AstroConstants<ScalarType>& astro_params) {
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
  State<ScalarType> ast_state_G = ast_state;
  State<ScalarType> p2_state_G = p2_state;

  const Vector3d<ScalarType> r_ast_G = {ast_state_G.x, ast_state_G.y, ast_state_G.z};
  const Vector3d<ScalarType> v_ast_G = {ast_state_G.vx, ast_state_G.vy, ast_state_G.vz};
  const Vector3d<ScalarType> r_p2_G = {p2_state_G.x, p2_state_G.y, p2_state_G.z};
  const Vector3d<ScalarType> v_p2_G = {p2_state_G.vx, p2_state_G.vy, p2_state_G.vz};
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

  return {ast_pos_R_nd.x(), ast_pos_R_nd.y(), ast_pos_R_nd.z(),
          ast_vel_R_nd.x(), ast_vel_R_nd.y(), ast_vel_R_nd.z()};
}

namespace {

/**
 * @brief 有効ポテンシャル V(q) の勾配 grad_V を計算
 *
 * @details
 * この関数はCR3BPのハミルトニアン
 *    H_A = T(p) + V(q)
 * のうち有効ポテンシャル項 V(q) の勾配 (偏導関数ベクトル) を計算
 *
 * ここで定義される有効ポテンシャル V(q) は
 *    V(q) = (1-mu)/r1 + mu/r2
 * (注意: これはCR3BPのヤコビ積分 C_J = 2*U - (vx^2+vy^2+vz^2) の
 * U = (1/2)(x^2+y^2) + V(q) とは異なる、シンプレクティック積分用に
 * 分離されたハミルトニアンの項である。)
 *
 * r1 と r2 は、それぞれ第三体 (x,y,z) から主天体 m1 (-mu, 0, 0) と
 * m2 (1-mu, 0, 0) までの距離
 *    r1 = sqrt( (x+mu)^2 + y^2 + z^2 )
 *    r2 = sqrt( (x-(1-mu))^2 + y^2 + z^2 )
 *
 * 計算される勾配 grad_V = (dV/dx, dV/dy, dV/dz) は
 * ApplyKick ステップで運動量を更新するために使用される
 *
 *
 * @param[in]  params   AstroConstants構造体 (gm_sun, gm_earth から mu を計算)
 * @param[in]  x        第三体の x 座標。
 * @param[in]  y        第三体の y 座標。
 * @param[in]  z        第三体の z 座標。
 * @param[out] grad_V   計算された勾配ベクトル (dV/dx, dV/dy, dV/dz) を格納する配列
 *
 * @cite Szebehely1967 "Theory of Orbits" (Chapter 9, The Potential U)
 */
template <typename ScalarType>
void CalculateGradientP(const ScalarType mu, ScalarType x, ScalarType y, ScalarType z,
                        ScalarType* grad_P) {
  const ScalarType mu1 = 1.0 - mu;

  const ScalarType x1 = x + mu;
  const ScalarType x2 = x - (1.0 - mu);

  const ScalarType r1_sq = x1 * x1 + y * y + z * z;
  const ScalarType r2_sq = x2 * x2 + y * y + z * z;

  if (r1_sq < kMinDistanceSq || r2_sq < kMinDistanceSq) {
    throw std::runtime_error("Position too close to primary in CalculateGradientP.");
  }

  const ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
  const ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

  // grad_P = dP/dq = d/dq ((1-mu)/r1 + mu/r2)
  grad_P[0] = -mu1 * x1 * r1_inv3 - mu * x2 * r2_inv3;
  grad_P[1] = -mu1 * y * r1_inv3 - mu * y * r2_inv3;
  grad_P[2] = -mu1 * z * r1_inv3 - mu * z * r2_inv3;
}

/**
 * @brief ポテンシャル V(q) の勾配 grad_V とヘッセ行列 H_V を計算
 *
 * @details
 * H_A = T(p) + V(q)
 * V(q) = -(1-mu)/r1 - mu/r2
 * * Kick (主軌道): p_dot = -grad_V(q)
 * Kick (偏差): dp_dot = -(H_V(q)) * dq
 *
 * この関数は grad_V = dV/dq と H_V = d^2 V / (dq_i dq_j) を計算する
 * (注: 以前の CalculateGradientV は grad(-V) を計算していた)
 *
 * @param[in]  params   AstroConstants構造体
 * @param[in]  x, y, z  主軌道の位置
 * @param[out] grad_V   計算された勾配ベクトル (dV/dx, dV/dy, dV/dz)
 * @param[out] hessian_V 計算されたヘッセ行列 (3x3)
 */
/**
 * @brief P(q) = (1-mu)/r1 + mu/r2 の勾配 grad_P とヘッセ行列 H_P を計算
 * @note V(q) = -P(q)
 * p_dot = +grad_P(q)
 * dp_dot = +H_P(q) * dq
 */
template <typename ScalarType>
void CalculateGradientAndHessianP(const ScalarType mu, ScalarType x, ScalarType y, ScalarType z,
                                  ScalarType* grad_P, HessianMatrix<ScalarType>* hessian_P) {
  const ScalarType mu1 = 1.0 - mu;

  const ScalarType x1 = x + mu;
  const ScalarType x2 = x - (1.0 - mu);

  const ScalarType r1_sq = x1 * x1 + y * y + z * z;
  const ScalarType r2_sq = x2 * x2 + y * y + z * z;

  if (r1_sq < kMinDistanceSq || r2_sq < kMinDistanceSq) {
    throw std::runtime_error("Position too close to primary in CalculateGradientAndHessianP.");
  }

  const ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
  const ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));
  const ScalarType r1_inv5 = r1_inv3 / r1_sq;
  const ScalarType r2_inv5 = r2_inv3 / r2_sq;

  // 1. 勾配 grad_P = dP/dq
  grad_P[0] = -mu1 * x1 * r1_inv3 - mu * x2 * r2_inv3;
  grad_P[1] = -mu1 * y * r1_inv3 - mu * y * r2_inv3;
  grad_P[2] = -mu1 * z * r1_inv3 - mu * z * r2_inv3;

  // 2. ヘッセ行列 H_P = d^2 P / (dq_i dq_j)
  // d/dx(-x1/r1^3) = -(r1^2 - 3*x1^2) / r1^5 = (3x1^2 - r1^2) / r1^5
  // (重大な符号エラーを修正: (r1_inv3 - ...) -> (3*... - r1_inv3))
  ScalarType term1, term2;

  // Hxx = d/dx(grad_P[0])
  term1 = mu1 * (3.0 * x1 * x1 * r1_inv5 - r1_inv3);
  term2 = mu * (3.0 * x2 * x2 * r2_inv5 - r2_inv3);
  hessian_P->hxx = term1 + term2;

  // Hyy = d/dy(grad_P[1])
  term1 = mu1 * (3.0 * y * y * r1_inv5 - r1_inv3);
  term2 = mu * (3.0 * y * y * r2_inv5 - r2_inv3);
  hessian_P->hyy = term1 + term2;

  // Hzz = d/dz(grad_P[2])
  term1 = mu1 * (3.0 * z * z * r1_inv5 - r1_inv3);
  term2 = mu * (3.0 * z * z * r2_inv5 - r2_inv3);
  hessian_P->hzz = term1 + term2;

  // Hxy = d/dy(grad_P[0])
  term1 = mu1 * (3.0 * x1 * y * r1_inv5);
  term2 = mu * (3.0 * x2 * y * r2_inv5);
  hessian_P->hxy = term1 + term2;

  // Hxz = d/dz(grad_P[0])
  term1 = mu1 * (3.0 * x1 * z * r1_inv5);
  term2 = mu * (3.0 * x2 * z * r2_inv5);
  hessian_P->hxz = term1 + term2;

  // Hyz = d/dz(grad_P[1])
  term1 = mu1 * (3.0 * y * z * r1_inv5);
  term2 = mu * (3.0 * y * z * r2_inv5);
  hessian_P->hyz = term1 + term2;
}

/**
 * @brief ハミルトニアン H_B (コリオリ項) による流れを適用
 *
 * @details
 * 円制限三体問題 (CR3BP) のハミルトニアン
 *    H =  (1/2)(px^2 + py^2 + pz^2) -(1-mu)/r1 - mu/r2 - (p_x * y - p_y * x)
 * は
 *　   T(p) = (1/2)(px^2 + py^2 + pz^2)
 *    V(q) = (1-mu)/r1 + mu/r2
 *    H_A = T(p) + V(q)
 *    H_B = - (p_x * y - p_y * x)
 * として
 *    H = H_A + H_B に分離できる
 * H_Bは、コリオリ項に由来するジャイロ項
 * この H_B のみから導出されるハミルトン方程式は以下
 *    q_dot = dH_B/dp  => (x_dot, y_dot) = (-y, x)
 *    p_dot = -dH_B/dq => (px_dot, py_dot) = (-py, px)
 * (z と pz は変化しない)
 *
 * これらの方程式の解は、(x, y) 平面と (px, py) 平面における、
 * 角速度 n=1 の単純な反時計回りの剛体回転運動みたいなかんじ
 *
 * したがって、H_B の流れを時間 `angle` だけ積分する操作は、
 * (x, y) と (px, py) のペアを `angle` ラジアンだけ回転させる
 * 解析解（回転行列の適用）と等価
 *
 * @param[in,out] state 更新対象の正準状態 (q, p)。
 * @param[in]     angle 回転角 (ラジアン)。
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II, II.5 Splitting Methods)
 * @cite yoshida,ハミルトン力学系のためのシンプレクティック数値積分法,
 */
template <typename ScalarType>
void ApplyRotation(CanonicalState<ScalarType>* state, ScalarType angle) {
  const ScalarType c = std::cos(angle);
  const ScalarType s = std::sin(angle);

  const ScalarType qx_old = state->qx;
  const ScalarType qy_old = state->qy;
  const ScalarType px_old = state->px;
  const ScalarType py_old = state->py;

  state->qx = qx_old * c - qy_old * s;
  state->qy = qx_old * s + qy_old * c;

  state->px = px_old * c - py_old * s;
  state->py = px_old * s + py_old * c;
}

/**
 * @brief H_A の Kick (運動量更新) ステップを適用
 *
 * @details
 * これはシンプレクティック積分 (Strangスプリッティング) における
 * リープフロッグ法の「Kick」ステップ
 * ハミルトニアン H_A = T(p) + V(q) ののうち
 * ポテンシャル V(q) による正準運動量 p の変化を計算
 *
 * ハミルトン方程式 p_dot = -dH_A/dq = -grad_V(q) に基づき、
 * 時間 `dt` だけ運動量を更新
 *    p(t + dt) = p(t) + dt * p_dot(t)
 *              = p(t) - dt * grad_V(q(t))
 *
 * この関数は `CalculateGradientV` で計算された勾配 `grad_V` を
 * 減算することでこの力積を適用
 *
 * @param[in]  params 質量比 mu
 * @param[in,out] state 更新対象の正準状態 (q, p)
 * @param[in]  dt     時間ステップ幅
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II.3, Leapfrog)
 */
template <typename ScalarType>
void ApplyKick(const ScalarType mu, CanonicalState<ScalarType>* state, ScalarType dt) {
  ScalarType grad_P[3];
  CalculateGradientP(mu, state->qx, state->qy, state->qz, grad_P);

  // p_dot = +grad_P  =>  p += dt * grad_P
  state->px += dt * grad_P[0];
  state->py += dt * grad_P[1];
  state->pz += dt * grad_P[2];
}

/**
 * @brief H_A の Drift (位置更新) ステップを適用
 *
 * @details
 * これはシンプレクティック積分 (Strangスプリッティング) における
 * リープフロッグ法の「Drift」ステップ
 * ハミルトニアン H_A = T(p) + V(q) の流れのうち
 * 運動エネルギー T(p) = (1/2)(px^2 + py^2 + pz^2) による
 * 正準座標 q の変化を計算
 *
 * ハミルトン方程式 q_dot = dH_A/dp = dT/dp = p に基づき
 * 時間 `dt` だけ位置を更新
 *    q(t + dt) = q(t) + dt * q_dot(t)
 *              = q(t) + dt * p(t)
 *
 * @note
 * このステップでの `q_dot = p` は、あくまで H_A の流れにおける関係式であり
 * 物理的な速度 `v = q_dot_phys` は `v_x = p_x - y` とは違う
 *
 * @param[in,out] state 更新対象の正準状態 (q, p)
 * @param[in]  dt     時間ステップ幅
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II.3, Leapfrog)
 */
template <typename ScalarType>
void ApplyDrift(CanonicalState<ScalarType>* state, ScalarType dt) {
  state->qx += dt * state->px;
  state->qy += dt * state->py;
  state->qz += dt * state->pz;
}
/**
 * @brief H_B (コリオリ項) の流れを SaliState 全体に適用
 */
template <typename ScalarType>
void ApplyRotationSALI(SaliState<ScalarType>* state, ScalarType angle) {
  ApplyRotation(&(state->state), angle);
  ApplyRotation(&(state->w1), angle);
  ApplyRotation(&(state->w2), angle);
}

/**
 * @brief H_A (Drift) の流れを SaliState 全体に適用
 */
template <typename ScalarType>
void ApplyDriftSALI(SaliState<ScalarType>* state, ScalarType dt) {
  ApplyDrift(&(state->state), dt);
  ApplyDrift(&(state->w1), dt);
  ApplyDrift(&(state->w2), dt);
}
/**
 * @brief H_A (Kick) の流れを SaliState 全体に適用
 */
template <typename ScalarType>
void ApplyKickSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType dt) {
  ScalarType grad_P[3];
  HessianMatrix<ScalarType> hessian_P;

  // 主軌道の位置 (q) で 勾配とヘッセ行列を計算
  CalculateGradientAndHessianP(mu, state->state.qx, state->state.qy, state->state.qz, grad_P,
                               &hessian_P);

  // 1. 主軌道 p の更新
  // p_dot = +grad_P  =>  p += dt * grad_P
  state->state.px += dt * grad_P[0];
  state->state.py += dt * grad_P[1];
  state->state.pz += dt * grad_P[2];

  // 2. 偏差ベクトル dp の更新
  // dp_dot = +H_P * dq  =>  dp += dt * (H_P * dq)
  // (ロジックの符号を修正)
  const ScalarType dq1x = state->w1.qx, dq1y = state->w1.qy, dq1z = state->w1.qz;
  state->w1.px += dt * (hessian_P.hxx * dq1x + hessian_P.hxy * dq1y + hessian_P.hxz * dq1z);
  state->w1.py += dt * (hessian_P.hxy * dq1x + hessian_P.hyy * dq1y + hessian_P.hyz * dq1z);
  state->w1.pz += dt * (hessian_P.hxz * dq1x + hessian_P.hyz * dq1y + hessian_P.hzz * dq1z);

  const ScalarType dq2x = state->w2.qx, dq2y = state->w2.qy, dq2z = state->w2.qz;
  state->w2.px += dt * (hessian_P.hxx * dq2x + hessian_P.hxy * dq2y + hessian_P.hxz * dq2z);
  state->w2.py += dt * (hessian_P.hxy * dq2x + hessian_P.hyy * dq2y + hessian_P.hyz * dq2z);
  state->w2.pz += dt * (hessian_P.hxz * dq2x + hessian_P.hyz * dq2y + hessian_P.hzz * dq2z);
}

}  // namespace

/**
 * @brief CR3BPの運動方程式 (修正版)
 *
 * 物理状態 State (q, v) を受け取り
 * その時間微分 d(State)/dt = (v, a) を State 型で返す
 */
template <typename ScalarType>
class EquationOfMotion {
 private:
  ScalarType mu_;

 public:
  explicit EquationOfMotion(const AstroConstants<ScalarType>& params) {
    mu_ = params.gm_earth / (params.gm_earth + params.gm_sun);
  }

  State<ScalarType> operator()(const State<ScalarType>& state, const ScalarType /* t */) const {
    const ScalarType mu1 = 1.0 - mu_;
    const ScalarType x = state.x, y = state.y, z = state.z;

    const ScalarType r1 = calc_r1(x, y, z, mu_);
    const ScalarType r2 = calc_r2(x, y, z, mu_);

    if (r1 == 0.0 || r2 == 0.0) {
      throw std::runtime_error("Error: Position coincides with a primary.");
    }

    const ScalarType r1_inv3 = 1.0 / (r1 * r1 * r1);
    const ScalarType r2_inv3 = 1.0 / (r2 * r2 * r2);

    State<ScalarType> dxdt;
    dxdt.x = state.vx;
    dxdt.y = state.vy;
    dxdt.z = state.vz;
    dxdt.vx = 2.0 * state.vy + x - mu1 * (x + mu_) * r1_inv3 - mu_ * (x - 1.0 + mu_) * r2_inv3;
    dxdt.vy = -2.0 * state.vx + y - mu1 * y * r1_inv3 - mu_ * y * r2_inv3;
    dxdt.vz = -mu1 * z * r1_inv3 - mu_ * z * r2_inv3;
    // std::cout << std::setprecision(15)
    //           << "2 * vy, x, -(1-mu)*(x+mu)/r1^3, -mu*(x-1+mu)/r2^3: " << 2.0 * state.vy << ", "
    //           << x << ", " << -(mu_) << ", " << -(mu_ * (x - 1.0 + mu_) * r2_inv3) << std::endl;
    return dxdt;
  }
};

/**
 * @brief 4次のルンゲ＝クッタ法 (RK4) の1ステップ
 *
 * @param eom 運動方程式 (State を返し State を受け取る)
 * @param state 現在の状態ベクトル
 * @param t 現在の時刻
 * @param h 1ステップの時間幅
 * @return 積分後の状態ベクトル (t + h)
 */
template <typename ScalarType, typename EomType>
State<ScalarType> RungeKutta4Step(const EomType& eom, const State<ScalarType>& state, ScalarType t,
                                  ScalarType h) {
  const State<ScalarType> k1 = eom(state, t);
  // std::cout << "state: " << state.x << ", " << state.y << ", " << state.z << ", " << state.vx
  //           << ", " << state.vy << ", " << state.vz << std::endl;
  // std::cout << "k1: " << k1.x << ", " << k1.y << ", " << k1.z << ", " << k1.vx << ", " << k1.vy
  //           << ", " << k1.vz << std::endl;
  const State<ScalarType> k2 = eom(state + k1 * (h / 2.0), t + h / 2.0);
  const State<ScalarType> k3 = eom(state + k2 * (h / 2.0), t + h / 2.0);
  const State<ScalarType> k4 = eom(state + k3 * h, t + h);
  // std::cout << "k2: " << k2.x << ", " << k2.y << ", " << k2.z << ", " << k2.vx << ", " << k2.vy
  //           << ", " << k2.vz << std::endl;
  // std::cout << "k3: " << k3.x << ", " << k3.y << ", " << k3.z << ", " << k3.vx << ", " << k3.vy
  //           << ", " << k3.vz << std::endl;
  // std::cout << "k4: " << k4.x << ", " << k4.y << ", " << k4.z << ", " << k4.vx << ", " << k4.vy
  //           << ", " << k4.vz << std::endl;
  return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
}

// -----------------------------------------------------------------------------
// 積分器 (シンプレクティック)
// -----------------------------------------------------------------------------

/**
 * @brief 2次のシンプレクティック積分ステップ (Strangスプリッティング)
 *
 * @details
 * この関数は Strang Splitting と呼ばれる手法に基づく
 * 円制限三体問題 (CR3BP) の軌道を1ステップ (時間 h) 積分
 * この手法は O(h^2) の2次精度をもつ
 *
 * 計算は
 * 円制限三体問題 (CR3BP) のハミルトニアン
 *    H =  (1/2)(px^2 + py^2 + pz^2) -(1-mu)/r1 - mu/r2 - (p_x * y - p_y * x)
 * が
 *　   T(p) = (1/2)(px^2 + py^2 + pz^2)
 *    V(q) = (1-mu)/r1 + mu/r2
 *    H_A = T(p) + V(q)
 *    H_B = - (p_x * y - p_y * x)
 * として
 *    H = H_A + H_B に分離できることに基づく
 *
 * 1.B(h/2) - 回転:
 * ハミルトニアンの H_B (コリオリ項) の流れを h/2 時間進めます
 * これは解析的に解ける回転写像として実装
 *
 * 2.A(h) - リープフロッグ (Kick-Drift-Kick)
 * 次に H_A (運動エネルギー + ポテンシャル) の流れを h 時間進める
 * - Kick(h/2): ポテンシャル項 V(q) による運動量の変化 (力積) を h/2 時間適用
 * - Drift(h):  運動エネルギー項 T(p) による位置の変化を h 時間適用
 * - Kick(h/2): 再度 ポテンシャルによる運動量の変化を h/2 時間適用
 *
 * 3.B(h/2) - 回転
 * 最後に 再び H_B (コリオリ項) の流れを h/2 時間進める
 *
 * @param state 積分開始時の物理状態
 * @param h 1ステップの時間幅
 * @return 積分後の物理状態 (t + h)
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep(const ScalarType mu, const State<ScalarType>& state,
                                 ScalarType h) {
  CanonicalState<ScalarType> canonical_state = ConvertToCanonical(state);

  // 1. B(h/2) - 回転
  ApplyRotation(&canonical_state, h / 2.0);

  // 2. A(h) = Kick(h/2) * Drift(h) * Kick(h/2)
  ApplyKick(mu, &canonical_state, h / 2.0);
  ApplyDrift(&canonical_state, h);
  ApplyKick(mu, &canonical_state, h / 2.0);

  // 3. B(h/2) - 回転
  ApplyRotation(&canonical_state, h / 2.0);

  return ConvertToPhysical(canonical_state);
}

/**
 * @brief 4次のシンプレクティック積分ステップ (吉田法, 1990)
 *
 * @details
 * この関数は 吉田法 (1990) を用いて積分精度を4次に高める
 *
 * 計算は S4(tau) = S2(x1 * tau) * S2(x0 * tau) * S2(x1 * tau) という構成
 *
 * 1.
 * 4次精度を達成するために吉田によって導出された
 * 特殊な係数 kX1 と kX0 を定義します
 *
 * 2.1回目の積分
 * `SymplecticStep` を使い kX1 * tau 時間だけ進めます
 *
 * 3.2回目の積分
 * 1の結果を初期値として `SymplecticStep` で kX0 * tau 時間進める
 *
 * 4.3回目の積分
 * 2の結果を初期値として `SymplecticStep` で kX1 * tau 時間進める
 *
 * @param params 質量比 mu を含むパラメータ
 * @param state 積分開始時の物理状態
 * @param tau 1ステップの時間幅
 * @return 積分後の物理状態 (t + tau)
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep4thOrder(const ScalarType mu, const State<ScalarType>& state,
                                         ScalarType tau) {
  // 吉田 (1990) による4次積分のための係数
  // x0 = -2^(1/3) / (2 - 2^(1/3))
  // x1 = 1 / (2 - 2^(1/3))
  const ScalarType kX1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
  const ScalarType kX0 = 1.0 - 2.0 * kX1;  // (x0 + 2*x1 = 1 の関係から)

  // S4(tau) = S2(x1 * tau) * S2(x0 * tau) * S2(x1 * tau)
  State<ScalarType> state1 = SymplecticStep(mu, state, kX1 * tau);
  State<ScalarType> state2 = SymplecticStep(mu, state1, kX0 * tau);
  State<ScalarType> state3 = SymplecticStep(mu, state2, kX1 * tau);

  return state3;
}

// -----------------------------------------------------------------------------
// 積分器 (SALI対応シンプレクティック)
// -----------------------------------------------------------------------------

/**
 * @brief 2次のSALI対応シンプレクティック積分ステップ (Strangスプリッティング)
 *
 * @param[in]     params 質量比 mu を含むパラメータ
 * @param[in,out] state  積分対象の SaliState
 * @param[in]     h      1ステップの時間幅
 */
template <typename ScalarType>
void SymplecticStepSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType h) {
  ApplyRotationSALI(state, h / 2.0);
  ApplyKickSALI(mu, state, h / 2.0);
  ApplyDriftSALI(state, h);
  ApplyKickSALI(mu, state, h / 2.0);
  ApplyRotationSALI(state, h / 2.0);
}

/**
 * @brief 4次のSALI対応シンプレクティック積分ステップ (吉田法)
 *
 * @param[in]     params 質量比 mu を含むパラメータ
 * @param[in,out] state  積分対象の SaliState
 * @param[in]     tau    1ステップの時間幅
 */
template <typename ScalarType>
void SymplecticStep4thOrderSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType tau) {
  // 吉田 (1990) による4次積分のための係数
  const ScalarType kX1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
  const ScalarType kX0 = 1.0 - 2.0 * kX1;

  // S4(tau) = S2(x1 * tau) * S2(x0 * tau) * S2(x1 * tau)
  SymplecticStepSALI(mu, state, kX1 * tau);
  SymplecticStepSALI(mu, state, kX0 * tau);
  SymplecticStepSALI(mu, state, kX1 * tau);
}
/**
 * @brief 汎用積分ドライバ
 *
 * @tparam StateType 状態ベクトルの型 (例: State<double>)
 * @tparam IntegratorType 1ステップ積分を実行する "Callable" (関数, ラムダ, Functor)
 * シグネチャ: StateType(const StateType&, ScalarType, ScalarType)
 * @tparam ObserverType 各ステップで呼び出される "Callable" (関数, ラムダ, Functor)
 * シグネチャ: void(const StateType&, ScalarType)
 */
template <typename StateType, typename IntegratorType, typename ObserverType, typename ScalarType>
void Integrate(StateType& current_state, IntegratorType integrator, ObserverType observer,
               ScalarType start_time, ScalarType time_step, int num_steps) {
  ScalarType current_time = start_time;
  observer(current_state, current_time);
  for (int i = 0; i < num_steps; i++) {
    // std::cout << "state before integration: " << current_state.x << ", " << current_state.y << ",
    // "
    //           << current_state.z << ", " << current_state.vx << ", " << current_state.vy << ", "
    //           << current_state.vz << std::endl;
    // // ちょっとディレイ
    // std::this_thread::sleep_for(std::chrono::milliseconds(100));
    // std::cout << "integrator = " << &integrator << std::endl;
    // std::cout << "i = " << i << std::endl;
    current_state = integrator(current_state, current_time, time_step);
    current_time += time_step;
    observer(current_state, current_time);
  }
}
/**
 * @brief SALI計算用の汎用積分ドライバ
 *
 * @details
 * 積分ステップの実行 (Integrator)
 * SALIの計算 (Observer)
 * 偏差ベクトルの正規化 (Renormalization)
 * の3つを管理します
 */
template <typename IntegratorType, typename ObserverType, typename ScalarType>
void IntegrateSALI(SaliState<ScalarType>& current_state, IntegratorType integrator,
                   ObserverType observer, ScalarType start_time, ScalarType time_step,
                   int num_steps) {
  ScalarType current_time = start_time;
  observer(current_state, current_time);

  for (int i = 0; i < num_steps; ++i) {
    integrator(&current_state, time_step);
    current_time += time_step;

    // 偏差ベクトルを正規化
    current_state.w1.Normalize();
    current_state.w2.Normalize();
    // (注: 厳密には Gram-Schmidt法で直交化も推奨

    observer(current_state, current_time);
  }
}

/**
 * @brief コンソールに状態を出力するオブザーバー (例)
 */
template <typename StateType, typename ScalarType>
class ConsoleObserver {
 public:
  void operator()(const StateType& state, ScalarType t) {
    std::cout << std::fixed << std::setprecision(8) << "Time: " << t << "\t x: " << state.x
              << "\t y: " << state.y << "\t z: " << state.z << "\t vx: " << state.vx
              << "\t vy: " << state.vy << "\t vz: " << state.vz << std::endl;
  }
};

template <typename ScalarType>
class StateFileObserver {
 private:
  std::ostream& os_;  // 出力ストリームへの参照

 public:
  /**
   * @brief コンストラクタ
   * @param os [in,out] main側で開かれた書き込み先のファイルストリーム
   */
  explicit StateFileObserver(std::ostream& os) : os_(os) {}

  /**
   * @brief 積分ステップごとに呼び出され SALIと状態をストリームに書き込む
   */
  void operator()(const State<ScalarType>& state, ScalarType t) {
#pragma omp critical

    if (!os_.good()) return;

    // CSV形式で出力 (精度は高めに設定)
    os_ << std::fixed << std::setprecision(15) << t << "," << state.x << "," << state.y << ","
        << state.z << "," << state.vx << "," << state.vy << "," << state.vz << "\n";
  }
};
template <typename ScalarType>

class StateBufferObserver {
 private:
  std::vector<std::array<ScalarType, 8>>& history_;
  ScalarType mu_;

 public:
  explicit StateBufferObserver(std::vector<std::array<ScalarType, 8>>& history, ScalarType mu)
      : history_(history), mu_(mu) {}

  void operator()(const State<ScalarType>& state, ScalarType t) {
    history_.push_back({t, calc_jacobi_integral(state, mu_), state.x, state.y, state.z, state.vx,
                        state.vy, state.vz});
  }

  const std::vector<std::array<ScalarType, 8>>& GetHistory() const { return history_; }
};
template <typename ScalarType>
class SaliObserver {
 private:
  std::vector<ScalarType> sali_history_;
  std::vector<ScalarType> time_history_;

 public:
  void operator()(const SaliState<ScalarType>& state, ScalarType t) {
    // w1 と w2 は正規化されている前提で SALI を計算
    // SALI = min( ||w1+w2||, ||w1-w2|| )
    const ScalarType norm_plus = (state.w1 + state.w2).Norm();
    const ScalarType norm_minus = (state.w1 - state.w2).Norm();
    const ScalarType sali = std::min(norm_plus, norm_minus);

    sali_history_.push_back(sali);
    time_history_.push_back(t);

    // コンソールに出力
    std::cout << "Time: " << std::fixed << std::setprecision(8) << t << "\t SALI: " << sali
              << "\t (x: " << state.state.x << ")" << std::endl;
  }

  // (結果を取得するためのゲッター関数をここに追加可能)
};

template <typename ScalarType>
class SaliFileObserver {
 private:
  std::ostream& os_;  // 出力ストリームへの参照

 public:
  /**
   * @brief コンストラクタ
   * @param os [in,out] main側で開かれた書き込み先のファイルストリーム
   */
  explicit SaliFileObserver(std::ostream& os) : os_(os) {}

  /**
   * @brief 積分ステップごとに呼び出され SALIと状態をストリームに書き込む
   */
  void operator()(const SaliState<ScalarType>& state, ScalarType t) {
#pragma omp critical

    if (!os_.good()) return;

    // w1 と w2 は正規化されている前提で SALI を計算
    const ScalarType norm_plus = (state.w1 + state.w2).Norm();
    const ScalarType norm_minus = (state.w1 - state.w2).Norm();
    const ScalarType sali = std::min(norm_plus, norm_minus);

    // CSV形式で出力 (精度は高めに設定)
    os_ << std::fixed << std::setprecision(15) << t << "," << sali << "," << state.state.qx << ","
        << state.state.qy << "," << state.state.qz << "," << state.state.px << "," << state.state.py
        << "," << state.state.pz << "\n";
  }
};

}  // namespace crtbp
#endif  // RTBP_HPP
