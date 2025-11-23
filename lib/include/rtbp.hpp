/**
 * @file rtbp.hpp
 * @author tabata
 * @brief 円制限三体問題の関数をまとめたライブラリ
 * @version 2.0 (SALI対応シンプレクティック積分器 修正版)
 * @date 2025-11-10
 * @note C++20
 * (conceptを使っているからc++20より昔のバージョンを指定してコンパイルするとエラーになる)
 * @par history
 * - 2025-01-24 tabata version 1.0
 * - 2025-11-06 refactored template parameter names　and modified SALI calc method
 */
#ifndef RTBP_HPP
#define RTBP_HPP
#include <algorithm>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <concepts>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
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

  State<ScalarType> operator+=(const State<ScalarType>& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    vx += other.vx;
    vy += other.vy;
    vz += other.vz;
    return *this;
  }
};

// Allow scalar * State for Boost.Odeint vector_space_algebra
template <typename ScalarType>
State<ScalarType> operator*(ScalarType scalar, const State<ScalarType>& state) {
  return state * scalar;
}
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

  CanonicalState<ScalarType> operator*(ScalarType scalar) const {
    return {qx * scalar, qy * scalar, qz * scalar, px * scalar, py * scalar, pz * scalar};
  }

  CanonicalState<ScalarType>& operator+=(const CanonicalState<ScalarType>& other) {
    qx += other.qx;
    qy += other.qy;
    qz += other.qz;
    px += other.px;
    py += other.py;
    pz += other.pz;
    return *this;
  }

  // 6次元位相空間でのノルム (SALI計算用)
  ScalarType Norm() const {
    return std::sqrt(qx * qx + qy * qy + qz * qz + px * px + py * py + pz * pz);
  }
  ScalarType Dot(const CanonicalState<ScalarType> other) const {
    return (qx * other.qx + qy * other.qy + qz * other.qz + px * other.px + py * other.py +
            pz * other.pz);
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

// Allow scalar * CanonicalState for Boost.Odeint vector_space_algebra
template <typename ScalarType>
CanonicalState<ScalarType> operator*(ScalarType scalar, const CanonicalState<ScalarType>& state) {
  return state * scalar;
}
/**
 * @brief ポテンシャルのヘッセ行列 (2階微分)
 * @details H_ij = d^2 U / (dq_i dq_j)
 *          対称行列なので独立成分のみ保持
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
  * @details 主軌道と2つの偏差ベクトルを含む

 */
template <typename ScalarType>
struct SaliState {
  CanonicalState<ScalarType> state;  // 主軌道の状態 (q, p)
  CanonicalState<ScalarType> w1;     // 偏差ベクトル1 (dq1, dp1)
  CanonicalState<ScalarType> w2;     // 偏差ベクトル2 (dq2, dp2)

  SaliState<ScalarType> operator+(const SaliState<ScalarType>& other) const {
    return {state + other.state, w1 + other.w1, w2 + other.w2};
  }

  SaliState<ScalarType> operator*(ScalarType scalar) const {
    return {state * scalar, w1 * scalar, w2 * scalar};
  }

  SaliState<ScalarType>& operator+=(const SaliState<ScalarType>& other) {
    state += other.state;
    w1 += other.w1;
    w2 += other.w2;
    return *this;
  }
};

// Allow scalar * SaliState for Boost.Odeint vector_space_algebra
template <typename ScalarType>
SaliState<ScalarType> operator*(ScalarType scalar, const SaliState<ScalarType>& state) {
  return state * scalar;
}
};  // namespace my_type

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

/**
 * @brief 円制限三体問題の関数をまとめた名前空間
 *
 */
namespace crtbp {
using namespace param;
using namespace my_type;
/**
 * @brief 円制限三体問題の系で物理状態 (q, q_dot) を正準状態 (q, p) に変換
 * p_x = vx - y
 * p_y = vy + x
 * p_z = vz
 */
template <typename ScalarType>
CanonicalState<ScalarType> ConvertToCanonical(const State<ScalarType>& state) {
  return CanonicalState<ScalarType>{
      state.x,
      state.y,
      state.z,
      state.vx - state.y,  // px = vx - y
      state.vy + state.x,  // py = vy + x
      state.vz             // pz = vz
  };
}

/**
 * @brief 円制限三体問題の系で正準状態 (q, p) を物理状態 (q, q_dot) に変換
 * vx = p_x + y
 * vy = p_y - x
 * vz = p_z
 */
template <typename ScalarType>
State<ScalarType> ConvertToPhysical(const CanonicalState<ScalarType>& canonical_state) {
  return State<ScalarType>{
      canonical_state.qx,
      canonical_state.qy,
      canonical_state.qz,
      canonical_state.px + canonical_state.qy,  // vx = px + y
      canonical_state.py - canonical_state.qx,  // vy = py - x
      canonical_state.pz                        // vz = pz
  };
}
/**
 * @brief 第三体から主天体m1への距離r1を計算
 * @param x, y, z 第三体の位置
 * @param mu 質量比 mu = m2/(m1+m2)
 * @return r1 = ||r - r1||
 */
template <typename ScalarType>
inline ScalarType calc_r1(ScalarType x, ScalarType y, ScalarType z, ScalarType mu) {
  const ScalarType x1 = x + mu;
  return std::sqrt(x1 * x1 + y * y + z * z);
}

/**
 * @brief 第三体から主天体m2への距離r2を計算
 * @param x, y, z 第三体の位置
 * @param mu 質量比 mu = m2/(m1+m2)
 * @return r2 = ||r - r2||
 */
template <typename ScalarType>
inline ScalarType calc_r2(ScalarType x, ScalarType y, ScalarType z, ScalarType mu) {
  const ScalarType x2 = x - (1.0 - mu);
  return std::sqrt(x2 * x2 + y * y + z * z);
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
  // return 0.5 * (x * x + y * y) + (1.0 - mu) / r1 + mu / r2;
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
  ScalarType potential = calc_potential_U(x, y, z, mu);

  // return std::sqrt(x * x + y * y + 2. * (1. - mu) / r1 + 2. * mu / r2 + mu * (1. - mu) -
  //                  JACOBI_INTEGRAL);
  return std::sqrt(2.0 * potential - JACOBI_INTEGRAL);
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

template <typename ScalarType>
const State<ScalarType> ConvertInertial2Rotating____(
    const State<ScalarType>& ast_state, const State<ScalarType>& p2_state,
    const AstroConstants<ScalarType>& astro_params) {
  std::array<Vector3d<ScalarType>, 2> init_ast_state_G;
  std::array<Vector3d<ScalarType>, 2> init_p2_state_G;
  init_ast_state_G = {{Vector3d<ScalarType>(ast_state.x, ast_state.y, ast_state.z),
                       Vector3d<ScalarType>(ast_state.vx, ast_state.vy, ast_state.vz)}};
  init_p2_state_G = {{Vector3d<ScalarType>(p2_state.x, p2_state.y, p2_state.z),
                      Vector3d<ScalarType>(p2_state.vx, p2_state.vy, p2_state.vz)}};
  const ScalarType mu =
      astro_params.gm_earth / (astro_params.gm_sun + astro_params.gm_earth);  // 質量比
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

  //  慣性座標系から見た地球軌道面の法線ベクトルと慣性系から見た中間座標系rot1のz軸がなす角度を計算
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
  Vector3d<ScalarType> ast_pos_R{ast_pos_R2.x() + 1.0 - epos_R2.x() - mu, ast_pos_R2.y(),
                                 ast_pos_R2.z()};

  /* calc the velocity of rotating frame */
  //  変換行列convert_G_to_R1の微分を計算
  // theta1の微分

  ScalarType ND_time_ref =
      std::sqrt(astro_params.au * astro_params.au / astro_params.gm_sun);  // Non-Dimensionaltime
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
  //                      e_G_x * (e_G_x * ND_e_G_vx + e_G_sy * ND_e_G_vy +
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
  return {State<ScalarType>(ast_pos_R.x(), ast_pos_R.y(), ast_pos_R.z(), ast_vel_R.x(),
                            ast_vel_R.y(), ast_vel_R.z())};
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

/*
 * @brief シンプレクティック積分法を用いた計算ステップのロジックを構成するパーツ
 */
namespace {

/**
 * @brief 有効ポテンシャルU(q)の勾配を計算
 *
 * @details
 *          U(q) = -(1-mu)/r1 - mu/r2
 *          grad_U = dU/dq
 *
 *          ハミルトニアンは H_A = T(p) + U(q) の形なので
 *          dH_A/dq = dU/dq = grad_U
 * 計算される勾配 grad_V  は
 * ApplyKick ステップで運動量を更新するために使用される
 *
 * @param[in] mu 質量比
 * @param[in] x, y, z 第三体の位置
 * @param[out] grad_U 出力: 勾配ベクトル [dU/dx, dU/dy, dU/dz]
 */
template <typename ScalarType>
void CalculateGradientU(const ScalarType mu, ScalarType x, ScalarType y, ScalarType z,
                        ScalarType* grad_U) {
  const ScalarType mu1 = 1.0 - mu;

  const ScalarType x1 = x + mu;
  const ScalarType x2 = x - (1.0 - mu);

  const ScalarType r1_sq = x1 * x1 + y * y + z * z;
  const ScalarType r2_sq = x2 * x2 + y * y + z * z;

  if (r1_sq < kMinDistanceSq || r2_sq < kMinDistanceSq) {
    throw std::runtime_error("Position too close to primary in CalculateGradientU.");
  }

  const ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
  const ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

  // grad_U = dU/dq = d/dq ((1-mu)/r1 + mu/r2)
  // dU/dx = (1-mu)*x1/r1^3 + mu*x2/r2^3
  grad_U[0] = -mu1 * (-x1 * r1_inv3) - mu * (-x2 * r2_inv3);
  grad_U[1] = -mu1 * (-y * r1_inv3) - mu * (-y * r2_inv3);
  grad_U[2] = -mu1 * (-z * r1_inv3) - mu * (-z * r2_inv3);
}

/**
 * @brief 有効ポテンシャルU(q)の勾配とヘッセ行列を計算
 *
 * @details U(q) = -(1-mu)/r1 - mu/r2
 *          grad_U = dU/dq
 *          H_U = d^2U/(dq_i dq_j)
 *
 *          主軌道: p_dot = -dH_A/dq = -dU/dq = -grad_U
 *          偏差: dp_dot = H_U * dq
 *
 * @param mu 質量比
 * @param x, y, z 主軌道の位置
 * @param grad_U 出力: 勾配ベクトル
 * @param hessian_U 出力: ヘッセ行列
 */
template <typename ScalarType>
void CalculateGradientAndHessianU(const ScalarType mu, ScalarType x, ScalarType y, ScalarType z,
                                  ScalarType* grad_U, HessianMatrix<ScalarType>* hessian_U) {
  const ScalarType mu1 = 1.0 - mu;

  const ScalarType x1 = x + mu;
  const ScalarType x2 = x - (1.0 - mu);

  const ScalarType r1_sq = x1 * x1 + y * y + z * z;
  const ScalarType r2_sq = x2 * x2 + y * y + z * z;

  if (r1_sq < kMinDistanceSq || r2_sq < kMinDistanceSq) {
    throw std::runtime_error("Position too close to primary in CalculateGradientAndHessianU.");
  }

  const ScalarType r1_inv3 = 1.0 / (r1_sq * std::sqrt(r1_sq));
  const ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));
  const ScalarType r1_inv5 = r1_inv3 / r1_sq;
  const ScalarType r2_inv5 = r2_inv3 / r2_sq;

  // 1. 勾配 grad_U = dU/dq
  grad_U[0] = -mu1 * (-x1 * r1_inv3) - mu * (-x2 * r2_inv3);
  grad_U[1] = -mu1 * (-y * r1_inv3) - mu * (-y * r2_inv3);
  grad_U[2] = -mu1 * (-z * r1_inv3) - mu * (-z * r2_inv3);

  // 2. ヘッセ行列 H_U = d^2U/(dq_i dq_j)
  // d/dx(grad_U[0]) = d/dx((1-mu)*x1/r1^3 + mu*x2/r2^3)
  //                 = (1-mu)*(1/r1^3 - 3*x1^2/r1^5) + mu*(1/r2^3 - 3*x2^2/r2^5)
  ScalarType term1, term2;

  // Hxx
  term1 = mu1 * (r1_inv3 - 3.0 * x1 * x1 * r1_inv5);
  term2 = mu * (r2_inv3 - 3.0 * x2 * x2 * r2_inv5);
  hessian_U->hxx = term1 + term2;

  // Hyy
  term1 = mu1 * (r1_inv3 - 3.0 * y * y * r1_inv5);
  term2 = mu * (r2_inv3 - 3.0 * y * y * r2_inv5);
  hessian_U->hyy = term1 + term2;

  // Hzz
  term1 = mu1 * (r1_inv3 - 3.0 * z * z * r1_inv5);
  term2 = mu * (r2_inv3 - 3.0 * z * z * r2_inv5);
  hessian_U->hzz = term1 + term2;

  // Hxy (対称なのでHyx = Hxy)
  // d/dy(grad_U[0]) = (1-mu)*(-3*x1*y/r1^5) + mu*(-3*x2*y/r2^5)
  term1 = mu1 * (-3.0 * x1 * y * r1_inv5);
  term2 = mu * (-3.0 * x2 * y * r2_inv5);
  hessian_U->hxy = term1 + term2;

  // Hxz
  term1 = mu1 * (-3.0 * x1 * z * r1_inv5);
  term2 = mu * (-3.0 * x2 * z * r2_inv5);
  hessian_U->hxz = term1 + term2;

  // Hyz
  term1 = mu1 * (-3.0 * y * z * r1_inv5);
  term2 = mu * (-3.0 * y * z * r2_inv5);
  hessian_U->hyz = term1 + term2;
}

/**
 * @brief ハミルトニアン H_B (コリオリ項) による流れを適用
 *
 * @details
 * 円制限三体問題 (CR3BP) のハミルトニアン
 *    H =  (1/2)(px^2 + py^2 + pz^2) -(1-mu)/r1 - mu/r2 + (p_x * y - p_y * x)
 * は
 *　   T(p) = (1/2)(px^2 + py^2 + pz^2)
 *    V(q) = -(1-mu)/r1 - mu/r2
 *    H_A = T(p) + V(q)
 *    H_B = (p_x * y - p_y * x)
 * として
 *    H = H_A + H_B に分離できる
 *
 * この H_B のみから導出されるハミルトン方程式は以下
 *    q_dot = dH_B/dp  => (x_dot, y_dot) = (y, -x)
 *    p_dot = -dH_B/dq => (px_dot, py_dot) = (py, -px)
 * (z と pz は変化しない)
 *
 * ここでたとえば
 *    x_dot　= y, y_dot = -x
 * からは
 *    x(t) = Acos(t)+Bsin(t), y(t) = Ccos(t)+Dsin(t),
 * という厳密解を求められるため、ここではシンプレクティック法関係なく厳密解を計算する
 *
 * @param[in,out] state 更新対象の正準状態 (q, p)。
 * @param[in]     angle 回転角 (ラジアン)。
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II, II.5 Splitting Methods)
 * @cite yoshida,ハミルトン力学系のためのシンプレクティック数値積分法,
 */
template <typename ScalarType>
void RotateState(CanonicalState<ScalarType>* state, ScalarType angle) {
  const ScalarType c = std::cos(angle);
  const ScalarType s = std::sin(angle);

  const ScalarType qx_old = state->qx;
  const ScalarType qy_old = state->qy;
  const ScalarType px_old = state->px;
  const ScalarType py_old = state->py;

  state->qx = qx_old * c + qy_old * s;
  state->qy = -qx_old * s + qy_old * c;

  state->px = px_old * c + py_old * s;
  state->py = -px_old * s + py_old * c;
}

/**
 * @brief H_AのKick（運動量更新）ステップを適用
 *
 * @details シンプレクティック積分におけるKickステップ
 *          H_A = T(p) + U(q) のポテンシャル項による運動量変化
 *
 *          dH_A/dq = dU/dq = grad_U
 *          p(t+dt) = p(t) - dt * dH_A/dq = p(t) - dt * grad_U
 *
 * @param[in]  params 質量比 mu
 * @param[in,out] state 更新対象の正準状態 (q, p)
 * @param[in]  dt     時間ステップ幅とシンプレクティック法の係数の積
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II.5, splitting method)
 * @cite Hairer2006 "Geometric Numerical Integration" (p.150 V.3.1 Symmetric Composition of First
 * Order Methods)
 * @cite 吉田 春男, 「ハミルトニアン力学系のためのシンプレクティック数値積分法」,
 * 共同研究「非線形現象の数理科学」湘南レクチャー論文集, p. 68-83, 1997.
 */
template <typename ScalarType>
void UpdateMomentum(const ScalarType mu, CanonicalState<ScalarType>* state, ScalarType dt) {
  ScalarType grad_U[3];
  CalculateGradientU(mu, state->qx, state->qy, state->qz, grad_U);

  state->px -= dt * grad_U[0];
  state->py -= dt * grad_U[1];
  state->pz -= dt * grad_U[2];
}

/**
 * @brief H_AのDrift（位置更新）ステップを適用
 *
 * @details シンプレクティック積分におけるDriftステップ
 *          H_A = T(p) + U(q) の運動エネルギー項の寄与を計算
 *
 *          dH_A/dp = dT/dp = p
 *          q(t+dt) = q(t) + dt * dH_A/dp = q(t) + dt * p
 *
 * @param[in,out] state 更新対象の正準状態 (q, p)
 * @param[in]  dt     時間ステップ幅とシンプレクティック法の係数の積
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II.5, splitting method)
 * @cite Hairer2006 "Geometric Numerical Integration" (p.150 V.3.1 Symmetric Composition of First
 * Order Methods)
 * @cite 吉田 春男, 「ハミルトニアン力学系のためのシンプレクティック数値積分法」,
 * 共同研究「非線形現象の数理科学」湘南レクチャー論文集, p. 68-83, 1997.
 */
template <typename ScalarType>
void UpdatePosition(CanonicalState<ScalarType>* state, ScalarType dt) {
  state->qx += dt * state->px;
  state->qy += dt * state->py;
  state->qz += dt * state->pz;
}

/**
 * @brief H_B（コリオリ項）の流れをSaliState全体に適用
 * @param[in,out] state SALI計算用の状態
 * @param[in] angle 回転角（ラジアン）
 */
template <typename ScalarType>
void RotateStateSALI(SaliState<ScalarType>* state, ScalarType angle) {
  RotateState(&(state->state), angle);
  RotateState(&(state->w1), angle);
  RotateState(&(state->w2), angle);
}

/**
 * @brief H_A（Drift）の流れをSaliState全体に適用
 * @param[in,out] state SALI計算用の状態
 * @param[in] dt 時間ステップ幅
 */
template <typename ScalarType>
void UpdatePositionSALI(SaliState<ScalarType>* state, ScalarType dt) {
  UpdatePosition(&(state->state), dt);
  UpdatePosition(&(state->w1), dt);
  UpdatePosition(&(state->w2), dt);
}

/**
 * @brief H_A（Kick）の流れをSaliState全体に適用
 *
 * @details 主軌道の位置でポテンシャルの勾配とヘッセ行列を計算し、
 *          主軌道と偏差ベクトルの両方を更新
 *
 *          主軌道: p_dot = -dH_A/dq = -dU/dq = -grad_U
 *          偏差: dp_dot = -H_U * dq
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II.5, splitting method)
 * @cite Hairer2006 "Geometric Numerical Integration" (p.150 V.3.1 Symmetric Composition of First
 * Order Methods)
 * @cite 吉田 春男, 「ハミルトニアン力学系のためのシンプレクティック数値積分法」,
 * 共同研究「非線形現象の数理科学」湘南レクチャー論文集, p. 68-83, 1997.
 *
 * @param mu 質量比
 * @param state SALI計算用の状態
 * @param dt 時間ステップ幅
 */
template <typename ScalarType>
void UpdateMomentumSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType dt) {
  ScalarType grad_U[3];
  HessianMatrix<ScalarType> hessian_U;

  // 主軌道の位置でポテンシャルの勾配とヘッセ行列を計算
  CalculateGradientAndHessianU(mu, state->state.qx, state->state.qy, state->state.qz, grad_U,
                               &hessian_U);

  // 1. 主軌道の運動量更新
  // p_dot = -grad_U  =>  p -= dt * grad_U
  state->state.px -= dt * grad_U[0];
  state->state.py -= dt * grad_U[1];
  state->state.pz -= dt * grad_U[2];

  // 2. 偏差ベクトルの運動量更新
  // dp_dot = -H_U * dq  =>  dp -= dt * (H_U * dq)

  // 偏差ベクトル1
  const ScalarType dq1x = state->w1.qx;
  const ScalarType dq1y = state->w1.qy;
  const ScalarType dq1z = state->w1.qz;

  state->w1.px -= dt * (hessian_U.hxx * dq1x + hessian_U.hxy * dq1y + hessian_U.hxz * dq1z);
  state->w1.py -= dt * (hessian_U.hxy * dq1x + hessian_U.hyy * dq1y + hessian_U.hyz * dq1z);
  state->w1.pz -= dt * (hessian_U.hxz * dq1x + hessian_U.hyz * dq1y + hessian_U.hzz * dq1z);

  // 偏差ベクトル2
  const ScalarType dq2x = state->w2.qx;
  const ScalarType dq2y = state->w2.qy;
  const ScalarType dq2z = state->w2.qz;

  state->w2.px -= dt * (hessian_U.hxx * dq2x + hessian_U.hxy * dq2y + hessian_U.hxz * dq2z);
  state->w2.py -= dt * (hessian_U.hxy * dq2x + hessian_U.hyy * dq2y + hessian_U.hyz * dq2z);
  state->w2.pz -= dt * (hessian_U.hxz * dq2x + hessian_U.hyz * dq2y + hessian_U.hzz * dq2z);
}
/**
 * @brief Gram-Schmidt法で2つの偏差ベクトルを直交化
 *
 * @details SALI計算では偏差ベクトルが直交していることが重要
 *          w1はそのまま、w2からw1成分を除去して直交化
 *
 *          w2_orth = w2 - (w2·w1) * w1
 *
 * @param state SALI計算用の状態
 */
template <typename ScalarType>
void OrthogonalizeDeviationVectors(SaliState<ScalarType>* state) {
  // w1はそのまま正規化
  state->w1.Normalize();

  // w2からw1方向の成分を除去
  ScalarType projection = state->w2.Dot(state->w1);

  state->w2.qx -= projection * state->w1.qx;
  state->w2.qy -= projection * state->w1.qy;
  state->w2.qz -= projection * state->w1.qz;
  state->w2.px -= projection * state->w1.px;
  state->w2.py -= projection * state->w1.py;
  state->w2.pz -= projection * state->w1.pz;

  // w2も正規化
  state->w2.Normalize();
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
    // std::cout << "dxdt.x: " << dxdt.x << ", dxdt.y: " << dxdt.y << ", dxdt.z: " << dxdt.z
    //           << ", dxdt.vx: " << dxdt.vx << ", dxdt.vy: " << dxdt.vy << ", dxdt.vz: " << dxdt.vz
    //           << std::endl;
    return dxdt;
  }
};

/**
 * @brief System wrapper for odeint Runge-Kutta-Dopri5 (orbit only).
 */
template <typename ScalarType>
class Dopri5OrbitSystem {
 private:
  EquationOfMotion<ScalarType> eom_;

 public:
  explicit Dopri5OrbitSystem(const AstroConstants<ScalarType>& params) : eom_(params) {}

  void operator()(const State<ScalarType>& state, State<ScalarType>& dxdt,
                  const ScalarType t) const {
    dxdt = eom_(state, t);
  }
};

/**
 * @brief System wrapper that advances the orbit and two deviation vectors for SALI with Dopri5.
 */
template <typename ScalarType>
class Dopri5SaliSystem {
 private:
  ScalarType mu_;

  static void ComputeVariationDerivative(const CanonicalState<ScalarType>& w,
                                         const HessianMatrix<ScalarType>& hessian,
                                         CanonicalState<ScalarType>* dw_dt) {
    const ScalarType dq_x = w.qx;
    const ScalarType dq_y = w.qy;
    const ScalarType dq_z = w.qz;

    // dq/dt = rotation(q) + dp
    dw_dt->qx = w.px + w.qy;
    dw_dt->qy = w.py - w.qx;
    dw_dt->qz = w.pz;

    // dp/dt = rotation(p) - Hessian(U) * dq
    dw_dt->px = w.py - (hessian.hxx * dq_x + hessian.hxy * dq_y + hessian.hxz * dq_z);
    dw_dt->py = -w.px - (hessian.hxy * dq_x + hessian.hyy * dq_y + hessian.hyz * dq_z);
    dw_dt->pz = -(hessian.hxz * dq_x + hessian.hyz * dq_y + hessian.hzz * dq_z);
  }

 public:
  explicit Dopri5SaliSystem(const AstroConstants<ScalarType>& params) {
    mu_ = params.gm_earth / (params.gm_earth + params.gm_sun);
  }

  void operator()(const SaliState<ScalarType>& state, SaliState<ScalarType>& dstate_dt,
                  const ScalarType /*t*/) const {
    ScalarType grad_U[3];
    HessianMatrix<ScalarType> hessian;
    CalculateGradientAndHessianU(mu_, state.state.qx, state.state.qy, state.state.qz, grad_U,
                                 &hessian);

    // Central orbit
    dstate_dt.state.qx = state.state.px + state.state.qy;
    dstate_dt.state.qy = state.state.py - state.state.qx;
    dstate_dt.state.qz = state.state.pz;

    dstate_dt.state.px = state.state.py - grad_U[0];
    dstate_dt.state.py = -state.state.px - grad_U[1];
    dstate_dt.state.pz = -grad_U[2];

    // Variational equations for two deviation vectors
    ComputeVariationDerivative(state.w1, hessian, &dstate_dt.w1);
    ComputeVariationDerivative(state.w2, hessian, &dstate_dt.w2);
  }
};
// -----------積分器----------------------------------------------------------

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
  const State<ScalarType> k2 = eom(state + k1 * (h / 2.0), t + h / 2.0);
  const State<ScalarType> k3 = eom(state + k2 * (h / 2.0), t + h / 2.0);
  const State<ScalarType> k4 = eom(state + k3 * h, t + h);
  // std::cout << "k1.vx: " << k1.vx << ", k1.vy: " << k1.vy << ", k1.vz: " << k1.vz << std::endl;
  // std::cout << "k2.vx: " << k2.vx << ", k2.vy: " << k2.vy << ", k2.vz: " << k2.vz << std::endl;
  // std::cout << "k3.vx: " << k3.vx << ", k3.vy: " << k3.vy << ", k3.vz: " << k3.vz << std::endl;
  // std::cout << "k4.vx: " << k4.vx << ", k4.vy: " << k4.vy << ", k4.vz: " << k4.vz << std::endl;
  // std::cout << "returning state.vx: "
  //           << state.vx + (k1.vx + 2.0 * k2.vx + 2.0 * k3.vx + k4.vx) * (h / 6.0)
  //           << ", state.vy: " << state.vy + (k1.vy + 2.0 * k2.vy + 2.0 * k3.vy + k4.vy) * (h
  //           / 6.0)
  //           << ", state.vz: " << state.vz + (k1.vz + 2.0 * k2.vz + 2.0 * k3.vz + k4.vz) * (h
  //           / 6.0)
  //           << std::endl;
  return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
}

/**
 * @brief 2次のシンプレクティック積分ステップ（Strangスプリッティング）
 *
 * @details H = H_A + H_B に対する2次精度のシンプレクティック積分
 *
 *          分離方法:
 *          H_A = T(p) - U(q) = (1/2)(px^2+py^2+pz^2) - ((1-mu)/r1 + mu/r2)
 *          H_B = -(px*y - py*x) (コリオリ項)
 *
 *          積分スキーム: B(h/2) * A(h) * B(h/2)
 *          ここでA(h) = K(h/2) * D(h) * K(h/2) (Leapfrog)
 *
 *          完全な手順:
 *          1. B(h/2) - 回転 h/2
 *          2. K(h/2) - Kick h/2
 *          3. D(h)   - Drift h
 *          4. K(h/2) - Kick h/2
 *          5. B(h/2) - 回転 h/2
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II.5, splitting method)
 * @cite Hairer2006 "Geometric Numerical Integration" (p.150 V.3.1 Symmetric Composition of First
 * Order Methods)
 * @cite 吉田 春男, 「ハミルトニアン力学系のためのシンプレクティック数値積分法」,
 * 共同研究「非線形現象の数理科学」湘南レクチャー論文集, p. 68-83, 1997.
 *
 * @param mu 質量比
 * @param state 積分開始時の物理状態
 * @param h 1ステップの時間幅
 * @return 積分後の物理状態（t+h）
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep(const ScalarType mu, const State<ScalarType>& state,
                                 ScalarType h) {
  CanonicalState<ScalarType> canonical_state = ConvertToCanonical(state);

  // 1. B(h/2) - 回転
  RotateState(&canonical_state, h / 2.0);
  // 2. A(h) = Kick(h/2) * Drift(h) * Kick(h/2)
  // qを先に更新　　更新したqを使って最終的なp,qを出す
  UpdatePosition(&canonical_state, h / 2.0);
  UpdateMomentum(mu, &canonical_state, h);
  UpdatePosition(&canonical_state, h / 2.0);

  // 3. B(h/2) - 回転
  RotateState(&canonical_state, h / 2.0);
  return ConvertToPhysical(canonical_state);
}

/**
 * @brief 4次のシンプレクティック積分ステップ（吉田法、1990）
 *
 * @details 吉田の4次シンプレクティック積分法を使用
 *
 *          S4(tau) = S2(x1*tau) * S2(x0*tau) * S2(x1*tau)
 *
 *          係数:
 *          x1 = 1 / (2 - 2^(1/3))
 *          x0 = -2^(1/3) / (2 - 2^(1/3)) = 1 - 2*x1
 *
 *          これにより4次精度（誤差 O(h^5)）を実現
 *
 * @param[in]  mu 質量比
 * @param[in]  state 積分開始時の物理状態
 * @param[in]  tau 1ステップの時間幅
 * @return 積分後の物理状態（t+tau）
 *
 * @cite Hairer2006 "Geometric Numerical Integration" (Chapter II.5, splitting method)
 * @cite Hairer2006 "Geometric Numerical Integration" (p.150 V.3.1 Symmetric Composition of First
 * Order Methods)
 * @cite 吉田 春男, 「ハミルトニアン力学系のためのシンプレクティック数値積分法」,
 * 共同研究「非線形現象の数理科学」湘南レクチャー論文集, p. 68-83, 1997.
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep4thOrder(const ScalarType mu, const State<ScalarType>& state,
                                         ScalarType tau) {
  // 吉田（1990）による4次積分の係数
  const ScalarType kX1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
  const ScalarType kX0 = 1.0 - 2.0 * kX1;

  // S4(tau) = S2(x1*tau) * S2(x0*tau) * S2(x1*tau)
  State<ScalarType> state1 = SymplecticStep(mu, state, kX1 * tau);
  State<ScalarType> state2 = SymplecticStep(mu, state1, kX0 * tau);
  State<ScalarType> state3 = SymplecticStep(mu, state2, kX1 * tau);

  return state3;
}

/**
 * @brief 2次のSALI対応シンプレクティック積分ステップ
 *
 * @details 主軌道と2つの偏差ベクトルを同時に積分
 *          Strangスプリッティング B(h/2)*A(h)*B(h/2) を適用
 * @param[in] mu 質量比
 * @param[in,out] state 積分対象のSaliState
 * @param[in] h 1ステップの時間幅
 */
template <typename ScalarType>
void SymplecticStepSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType h) {
  RotateStateSALI(state, h / 2.0);
  // 先にqを更新
  UpdatePositionSALI(state, h / 2.0);
  UpdateMomentumSALI(mu, state, h);
  UpdatePositionSALI(state, h / 2.0);

  RotateStateSALI(state, h / 2.0);
}

/**
 * @brief 4次のSALI対応シンプレクティック積分ステップ（吉田法）
 *
 * @details 吉田の4次法をSALI計算に適用
 *          S4(tau) = S2(x1*tau) * S2(x0*tau) * S2(x1*tau)
 *
 * @param[in]  mu 質量比
 * @param[in,out] state 積分対象のSaliState
 * @param[in]  tau 1ステップの時間幅
 */
template <typename ScalarType>
void SymplecticStep4thOrderSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType tau) {
  const ScalarType kX1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
  const ScalarType kX0 = 1.0 - 2.0 * kX1;

  SymplecticStepSALI(mu, state, kX1 * tau);
  SymplecticStepSALI(mu, state, kX0 * tau);
  SymplecticStepSALI(mu, state, kX1 * tau);
}

// -----------積分器のラッパー関数----------------------------------------------------------

/**
 * @brief 汎用積分ドライバ
 *
 * @details 指定されたIntegratorとObserverを使用して
 *          num_stepsステップの積分を実行
 *
 * @tparam StateType 状態ベクトルの型（例: State<double>）
 * @tparam IntegratorType 1ステップ積分を実行する関数オブジェクト
 *         シグネチャ: StateType(const StateType&, ScalarType, ScalarType)
 * @tparam ObserverType 各ステップで呼び出される関数オブジェクト
 *         シグネチャ: void(const StateType&, ScalarType)
 * @tparam ScalarType スカラー型（通常はdouble）
 *
 * @param[in,out] current_state　初期状態
 * （更新されるため注意！！！計算の各ステップを保存したければobserverを用意すること）
 * @param[in] integrator 積分関数
 * @param[in] observer 観測関数
 * @param[in] start_time 積分開始時刻
 * @param[in] time_step 時間ステップ幅
 * @param[in] num_steps ステップ数
 */
template <typename IntegratorType, typename StateType, typename ScalarType>
concept Integrator =
    requires(IntegratorType func, const StateType& state, ScalarType t, ScalarType dt) {
      { func(state, t, dt) } -> std::convertible_to<StateType>;
    };

template <typename ObserverType, typename StateType, typename ScalarType>
concept Observer = requires(ObserverType func, const StateType& state, ScalarType t) {
  { func(state, t) } -> std::same_as<void>;
};

template <typename StateType, typename IntegratorType, typename ObserverType, typename ScalarType>
  requires Integrator<IntegratorType, StateType, ScalarType> &&
           Observer<ObserverType, StateType, ScalarType>
void Integrate(StateType& initial_state, IntegratorType integrator, ObserverType observer,
               ScalarType start_time, ScalarType time_step, int num_steps) {
  ScalarType current_time = start_time;
  StateType state_after_step = initial_state;
  observer(initial_state, current_time);
  for (int i = 0; i < num_steps; ++i) {
    // std::cout << "Integrate step " << i + 1 << "/" << num_steps << std::endl;
    // std::cout << "Current time: " << current_time << std::endl;
    // std::cout << "State before step: x=" << state_after_step.x << ", y=" << state_after_step.y
    //           << ", z=" << state_after_step.z << ", vx=" << state_after_step.vx
    //           << ", vy=" << state_after_step.vy << ", vz=" << state_after_step.vz << std::endl;
    state_after_step = integrator(state_after_step, current_time, time_step);
    current_time += time_step;
    observer(state_after_step, current_time);
  }
}
/**
 * @brief SALI計算用の汎用積分ドライバ
 *
 * @details 積分ステップの実行、SALIの計算、偏差ベクトルの
 *          正規化と直交化を管理
 *
 *          各ステップ後に:
 *          1. 偏差ベクトルをGram-Schmidt法で直交化
 *          2. 両ベクトルを正規化
 *          3. Observer関数を呼び出してSALIを計算
 *
 * @tparam IntegratorType void(*)(SaliState<ScalarType>*, ScalarType)
 * @tparam ObserverType void(*)(const SaliState<ScalarType>&, ScalarType)
 * @tparam ScalarType スカラー型（通常はdouble）
 *
 * @param[in,out] current_state 初期SALI状態（更新される）
 * @param[in] integrator 積分関数
 * @param[in] observer 観測関数（SALI計算を含む）
 * @param[in] start_time 積分開始時刻
 * @param[in] time_step 時間ステップ幅
 * @param[in] num_steps ステップ数
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

    // 偏差ベクトルを直交化・正規化
    // （SALI計算の精度向上のため重要）
    OrthogonalizeDeviationVectors(&current_state);

    observer(current_state, current_time);
  }
}

/**
 * @brief Fixed-step Dopri5 integration for the orbit only.
 */
template <typename ScalarType, typename ObserverType,
          typename Stepper = boost::numeric::odeint::runge_kutta_dopri5<
              State<ScalarType>, ScalarType, State<ScalarType>, ScalarType,
              boost::numeric::odeint::vector_space_algebra>>
  requires Observer<ObserverType, State<ScalarType>, ScalarType>
void IntegrateDopri5Orbit(const AstroConstants<ScalarType>& params, State<ScalarType>& state,
                          ScalarType start_time, ScalarType end_time, ScalarType time_step,
                          ObserverType observer, Stepper stepper = Stepper()) {
  Dopri5OrbitSystem<ScalarType> system(params);

  ScalarType current_time = start_time;
  observer(state, current_time);
  while (current_time < end_time) {
    const ScalarType dt = std::min(time_step, end_time - current_time);
    stepper.do_step(system, state, current_time, dt);
    current_time += dt;
    observer(state, current_time);
  }
}

/**
 * @brief Fixed-step Dopri5 integration for orbit + variational equations (SALI).
 *
 * @param orthonormalize_deviation true to orthogonalize w1/w2 after every step.
 */
template <typename ScalarType, typename ObserverType,
          typename Stepper = boost::numeric::odeint::runge_kutta_dopri5<
              SaliState<ScalarType>, ScalarType, SaliState<ScalarType>, ScalarType,
              boost::numeric::odeint::vector_space_algebra>>
  requires Observer<ObserverType, SaliState<ScalarType>, ScalarType>
void IntegrateDopri5SALI(const AstroConstants<ScalarType>& params, SaliState<ScalarType>& state,
                         ScalarType start_time, ScalarType end_time, ScalarType time_step,
                         ObserverType observer, bool orthonormalize_deviation = true,
                         Stepper stepper = Stepper()) {
  Dopri5SaliSystem<ScalarType> system(params);

  ScalarType current_time = start_time;
  observer(state, current_time);
  while (current_time < end_time) {
    const ScalarType dt = std::min(time_step, end_time - current_time);
    stepper.do_step(system, state, current_time, dt);
    current_time += dt;
    if (orthonormalize_deviation) {
      OrthogonalizeDeviationVectors(&state);
    }
    observer(state, current_time);
  }
}

/**
 * @brief integrate関数に渡すためのオブザーバー　コンソールに状態を出力する
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

/**
 * @brief integrate関数に渡すためのオブザーバー　各ステップを指定ファイルに出力する
 */
template <typename ScalarType>
class StateFileObserver {
 private:
  std::ostream& os_;
  // 出力ストリームへの参照

 public:
  /**
   * @brief コンストラクタ
   * @param os [in,out] main側で開かれた書き込み先のファイルストリーム
   * @note 呼び出し側で std::ofstream ofs("filename");って宣言したofsをコンストラクタに渡す
   */
  explicit StateFileObserver(std::ostream& os) : os_(os) {}

  /**
   * @brief integrate関数に渡すためのオブザーバー　各ステップを指定ファイルに出力する
   */
  void operator()(const State<ScalarType>& state, ScalarType t) {
#pragma omp critical

    if (!os_.good()) return;

    os_ << std::fixed << std::setprecision(15) << t << "," << state.x << "," << state.y << ","
        << state.z << "," << state.vx << "," << state.vy << "," << state.vz << "\n";
  }
};

/**
 * @brief
 * integrate関数に渡すためのオブザーバー　各ステップをバッファにコピーしつつついてにヤコビ積分も計算する
 */
template <typename ScalarType>
class StateBufferObserver {
 private:
  std::vector<std::array<ScalarType, 8>>& history_;
  ScalarType mu_;
  ScalarType JacobiIntegral_;

 public:
  explicit StateBufferObserver(std::vector<std::array<ScalarType, 8>>& history, ScalarType mu,
                               ScalarType JacobiIntegral)
      : history_(history), mu_(mu), JacobiIntegral_(JacobiIntegral) {}

  void operator()(const State<ScalarType>& state, ScalarType t) {
    // std::cout << "x: " << state.x << ", y: " << state.y << ", z: " << state.z
    //           << ", vx: " << state.vx << ", vy: " << state.vy << ", vz: " << state.vz <<
    //           std::endl;
    history_.push_back({t, JacobiIntegral_ - calc_jacobi_integral(state, mu_), state.x, state.y,
                        state.z, state.vx, state.vy, state.vz});
  }

  const std::vector<std::array<ScalarType, 8>>& GetHistory() const { return history_; }
};

/**
 * @brief GALI (Generalized Alignment Index) を偏差ベクトル集合から計算するユーティリティ.
 *
 * @details 入力ベクトルは自動的に正規化され、修正グラム・シュミット法で
 *          k体並進体積を求める. GALI_k = ||\hat{w}_1 \wedge ... \wedge \hat{w}_k||
 *          の定義と等価であり、値が急減するかどうかで軌道のカオス性を評価できる.
 *
 * @param deviation_vectors 正規化前でもよい偏差ベクトル群 (2 <= k <= 6 を想定)
 * @param zero_threshold 正規直交化の途中で線形独立性を失ったと見做す閾値
 * @return ScalarType GALI_k の値 (k = deviation_vectors.size())
 *
 * @throws std::invalid_argument 入力が空、または次元上限を超える場合に送出
 */
template <typename ScalarType>
ScalarType ComputeGALI(std::span<const CanonicalState<ScalarType>> deviation_vectors,
                       ScalarType zero_threshold = static_cast<ScalarType>(1e-12)) {
  const std::size_t k = deviation_vectors.size();
  if (k == 0) {
    throw std::invalid_argument("ComputeGALI requires at least one deviation vector.");
  }
  constexpr std::size_t kPhaseSpaceDim = 6;
  if (k > kPhaseSpaceDim) {
    throw std::invalid_argument("ComputeGALI supports up to 6 deviation vectors in CRTBP.");
  }
  if (zero_threshold <= static_cast<ScalarType>(0)) {
    throw std::invalid_argument("ComputeGALI expects a strictly positive threshold.");
  }

  std::vector<CanonicalState<ScalarType>> orthonormal_basis;
  orthonormal_basis.reserve(k);

  ScalarType gali = static_cast<ScalarType>(1);
  for (const auto& deviation : deviation_vectors) {
    CanonicalState<ScalarType> v = deviation;
    v.Normalize();

    for (const auto& basis_vec : orthonormal_basis) {
      const ScalarType projection = v.Dot(basis_vec);
      v.qx -= projection * basis_vec.qx;
      v.qy -= projection * basis_vec.qy;
      v.qz -= projection * basis_vec.qz;
      v.px -= projection * basis_vec.px;
      v.py -= projection * basis_vec.py;
      v.pz -= projection * basis_vec.pz;
    }

    const ScalarType orth_norm = v.Norm();
    if (orth_norm < zero_threshold) {
      return static_cast<ScalarType>(0);
    }

    gali *= orth_norm;

    const ScalarType inv_norm = static_cast<ScalarType>(1) / orth_norm;
    v.qx *= inv_norm;
    v.qy *= inv_norm;
    v.qz *= inv_norm;
    v.px *= inv_norm;
    v.py *= inv_norm;
    v.pz *= inv_norm;

    orthonormal_basis.push_back(v);
  }

  return gali;
}

template <typename ScalarType>
ScalarType ComputeGALI(const std::vector<CanonicalState<ScalarType>>& deviation_vectors,
                       ScalarType zero_threshold = static_cast<ScalarType>(1e-12)) {
  return ComputeGALI<ScalarType>(
      std::span<const CanonicalState<ScalarType>>(deviation_vectors.data(), deviation_vectors.size()),
      zero_threshold);
}

template <typename ScalarType, std::size_t N>
ScalarType ComputeGALI(const std::array<CanonicalState<ScalarType>, N>& deviation_vectors,
                       ScalarType zero_threshold = static_cast<ScalarType>(1e-12)) {
  return ComputeGALI<ScalarType>(
      std::span<const CanonicalState<ScalarType>>(deviation_vectors.data(), deviation_vectors.size()),
      zero_threshold);
}

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
// ---------------------------------------------------------------------------
// Additional high-order symplectic schemes (state and SALI variants)
// ---------------------------------------------------------------------------

/**
 * @brief Blanes-Moan 4th-order symmetric (ABAH4/SS4) with reduced error constant.
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep4thOrderBM(const ScalarType mu, const State<ScalarType>& state,
                                           ScalarType h) {
  constexpr ScalarType a1 = static_cast<ScalarType>(0.0792036964311957);
  constexpr ScalarType a2 = static_cast<ScalarType>(0.353172906049774);
  constexpr ScalarType a3 = static_cast<ScalarType>(-0.0420650803577195);
  constexpr ScalarType a4 = static_cast<ScalarType>(1.0 - 2.0 * (a1 + a2 + a3));
  constexpr ScalarType b1 = static_cast<ScalarType>(0.209515106613362);
  constexpr ScalarType b2 = static_cast<ScalarType>(-0.143851773179818);
  constexpr ScalarType b3 = static_cast<ScalarType>(0.434336666566456);

  CanonicalState<ScalarType> cs = ConvertToCanonical(state);
  auto A = [&](ScalarType dt) {
    ApplyDrift(&cs, dt / 2.0);
    ApplyKick(mu, &cs, dt);
    ApplyDrift(&cs, dt / 2.0);
  };
  auto B = [&](ScalarType dt) { ApplyRotation(&cs, dt); };

  A(a1 * h);
  B(b1 * h);
  A(a2 * h);
  B(b2 * h);
  A(a3 * h);
  B(b3 * h);
  A(a4 * h);
  B(b3 * h);
  A(a3 * h);
  B(b2 * h);
  A(a2 * h);
  B(b1 * h);
  A(a1 * h);

  return ConvertToPhysical(cs);
}

/**
 * @brief 6th-order Yoshida (1990) using 7-fold composition of the 2nd-order Strang step.
 */
template <typename ScalarType>
State<ScalarType> SymplecticStep6thOrder(const ScalarType mu, const State<ScalarType>& state,
                                         ScalarType h) {
  constexpr ScalarType w1 = static_cast<ScalarType>(0.784513610477560);
  constexpr ScalarType w2 = static_cast<ScalarType>(0.235573213359357);
  constexpr ScalarType w3 = static_cast<ScalarType>(-1.17767998417887);
  constexpr ScalarType w4 = static_cast<ScalarType>(1.31518632068391);
  constexpr ScalarType weights[7] = {w1, w2, w3, w4, w3, w2, w1};

  State<ScalarType> s = state;
  for (auto c : weights) {
    s = SymplecticStep(mu, s, c * h);
  }
  return s;
}

/**
 * @brief Blanes-Moan 4th order for SALI (same coefficients as SymplecticStep4thOrderBM).
 */
template <typename ScalarType>
void SymplecticStep4thOrderBMSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType h) {
  constexpr ScalarType a1 = static_cast<ScalarType>(0.0792036964311957);
  constexpr ScalarType a2 = static_cast<ScalarType>(0.353172906049774);
  constexpr ScalarType a3 = static_cast<ScalarType>(-0.0420650803577195);
  constexpr ScalarType a4 = static_cast<ScalarType>(1.0 - 2.0 * (a1 + a2 + a3));
  constexpr ScalarType b1 = static_cast<ScalarType>(0.209515106613362);
  constexpr ScalarType b2 = static_cast<ScalarType>(-0.143851773179818);
  constexpr ScalarType b3 = static_cast<ScalarType>(0.434336666566456);

  auto A = [&](ScalarType dt) {
    ApplyDriftSALI(state, dt / 2.0);
    ApplyKickSALI(mu, state, dt);
    ApplyDriftSALI(state, dt / 2.0);
  };
  auto B = [&](ScalarType dt) { ApplyRotationSALI(state, dt); };

  A(a1 * h);
  B(b1 * h);
  A(a2 * h);
  B(b2 * h);
  A(a3 * h);
  B(b3 * h);
  A(a4 * h);
  B(b3 * h);
  A(a3 * h);
  B(b2 * h);
  A(a2 * h);
  B(b1 * h);
  A(a1 * h);
}

/**
 * @brief 6th-order Yoshida for SALI (7-fold composition of 2nd-order step).
 */
template <typename ScalarType>
void SymplecticStep6thOrderSALI(const ScalarType mu, SaliState<ScalarType>* state, ScalarType h) {
  constexpr ScalarType w1 = static_cast<ScalarType>(0.784513610477560);
  constexpr ScalarType w2 = static_cast<ScalarType>(0.235573213359357);
  constexpr ScalarType w3 = static_cast<ScalarType>(-1.17767998417887);
  constexpr ScalarType w4 = static_cast<ScalarType>(1.31518632068391);
  constexpr ScalarType weights[7] = {w1, w2, w3, w4, w3, w2, w1};

  for (auto c : weights) {
    SymplecticStepSALI(mu, state, c * h);
  }
}

}  // namespace crtbp

/**
 * @brief Helper namespace that groups both circular (crtbp) and elliptic (ertbp) RTBP utilities.
 *
 * The existing circular toolkit is exposed as rtbp::crtbp (namespace alias). New elliptic helpers
 * live under rtbp::ertbp.
 */
namespace rtbp {
namespace crtbp = ::crtbp;

namespace ertbp {
using my_type::State;
using my_type::State3d;

/**
 * @brief Elliptic two-body elements used to place the primaries.
 *
 * All quantities are non-dimensional: the semi-major axis and mean motion default to 1, matching
 * the scaling of the circular toolkit.
 */
template <typename ScalarType>
struct Elements {
  ScalarType mu;                   ///< Mass parameter m2/(m1+m2)
  ScalarType eccentricity;         ///< Orbital eccentricity (0 <= e < 1)
  ScalarType mean_motion = static_cast<ScalarType>(1);    ///< Mean motion (rad / time unit)
  ScalarType semi_major_axis = static_cast<ScalarType>(1);  ///< Relative semi-major axis of the primaries

  Elements() = default;
  Elements(ScalarType mu_in, ScalarType e_in, ScalarType n_in = static_cast<ScalarType>(1),
           ScalarType a_in = static_cast<ScalarType>(1))
      : mu(mu_in), eccentricity(e_in), mean_motion(n_in), semi_major_axis(a_in) {
    if (eccentricity < 0.0 || eccentricity >= 1.0) {
      throw std::invalid_argument("eccentricity must satisfy 0 <= e < 1 for elliptic RTBP.");
    }
  }
};

/**
 * @brief Ephemeris of the two primaries in the inertial barycentric frame.
 */
template <typename ScalarType>
struct PrimariesState {
  Vector3d<ScalarType> primary1_position;  ///< m1 (larger body) position
  Vector3d<ScalarType> primary2_position;  ///< m2 (smaller body) position
  Vector3d<ScalarType> primary1_velocity;
  Vector3d<ScalarType> primary2_velocity;
  ScalarType separation = 0.0;  ///< |r2 - r1|
};

/**
 * @brief Solve Kepler's equation M = E - e sin(E) for the eccentric anomaly.
 */
template <typename ScalarType>
ScalarType SolveKeplerEquation(ScalarType mean_anomaly, ScalarType eccentricity,
                               int max_iterations = 30,
                               ScalarType tolerance = static_cast<ScalarType>(1e-12)) {
  ScalarType E = mean_anomaly;
  for (int i = 0; i < max_iterations; ++i) {
    const ScalarType f = E - eccentricity * std::sin(E) - mean_anomaly;
    const ScalarType f_prime = static_cast<ScalarType>(1) - eccentricity * std::cos(E);
    if (std::abs(f_prime) < tolerance) {
      break;
    }
    const ScalarType delta = -f / f_prime;
    E += delta;
    if (std::abs(delta) < tolerance) {
      break;
    }
  }
  return E;
}

/**
 * @brief Compute barycentric positions/velocities of the primaries on an ellipse in the xy-plane.
 *
 * The relative motion of the primaries follows the standard two-body ellipse with semi-major axis
 * `a` and eccentricity `e`. The returned coordinates place the system barycenter at the origin.
 */
template <typename ScalarType>
PrimariesState<ScalarType> ComputePrimariesEphemeris(const Elements<ScalarType>& elements,
                                                     ScalarType t) {
  const ScalarType mean_anomaly = elements.mean_motion * t;
  const ScalarType E = SolveKeplerEquation(mean_anomaly, elements.eccentricity);
  const ScalarType cos_E = std::cos(E);
  const ScalarType sin_E = std::sin(E);

  const ScalarType one_minus_e2 =
      static_cast<ScalarType>(1) - elements.eccentricity * elements.eccentricity;
  const ScalarType sqrt_one_minus_e2 = std::sqrt(std::max(static_cast<ScalarType>(0), one_minus_e2));
  const ScalarType denom = static_cast<ScalarType>(1) - elements.eccentricity * cos_E;

  // Relative m1->m2 state in perifocal (pericenter aligned with +x).
  const ScalarType rel_x = elements.semi_major_axis * (cos_E - elements.eccentricity);
  const ScalarType rel_y = elements.semi_major_axis * sqrt_one_minus_e2 * sin_E;
  const ScalarType rel_vx =
      -elements.semi_major_axis * elements.mean_motion * sin_E / denom;
  const ScalarType rel_vy =
      elements.semi_major_axis * elements.mean_motion * sqrt_one_minus_e2 * cos_E / denom;

  PrimariesState<ScalarType> ephem;
  ephem.primary1_position = Vector3d<ScalarType>(-elements.mu * rel_x, -elements.mu * rel_y, 0);
  ephem.primary2_position =
      Vector3d<ScalarType>((static_cast<ScalarType>(1) - elements.mu) * rel_x,
                           (static_cast<ScalarType>(1) - elements.mu) * rel_y, 0);
  ephem.primary1_velocity = Vector3d<ScalarType>(-elements.mu * rel_vx, -elements.mu * rel_vy, 0);
  ephem.primary2_velocity =
      Vector3d<ScalarType>((static_cast<ScalarType>(1) - elements.mu) * rel_vx,
                           (static_cast<ScalarType>(1) - elements.mu) * rel_vy, 0);
  ephem.separation = elements.semi_major_axis * denom;

  return ephem;
}

template <typename ScalarType>
inline ScalarType calc_r1(const State<ScalarType>& state,
                          const PrimariesState<ScalarType>& primaries) {
  const Vector3d<ScalarType> delta(state.x - primaries.primary1_position.x(),
                                   state.y - primaries.primary1_position.y(),
                                   state.z - primaries.primary1_position.z());
  return delta.magnitude();
}

template <typename ScalarType>
inline ScalarType calc_r2(const State<ScalarType>& state,
                          const PrimariesState<ScalarType>& primaries) {
  const Vector3d<ScalarType> delta(state.x - primaries.primary2_position.x(),
                                   state.y - primaries.primary2_position.y(),
                                   state.z - primaries.primary2_position.z());
  return delta.magnitude();
}

template <typename ScalarType>
ScalarType calc_potential_U(const State<ScalarType>& state,
                            const PrimariesState<ScalarType>& primaries,
                            ScalarType mu) {
  const ScalarType r1 = calc_r1(state, primaries);
  const ScalarType r2 = calc_r2(state, primaries);
  if (r1 == 0.0 || r2 == 0.0) {
    throw std::runtime_error("Position coincides with a primary in calc_potential_U (ERTBP).");
  }
  return -(static_cast<ScalarType>(1) - mu) / r1 - mu / r2;
}

template <typename ScalarType>
ScalarType calc_mechanical_energy(const State<ScalarType>& state,
                                  const PrimariesState<ScalarType>& primaries,
                                  ScalarType mu) {
  const ScalarType v_sq = state.vx * state.vx + state.vy * state.vy + state.vz * state.vz;
  return static_cast<ScalarType>(0.5) * v_sq + calc_potential_U(state, primaries, mu);
}

/**
 * @brief Elliptic restricted three-body equations in the inertial barycentric frame.
 *
 * The primaries evolve on an ellipse; the third body feels their instantaneous Newtonian gravity.
 */
template <typename ScalarType>
class EquationOfMotion {
 private:
  Elements<ScalarType> elements_;

 public:
  explicit EquationOfMotion(const Elements<ScalarType>& elements) : elements_(elements) {}

  State<ScalarType> operator()(const State<ScalarType>& state, const ScalarType t) const {
    const PrimariesState<ScalarType> primaries = ComputePrimariesEphemeris(elements_, t);
    const Vector3d<ScalarType> r(state.x, state.y, state.z);
    const Vector3d<ScalarType> r1 = r - primaries.primary1_position;
    const Vector3d<ScalarType> r2 = r - primaries.primary2_position;

    const ScalarType dist1 = r1.magnitude();
    const ScalarType dist2 = r2.magnitude();
    if (dist1 == 0.0 || dist2 == 0.0) {
      throw std::runtime_error("Position coincides with a primary in EquationOfMotion (ERTBP).");
    }

    const ScalarType inv_r1_3 = static_cast<ScalarType>(1) / (dist1 * dist1 * dist1);
    const ScalarType inv_r2_3 = static_cast<ScalarType>(1) / (dist2 * dist2 * dist2);

    const ScalarType mu1 = static_cast<ScalarType>(1) - elements_.mu;
    const Vector3d<ScalarType> accel =
        r1 * (-mu1 * inv_r1_3) + r2 * (-elements_.mu * inv_r2_3);

    State<ScalarType> dxdt;
    dxdt.x = state.vx;
    dxdt.y = state.vy;
    dxdt.z = state.vz;
    dxdt.vx = accel.x();
    dxdt.vy = accel.y();
    dxdt.vz = accel.z();
    return dxdt;
  }
};

/**
 * @brief System wrapper for odeint Runge-Kutta-Dopri5 (orbit only) in ERTBP.
 */
template <typename ScalarType>
class Dopri5OrbitSystem {
 private:
  EquationOfMotion<ScalarType> eom_;

 public:
  explicit Dopri5OrbitSystem(const Elements<ScalarType>& elements) : eom_(elements) {}

  void operator()(const State<ScalarType>& state, State<ScalarType>& dxdt,
                  const ScalarType t) const {
    dxdt = eom_(state, t);
  }
};

/**
 * @brief Fixed-step Dopri5 integration helper for the elliptic RTBP.
 */
template <typename ScalarType, typename ObserverType,
          typename Stepper = boost::numeric::odeint::runge_kutta_dopri5<
              State<ScalarType>, ScalarType, State<ScalarType>, ScalarType,
              boost::numeric::odeint::vector_space_algebra>>
void IntegrateDopri5Orbit(const Elements<ScalarType>& elements, State<ScalarType>& state,
                          ScalarType start_time, ScalarType end_time, ScalarType time_step,
                          ObserverType observer, Stepper stepper = Stepper()) {
  Dopri5OrbitSystem<ScalarType> system(elements);

  ScalarType current_time = start_time;
  observer(state, current_time);
  while (current_time < end_time) {
    const ScalarType dt = std::min(time_step, end_time - current_time);
    stepper.do_step(system, state, current_time, dt);
    current_time += dt;
    observer(state, current_time);
  }
}

// Reuse the generic RK4/Integrate helpers from the circular namespace for convenience.
using ::crtbp::Integrate;
using ::crtbp::RungeKutta4Step;

}  // namespace ertbp
}  // namespace rtbp
#endif  // RTBP_HPP
