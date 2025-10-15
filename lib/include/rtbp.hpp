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

#include "vector3d.hpp"
#include <array>
#include <vector>
inline static constexpr double G = 6.67430e-11;
inline static constexpr double AU = 1.49597870700e11;
inline static constexpr double M_SUN = 1.989400e30;
inline static constexpr double M_EARTH = 5.97219e24;
inline static constexpr double DEFAULT_MU = 3.003e-6;

// 6D状態を表す型のエイリアス
using state = std::array<double, 6>;
// 3D座標を表す型のエイリアス
using Coord3D = std::array<double, 3>;

enum class SolverType {
  DOPRI5,
  RK4,
};

// EquationOfMotionのインターフェースクラス
class EquationOfMotionInterface {
public:
  // virtual ~EquationOfMotionInterface() = default; // 仮想デストラクタ
  virtual void operator()(const state &x, state &dxdt,
                          const double t) const = 0; // 純粋仮想関数
};
namespace rtbp {

/**
 * @brief Create a frame conversion matrix object
 * @param unit_n 規格化されたか回転軸
 * @param theta
 * @return std::array<std::array<double, 3>, 3>
 */
static std::array<std::array<double, 3>, 3>
create_frame_conversion_matrix(const Vector3d &unit_n, double theta);

/**
 * @brief Create a coordinate rotation matrix object
 *
 * @param unit_n 規格化されたか回転軸
 * @param theta
 * @return std::array<std::array<double, 3>, 3>
 */
static std::array<std::array<double, 3>, 3>
create_coordinate_rotation_matrix(const Vector3d &unit_n, double theta);

/**
 * @brief Create a Rodrigues rotation matrix object
 *
 * @param unit_n 規格化された回転軸
 * @param theta
 * @return std::array<std::array<double, 3>, 3>
 */
static std::array<std::array<double, 3>, 3>
create_rodrigues_matrix(const Vector3d &unit_n, double theta);

/**
 * @brief Apply a rotation matrix to a vector *
 * @param rot_matrix matrix to apply
 * @param v vector to rotate
 * @return Vector3d
 */
static Vector3d apply_matrix(std::array<std::array<double, 3>, 3> rot_matrix,
                             const Vector3d &v);

/**
 * @brief 円制限三体問題の関数をまとめた名前空間
 *
 */
namespace crtbp {

template <typename T, std::size_t N>
const std::array<T, N> frame_transformation(const std::array<T, N> &ast_state,
                                            const std::array<T, N> &p2_state,
                                            double mu);

template <typename T>
double calc_v_abs(const T &r, const double JACOBI_INTEGRAL, const double mu);
// CRTBPの運動方程式
// EquationOfMotionを継承
class EquationOfMotion : public EquationOfMotionInterface {
private:
  double mu_;

public:
  explicit EquationOfMotion(double mu);
  void operator()(const state &x, state &dxdt, const double /* t */) const;
};
class observer {
private:
  double init_jacobi_;
  double mu_;

public:
  std::vector<std::array<double, 8>> &m_out;
  explicit observer(double init_jacobi, double mu,
                    std::vector<std::array<double, 8>> &out);
  void operator()(const state &x, double t) const;
};
class Trajectory {
private:
  double jacobi_integral_;
  double mu_;

public:
  std::vector<std::array<double, 8>> trajectory_;

  Trajectory(double jacobi_integral, double mu);
  std::vector<std::array<double, 8>> integrate(const state &x0, double t0,
                                               double tf, double dt);
};

/**
 * @brief 円制限三体問題の運動方程式を積分する関数
 * @param x0 初期状態 [x, y, z, vx, vy, vz]
 * @param t0 初期時刻
 * @param tf 終了時刻
 * @param dt 時間刻み
 * @param mu 質量パラメータ（デフォルト: DEFAULT_MU）
 * @param solver_type ソルバータイプ（デフォルト: DOPRI5）
 * @return std::vector<std::array<double, 8>> 軌道データ [t, x, y, z, vx, vy,
 * vz, jacobi_error]
 */
std::vector<std::array<double, 8>>
integrate(const state &x0, double t0, double tf, double dt, double mu,
          SolverType solver_type = SolverType::DOPRI5);

/**
 * @brief 適応的ステップサイズで円制限三体問題の運動方程式を積分する関数
 * @param x0 初期状態 [x, y, z, vx, vy, vz]
 * @param t0 初期時刻
 * @param tf 終了時刻
 * @param dt_init 初期時間刻み
 * @param abs_err 絶対誤差許容値
 * @param rel_err 相対誤差許容値
 * @param mu 質量パラメータ（デフォルト: DEFAULT_MU）
 * @return std::vector<std::array<double, 8>> 軌道データ [t, x, y, z, vx, vy,
 * vz, jacobi_error]
 */
std::vector<std::array<double, 8>>
integrate_crtbp_adaptive(const state &x0, double t0, double tf, double mu,
                         double dt_init = 0.01, double abs_err = 1.0e-10,
                         double rel_err = 1.0e-10);

/**
 * @brief 高精度で円制限三体問題の運動方程式を積分する関数（Fehlberg78法）
 * @param x0 初期状態 [x, y, z, vx, vy, vz]
 * @param t0 初期時刻
 * @param tf 終了時刻
 * @param dt 時間刻み
 * @param abs_err 絶対誤差許容値
 * @param rel_err 相対誤差許容値
 * @param mu 質量パラメータ（デフォルト: DEFAULT_MU）
 * @return std::vector<std::array<double, 8>> 軌道データ [t, x, y, z, vx, vy,
 * vz, jacobi_error]
 */
std::vector<std::array<double, 8>>
integrate_crtbp_fehlberg78(const state &x0, double t0, double tf, double dt,
                           double mu, double abs_err = 1.0e-12,
                           double rel_err = 1.0e-12);

} // namespace crtbp
} // namespace rtbp

#endif // RTBP_HPP