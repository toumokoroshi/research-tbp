/**
 * @file periodic_orbit_stability.hpp
 * @brief モノドロミー行列計算と安定性解析の実装
 *
 * @par 理論的背景
 * 周期軌道の安定性は、状態遷移行列 (STM) \f$ \Phi(t, 0) \f$ を1周期分
 * 積分したモノドロミー行列 \f$ M = \Phi(T, 0) \f$ の固有値解析により判定する。
 * STMは変分方程式 \f$ \dot{\Phi} = A(t)\Phi \f$ を満たし、\f$ \Phi(0,0) = I \f$ である。
 *
 * @par 参考文献
 * - [Parker2006] Parker, J.S., Lo, M.W. "Practical Methods for Libration Point
 *   Orbit Analysis" - STMの数値計算手法
 * - [Koon2011] Koon, W.S. et al. "Dynamical Systems, the Three-Body Problem
 *   and Space Mission Design" (Marsden Books, 2011) - Ch.6: 安定性解析
 * - [Szebehely1967] Szebehely, V. "Theory of Orbits" (Academic Press, 1967)
 *   Ch.5: 変分方程式の導出
 * - [Meyer2009] Meyer, K.R. et al. "Introduction to Hamiltonian Dynamical
 *   Systems" (Springer, 2009) - シンプレクティック構造と固有値制約
 */

#ifndef PERIODIC_ORBIT_STABILITY_HPP
#define PERIODIC_ORBIT_STABILITY_HPP

#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <complex>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "periodic_orbit.hpp"

namespace periodic_orbit {

using namespace my_type;

// ---------------------------------------------------------------------------
// モノドロミー行列の計算（STM: State Transition Matrix）
// ---------------------------------------------------------------------------

/**
 * @brief 変分方程式の状態（軌道 + STM）
 */
template <typename ScalarType>
struct STMState {
  State<ScalarType> orbit;                       // 主軌道
  std::array<std::array<ScalarType, 6>, 6> stm;  // State Transition Matrix

  STMState() {
    // STMを単位行列で初期化
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        stm[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }
  }

  // Boost.Odeint用の演算子
  STMState operator+(const STMState& other) const {
    STMState result;
    result.orbit = orbit + other.orbit;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        result.stm[i][j] = stm[i][j] + other.stm[i][j];
      }
    }
    return result;
  }

  STMState operator*(ScalarType scalar) const {
    STMState result;
    result.orbit = orbit * scalar;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        result.stm[i][j] = stm[i][j] * scalar;
      }
    }
    return result;
  }

  STMState operator/(ScalarType scalar) const {
    STMState result;
    result.orbit = orbit * (1.0 / scalar);
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        result.stm[i][j] = stm[i][j] / scalar;
      }
    }
    return result;
  }

  // 要素ごとの除算（Boost.Odeintの誤差推定で必要）
  STMState operator/(const STMState& other) const {
    STMState result;
    result.orbit.x = orbit.x / other.orbit.x;
    result.orbit.y = orbit.y / other.orbit.y;
    result.orbit.z = orbit.z / other.orbit.z;
    result.orbit.vx = orbit.vx / other.orbit.vx;
    result.orbit.vy = orbit.vy / other.orbit.vy;
    result.orbit.vz = orbit.vz / other.orbit.vz;

    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        result.stm[i][j] = stm[i][j] / other.stm[i][j];
      }
    }

    return result;
  }

  STMState& operator+=(const STMState& other) {
    orbit += other.orbit;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        stm[i][j] += other.stm[i][j];
      }
    }
    return *this;
  }

  // ノルム計算（Boost.Odeintの適応ステップ制御用）
  ScalarType norm_inf() const {
    ScalarType max_val = 0.0;

    // 軌道部分のノルム
    max_val = std::max(max_val, std::abs(orbit.x));
    max_val = std::max(max_val, std::abs(orbit.y));
    max_val = std::max(max_val, std::abs(orbit.z));
    max_val = std::max(max_val, std::abs(orbit.vx));
    max_val = std::max(max_val, std::abs(orbit.vy));
    max_val = std::max(max_val, std::abs(orbit.vz));

    // STM部分のノルム
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        max_val = std::max(max_val, std::abs(stm[i][j]));
      }
    }

    return max_val;
  }
};

// abs関数の定義（Boost.Odeintの誤差推定用）
template <typename ScalarType>
STMState<ScalarType> abs(const STMState<ScalarType>& state) {
  STMState<ScalarType> result;
  result.orbit.x = std::abs(state.orbit.x);
  result.orbit.y = std::abs(state.orbit.y);
  result.orbit.z = std::abs(state.orbit.z);
  result.orbit.vx = std::abs(state.orbit.vx);
  result.orbit.vy = std::abs(state.orbit.vy);
  result.orbit.vz = std::abs(state.orbit.vz);

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      result.stm[i][j] = std::abs(state.stm[i][j]);
    }
  }

  return result;
}

template <typename ScalarType>
STMState<ScalarType> operator*(ScalarType scalar, const STMState<ScalarType>& state) {
  return state * scalar;
}

// スカラー加算演算子（Boost.Odeintの誤差推定で必要）
template <typename ScalarType>
STMState<ScalarType> operator+(const STMState<ScalarType>& state, ScalarType scalar) {
  STMState<ScalarType> result;
  result.orbit.x = state.orbit.x + scalar;
  result.orbit.y = state.orbit.y + scalar;
  result.orbit.z = state.orbit.z + scalar;
  result.orbit.vx = state.orbit.vx + scalar;
  result.orbit.vy = state.orbit.vy + scalar;
  result.orbit.vz = state.orbit.vz + scalar;

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      result.stm[i][j] = state.stm[i][j] + scalar;
    }
  }

  return result;
}

template <typename ScalarType>
STMState<ScalarType> operator+(ScalarType scalar, const STMState<ScalarType>& state) {
  return state + scalar;
}

/**
 * @brief CR3BPの変分方程式（軌道 + STM）
 *
 * @par 変分方程式の構造
 * 状態遷移行列 \f$ \Phi \f$ の時間発展は、ヤコビ行列 \f$ A(t) \f$ を用いて:
 * \f[ \frac{d\Phi}{dt} = A(t) \Phi, \quad \Phi(0) = I \f]
 * ここで \f$ A(t) \f$ は運動方程式のヤコビ行列:
 * \f[
 *   A = \begin{pmatrix}
 *     0_{3\times 3} & I_{3\times 3} \\
 *     \nabla^2 U & 2\Omega
 *   \end{pmatrix}
 * \f]
 *
 * @par ポテンシャルのHessian
 * 有効ポテンシャルの2階偏微分 \f$ \nabla^2 U \f$ の各成分は、
 * 主天体・副天体からの距離 \f$ r_1, r_2 \f$ を用いて:
 * \f[
 *   U_{xx} = 1 - \frac{1-\mu}{r_1^3} - \frac{\mu}{r_2^3}
 *          + 3(1-\mu)\frac{(x+\mu)^2}{r_1^5} + 3\mu\frac{(x-1+\mu)^2}{r_2^5}
 * \f]
 * 等。\f$ 2\Omega \f$ はコリオリ行列 \f$ \Omega = \begin{pmatrix} 0&1&0 \\ -1&0&0 \\ 0&0&0 \end{pmatrix} \f$ である。
 *
 * @see [Szebehely1967] Ch.5: 変分方程式の導出
 * @see [Parker2006] STMの数値積分手法
 */
template <typename ScalarType>
class VariationalEquation {
 private:
  ScalarType mu_;

 public:
  VariationalEquation(ScalarType mu) : mu_(mu) {}

  void operator()(const STMState<ScalarType>& s, STMState<ScalarType>& dsdt, ScalarType /*t*/) {
    const auto& orbit = s.orbit;

    // 主軌道の運動方程式
    ScalarType x1 = orbit.x + mu_;
    ScalarType x2 = orbit.x - (1.0 - mu_);
    ScalarType r1_sq = x1 * x1 + orbit.y * orbit.y + orbit.z * orbit.z;
    ScalarType r2_sq = x2 * x2 + orbit.y * orbit.y + orbit.z * orbit.z;
    ScalarType r1 = std::sqrt(r1_sq);
    ScalarType r2 = std::sqrt(r2_sq);
    ScalarType r1_inv3 = 1.0 / (r1_sq * r1);
    ScalarType r2_inv3 = 1.0 / (r2_sq * r2);
    ScalarType r1_inv5 = 1.0 / (r1_sq * r1_sq * r1);
    ScalarType r2_inv5 = 1.0 / (r2_sq * r2_sq * r2);

    dsdt.orbit.x = orbit.vx;
    dsdt.orbit.y = orbit.vy;
    dsdt.orbit.z = orbit.vz;
    dsdt.orbit.vx = 2.0 * orbit.vy + orbit.x - (1.0 - mu_) * x1 * r1_inv3 - mu_ * x2 * r2_inv3;
    dsdt.orbit.vy =
        -2.0 * orbit.vx + orbit.y - (1.0 - mu_) * orbit.y * r1_inv3 - mu_ * orbit.y * r2_inv3;
    dsdt.orbit.vz = -(1.0 - mu_) * orbit.z * r1_inv3 - mu_ * orbit.z * r2_inv3;

    // ポテンシャルの2階微分（Hessian \f$ \nabla^2 U \f$ ）
    ScalarType Uxx = 1.0 - (1.0 - mu_) * r1_inv3 - mu_ * r2_inv3 +
                     3.0 * (1.0 - mu_) * x1 * x1 * r1_inv5 + 3.0 * mu_ * x2 * x2 * r2_inv5;
    ScalarType Uyy = 1.0 - (1.0 - mu_) * r1_inv3 - mu_ * r2_inv3 +
                     3.0 * (1.0 - mu_) * orbit.y * orbit.y * r1_inv5 +
                     3.0 * mu_ * orbit.y * orbit.y * r2_inv5;
    ScalarType Uzz = -(1.0 - mu_) * r1_inv3 - mu_ * r2_inv3 +
                     3.0 * (1.0 - mu_) * orbit.z * orbit.z * r1_inv5 +
                     3.0 * mu_ * orbit.z * orbit.z * r2_inv5;
    ScalarType Uxy =
        3.0 * (1.0 - mu_) * x1 * orbit.y * r1_inv5 + 3.0 * mu_ * x2 * orbit.y * r2_inv5;
    ScalarType Uxz =
        3.0 * (1.0 - mu_) * x1 * orbit.z * r1_inv5 + 3.0 * mu_ * x2 * orbit.z * r2_inv5;
    ScalarType Uyz =
        3.0 * (1.0 - mu_) * orbit.y * orbit.z * r1_inv5 + 3.0 * mu_ * orbit.y * orbit.z * r2_inv5;

    // STMの時間微分: \f$ d\Phi/dt = A(t) \cdot \Phi \f$
    // \f$ A = \begin{pmatrix} 0 & I \\ \nabla^2 U & 2\Omega \end{pmatrix} \f$
    // \f$ \Omega = \begin{pmatrix} 0&1&0 \\ -1&0&0 \\ 0&0&0 \end{pmatrix} \f$

    std::array<std::array<ScalarType, 6>, 6> A_Phi;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        A_Phi[i][j] = 0.0;
      }
    }

    // A * Φの計算
    for (int j = 0; j < 6; ++j) {
      // 上半分: d(x)/dt部分 → Φの下3行をコピー
      A_Phi[0][j] = s.stm[3][j];
      A_Phi[1][j] = s.stm[4][j];
      A_Phi[2][j] = s.stm[5][j];

      // 下半分: d(v)/dt部分 → ∇²U * Φ(上3行) + 2Ω * Φ(下3行)
      A_Phi[3][j] = Uxx * s.stm[0][j] + Uxy * s.stm[1][j] + Uxz * s.stm[2][j] + 2.0 * s.stm[4][j];
      A_Phi[4][j] = Uxy * s.stm[0][j] + Uyy * s.stm[1][j] + Uyz * s.stm[2][j] - 2.0 * s.stm[3][j];
      A_Phi[5][j] = Uxz * s.stm[0][j] + Uyz * s.stm[1][j] + Uzz * s.stm[2][j];
    }

    // 結果をコピー
    dsdt.stm = A_Phi;
  }
};

/**
 * @brief モノドロミー行列の計算実装
 */
template <typename ScalarType>
std::array<std::array<ScalarType, 6>, 6> ComputeMonodromyMatrixImpl(
    const PeriodicOrbit<ScalarType>& orbit, ScalarType mu, ScalarType dt) {
  using namespace boost::numeric::odeint;

  // 変分方程式の積分器
  using Stepper = runge_kutta_dopri5<STMState<ScalarType>, ScalarType, STMState<ScalarType>,
                                     ScalarType, vector_space_algebra>;

  VariationalEquation<ScalarType> var_eq(mu);

  // 初期状態: 軌道 + 単位行列
  STMState<ScalarType> state;
  state.orbit = orbit.initial_state;
  // state.stm は constructor で単位行列に初期化済み

  // 1周期分積分（integrate_adaptiveを使用）
  ScalarType period = orbit.period;
  integrate_adaptive(make_controlled(1e-10, 1e-10, Stepper()),
                     var_eq, state, 0.0, period, dt);

  // 最終的なSTMがモノドロミー行列
  return state.stm;
}

/**
 * @brief Power法で最大固有値を計算
 */
template <typename ScalarType>
ScalarType PowerMethod(const std::array<std::array<ScalarType, 6>, 6>& matrix, int max_iter = 100) {
  std::array<ScalarType, 6> v = {1, 0, 0, 0, 0, 0};

  for (int iter = 0; iter < max_iter; ++iter) {
    std::array<ScalarType, 6> Av = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        Av[i] += matrix[i][j] * v[j];
      }
    }

    // 正規化
    ScalarType norm = 0;
    for (int i = 0; i < 6; ++i) {
      norm += Av[i] * Av[i];
    }
    norm = std::sqrt(norm);

    for (int i = 0; i < 6; ++i) {
      v[i] = Av[i] / norm;
    }
  }

  // Rayleigh商で固有値を計算
  std::array<ScalarType, 6> Av = {0, 0, 0, 0, 0, 0};
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      Av[i] += matrix[i][j] * v[j];
    }
  }

  ScalarType numerator = 0, denominator = 0;
  for (int i = 0; i < 6; ++i) {
    numerator += v[i] * Av[i];
    denominator += v[i] * v[i];
  }

  return numerator / denominator;
}

/**
 * @brief 安定性解析の実装（Eigenライブラリによる完全固有値計算）
 *
 * @par 安定性判定基準
 * ハミルトン系ではモノドロミー行列はシンプレクティックであり、
 * 固有値は逆数ペア \f$ (\lambda, 1/\lambda) \f$ で現れる。
 * 安定性の判定: 全ての固有値が \f$ |\lambda_i| \leq 1 \f$ を満たすとき安定 (SPO)。
 * \f$ |\lambda_i| > 1 \f$ の固有値が存在するとき不安定 (UPO) となり、
 * 対応する固有ベクトル方向に不安定多様体 \f$ W^u \f$ が存在する。
 *
 * @see [Meyer2009] ハミルトン系のシンプレクティック構造と固有値制約
 * @see [Koon2011] Ch.6: 周期軌道の安定性指数
 */
template <typename ScalarType>
void AnalyzeStabilityImpl(PeriodicOrbit<ScalarType>* orbit, ScalarType eigenvalue_threshold) {
  if (orbit->monodromy_matrix[0][0] == 0.0 && orbit->monodromy_matrix[1][1] == 0.0) {
    throw std::runtime_error("AnalyzeStability: Monodromy matrix not computed");
  }

  // std::array → Eigen::Matrix への変換
  Eigen::Matrix<ScalarType, 6, 6> M;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      M(i, j) = orbit->monodromy_matrix[i][j];
    }
  }

  // 固有値計算
  Eigen::EigenSolver<Eigen::Matrix<ScalarType, 6, 6>> solver(M);
  auto eigenvalues = solver.eigenvalues();

  // 固有値をPeriodicOrbitに格納
  orbit->eigenvalues.clear();
  orbit->eigenvalues.reserve(6);
  ScalarType max_abs = 0.0;

  for (int i = 0; i < 6; ++i) {
    std::complex<ScalarType> ev(eigenvalues[i].real(), eigenvalues[i].imag());
    orbit->eigenvalues.push_back(ev);

    ScalarType abs_ev = std::abs(ev);
    if (abs_ev > max_abs) {
      max_abs = abs_ev;
    }
  }

  orbit->stability_index = max_abs;

  // 安定性判定: 全ての固有値の絶対値が1以下なら安定
  orbit->is_stable = (max_abs <= eigenvalue_threshold);
  orbit->stability_computed = true;
}

/**
 * @brief 分岐方向ベクトルを抽出
 * @details ピッチフォーク分岐（\f$ \lambda = +1 \f$ 通過）の固有ベクトルを返す。
 *
 * @par 物理的意味
 * Lyapunov軌道族の継続過程でモノドロミー行列の固有値が \f$ +1 \f$ を
 * 通過するとき、ピッチフォーク分岐が発生し、新たな軌道族（Halo軌道）が
 * 分岐する。\f$ \lambda = +1 \f$ に対応する固有ベクトルはz成分を持ち、
 * この方向への摂動がHalo軌道の生成に用いられる。
 *
 * @param orbit 分岐点の周期軌道（モノドロミー行列計算済み）
 * @return z成分を持つ分岐方向ベクトル (6次元、正規化済み)
 *
 * @see [Doedel2007] Doedel, E.J. et al. "Elemental Periodic Orbits Associated
 *   with the Libration Points in the CR3BP" (Int. J. Bifurcation Chaos, 17, 2007)
 * @see [Howell1984] Howell, K.C. "Three-Dimensional, Periodic, 'Halo' Orbits"
 *   (Celestial Mechanics, 32, 53-71, 1984) DOI:10.1007/BF01358403
 */
template <typename ScalarType>
std::array<ScalarType, 6> ComputeBifurcationEigenvector(const PeriodicOrbit<ScalarType>& orbit) {
  // std::array → Eigen::Matrix への変換
  Eigen::Matrix<ScalarType, 6, 6> M;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      M(i, j) = orbit.monodromy_matrix[i][j];
    }
  }

  // 固有値・固有ベクトル計算
  Eigen::EigenSolver<Eigen::Matrix<ScalarType, 6, 6>> solver(M);
  auto eigenvalues = solver.eigenvalues();
  auto eigenvectors = solver.eigenvectors();

  // λ ≈ +1 に最も近い固有値を探す
  int closest_idx = 0;
  ScalarType min_dist = std::abs(eigenvalues[0] - std::complex<ScalarType>(1.0, 0.0));

  for (int i = 1; i < 6; ++i) {
    ScalarType dist = std::abs(eigenvalues[i] - std::complex<ScalarType>(1.0, 0.0));
    if (dist < min_dist) {
      min_dist = dist;
      closest_idx = i;
    }
  }

  // 対応する固有ベクトルを抽出（実部のみ）
  std::array<ScalarType, 6> bifurcation_vector;
  for (int i = 0; i < 6; ++i) {
    bifurcation_vector[i] = eigenvectors.col(closest_idx)(i).real();
  }

  // 正規化
  ScalarType norm = 0.0;
  for (int i = 0; i < 6; ++i) {
    norm += bifurcation_vector[i] * bifurcation_vector[i];
  }
  norm = std::sqrt(norm);
  if (norm > 1e-14) {
    for (int i = 0; i < 6; ++i) {
      bifurcation_vector[i] /= norm;
    }
  }

  return bifurcation_vector;
}

}  // namespace periodic_orbit

// Boost.Odeint用の特殊化（vector_space_algebraサポート）
namespace boost {
namespace numeric {
namespace odeint {

template <typename ScalarType>
struct vector_space_norm_inf<periodic_orbit::STMState<ScalarType>> {
  typedef ScalarType result_type;
  result_type operator()(const periodic_orbit::STMState<ScalarType>& s) const {
    return s.norm_inf();
  }
};

}  // namespace odeint
}  // namespace numeric
}  // namespace boost

#endif  // PERIODIC_ORBIT_STABILITY_HPP
