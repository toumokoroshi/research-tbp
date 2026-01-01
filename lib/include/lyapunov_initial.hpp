/**
 * @file lyapunov_initial.hpp
 * @brief Lyapunov軌道の初期条件生成ライブラリ
 * @details L1/L2/L3点周りのLyapunov軌道の解析的初期条件を生成
 *          研究計画書 Step 1-1（線形化による解析的近似）に対応
 * @date 2026-01-01
 */

#ifndef LYAPUNOV_INITIAL_HPP
#define LYAPUNOV_INITIAL_HPP

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <stdexcept>

#include "rtbp.hpp"

namespace lyapunov {

using namespace my_type;

// ---------------------------------------------------------------------------
// 列挙型定義
// ---------------------------------------------------------------------------

/**
 * @brief ラグランジュ点の種類
 */
enum class LagrangePoint { L1, L2, L3 };

/**
 * @brief 軌道族の種類（将来の拡張用）
 */
enum class OrbitFamilyType {
  LYAPUNOV_L1,
  LYAPUNOV_L2,
  LYAPUNOV_L3,
  HALO_L1_NORTH,
  HALO_L1_SOUTH,
  HALO_L2_NORTH,
  HALO_L2_SOUTH,
  DRO,       // Distant Retrograde Orbit - 将来拡張
  RESONANT,  // 共鳴軌道 - 将来拡張
};

// ---------------------------------------------------------------------------
// 構造体定義
// ---------------------------------------------------------------------------

/**
 * @brief Lyapunov軌道の係数を保持する構造体
 * @details L点周りの線形化から得られるパラメータ
 */
template <typename ScalarType>
struct LyapunovCoefficients {
  ScalarType x_L;       ///< L点のx座標
  ScalarType c2;        ///< 線形化係数 c2
  ScalarType omega_xy;  ///< 面内周波数
  ScalarType omega_z;   ///< 面外周波数
  ScalarType kappa;     ///< 振幅比係数 κ = (ω² + 1 + 2c₂) / (2ω)
  ScalarType lambda;    ///< サドル型固有値（実部）
};

/**
 * @brief 継続法の初期軌道を生成するための抽象基底クラス
 * @details DRO等への拡張のためのインターフェース
 */
template <typename ScalarType>
class InitialOrbitGenerator {
 public:
  virtual ~InitialOrbitGenerator() = default;

  /**
   * @brief 初期条件を生成
   * @param amplitude 振幅パラメータ
   * @return 初期状態ベクトル (x, y, z, vx, vy, vz)
   */
  virtual State<ScalarType> GenerateInitialGuess(ScalarType amplitude) = 0;

  /**
   * @brief 推定周期を取得
   * @param amplitude 振幅パラメータ
   * @return 推定周期
   */
  virtual ScalarType EstimatePeriod(ScalarType amplitude) = 0;

  /**
   * @brief 族の種類を取得
   */
  virtual OrbitFamilyType GetFamilyType() const = 0;
};

// ---------------------------------------------------------------------------
// 関数宣言
// ---------------------------------------------------------------------------

/**
 * @brief ラグランジュ点L1のx座標を計算（ニュートン法）
 * @param mu 質量比パラメータ
 * @param tolerance 収束許容誤差
 * @param max_iter 最大反復回数
 * @return L1点のx座標
 */
template <typename ScalarType>
ScalarType ComputeL1Position(ScalarType mu, ScalarType tolerance = 1e-14, int max_iter = 100);

/**
 * @brief ラグランジュ点L2のx座標を計算（ニュートン法）
 * @param mu 質量比パラメータ
 * @param tolerance 収束許容誤差
 * @param max_iter 最大反復回数
 * @return L2点のx座標
 */
template <typename ScalarType>
ScalarType ComputeL2Position(ScalarType mu, ScalarType tolerance = 1e-14, int max_iter = 100);

/**
 * @brief ラグランジュ点L3のx座標を計算（ニュートン法）
 * @param mu 質量比パラメータ
 * @param tolerance 収束許容誤差
 * @param max_iter 最大反復回数
 * @return L3点のx座標
 */
template <typename ScalarType>
ScalarType ComputeL3Position(ScalarType mu, ScalarType tolerance = 1e-14, int max_iter = 100);

/**
 * @brief 指定されたラグランジュ点のx座標を計算
 * @param point ラグランジュ点
 * @param mu 質量比パラメータ
 * @return L点のx座標
 */
template <typename ScalarType>
ScalarType ComputeLagrangePointPosition(LagrangePoint point, ScalarType mu);

/**
 * @brief Lyapunov軌道の係数を計算
 * @param point ラグランジュ点
 * @param mu 質量比パラメータ
 * @return Lyapunov係数構造体
 */
template <typename ScalarType>
LyapunovCoefficients<ScalarType> ComputeLyapunovCoefficients(LagrangePoint point, ScalarType mu);

/**
 * @brief Lyapunov軌道の初期条件を生成（0次近似）
 * @param point ラグランジュ点
 * @param mu 質量比パラメータ
 * @param amplitude 振幅 A_x（無次元）
 * @return 初期状態ベクトル (x, y, z, vx, vy, vz)
 */
template <typename ScalarType>
State<ScalarType> GenerateLyapunovInitialGuess(LagrangePoint point, ScalarType mu,
                                               ScalarType amplitude);

/**
 * @brief Lyapunov軌道の推定周期を計算
 * @param point ラグランジュ点
 * @param mu 質量比パラメータ
 * @return 推定周期（無次元時間）
 */
template <typename ScalarType>
ScalarType EstimateLyapunovPeriod(LagrangePoint point, ScalarType mu);

// ---------------------------------------------------------------------------
// 実装部
// ---------------------------------------------------------------------------

/**
 * @brief L1点位置計算（ニュートン法）
 * @details L1: 太陽と地球の間に位置
 *          平衡点条件: x - (1-μ)/r1² + μ/r2² = 0
 *          r1 = x + μ (太陽からの距離), r2 = (1-μ) - x (地球からの距離)
 */
template <typename ScalarType>
ScalarType ComputeL1Position(ScalarType mu, ScalarType tolerance, int max_iter) {
  // 初期推定値: 地球からHill球半径程度内側
  ScalarType gamma = std::pow(mu / 3.0, 1.0 / 3.0);
  ScalarType x = 1.0 - mu - gamma;  // 地球より少し太陽側

  for (int i = 0; i < max_iter; ++i) {
    ScalarType r1 = x + mu;          // 太陽からの距離（正）
    ScalarType r2 = (1.0 - mu) - x;  // 地球からの距離（正）

    if (r1 < 1e-14 || r2 < 1e-14) {
      throw std::runtime_error("ComputeL1Position: Too close to primary");
    }

    ScalarType r1_2 = r1 * r1;
    ScalarType r2_2 = r2 * r2;
    ScalarType r1_3 = r1_2 * r1;
    ScalarType r2_3 = r2_2 * r2;

    // 平衡点条件: f = x - (1-μ)/r1² + μ/r2² = 0
    ScalarType f = x - (1.0 - mu) / r1_2 + mu / r2_2;

    // 導関数: df/dx = 1 + 2(1-μ)/r1³ + 2μ/r2³
    ScalarType df = 1.0 + 2.0 * (1.0 - mu) / r1_3 + 2.0 * mu / r2_3;

    ScalarType dx = -f / df;
    x += dx;

    if (std::abs(dx) < tolerance) {
      return x;
    }
  }

  throw std::runtime_error("ComputeL1Position: Failed to converge");
}

/**
 * @brief L2点位置計算（ニュートン法）
 * @details L2: 月の外側に位置
 */
template <typename ScalarType>
ScalarType ComputeL2Position(ScalarType mu, ScalarType tolerance, int max_iter) {
  // 初期推定値: 月の外側に Hill 球半径程度
  ScalarType gamma = std::pow(mu / 3.0, 1.0 / 3.0);
  ScalarType x = 1.0 - mu + gamma;

  for (int i = 0; i < max_iter; ++i) {
    ScalarType r1 = x + mu;
    ScalarType r2 = x - (1.0 - mu);

    ScalarType f = x - (1.0 - mu) / (r1 * r1) * std::copysign(1.0, r1) -
                   mu / (r2 * r2) * std::copysign(1.0, r2);

    ScalarType df =
        1.0 + 2.0 * (1.0 - mu) / std::abs(r1 * r1 * r1) + 2.0 * mu / std::abs(r2 * r2 * r2);

    ScalarType dx = -f / df;
    x += dx;

    if (std::abs(dx) < tolerance) {
      return x;
    }
  }

  throw std::runtime_error("ComputeL2Position: Failed to converge");
}

/**
 * @brief L3点位置計算（ニュートン法）
 * @details L3: 地球の反対側（太陽側）に位置
 */
template <typename ScalarType>
ScalarType ComputeL3Position(ScalarType mu, ScalarType tolerance, int max_iter) {
  // 初期推定値: x ≈ -1 - 5μ/12
  ScalarType x = -1.0 - 5.0 * mu / 12.0;

  for (int i = 0; i < max_iter; ++i) {
    ScalarType r1 = x + mu;
    ScalarType r2 = x - (1.0 - mu);

    ScalarType f = x + (1.0 - mu) / (r1 * r1) + mu / (r2 * r2);

    ScalarType df =
        1.0 + 2.0 * (1.0 - mu) / std::abs(r1 * r1 * r1) + 2.0 * mu / std::abs(r2 * r2 * r2);

    ScalarType dx = -f / df;
    x += dx;

    if (std::abs(dx) < tolerance) {
      return x;
    }
  }

  throw std::runtime_error("ComputeL3Position: Failed to converge");
}

template <typename ScalarType>
ScalarType ComputeLagrangePointPosition(LagrangePoint point, ScalarType mu) {
  switch (point) {
    case LagrangePoint::L1:
      return ComputeL1Position(mu);
    case LagrangePoint::L2:
      return ComputeL2Position(mu);
    case LagrangePoint::L3:
      return ComputeL3Position(mu);
    default:
      throw std::invalid_argument("Unknown Lagrange point");
  }
}

/**
 * @brief Lyapunov係数の計算
 * @details 研究計画書のStep 1-1に対応
 *          L点周りの線形化: ẍ - 2ẏ = (1 + 2c₂)x, ÿ + 2ẋ = (1 - c₂)y
 *          c₂ = (1-μ)/|L-μ|³ + μ/|L-(1-μ)|³
 */
template <typename ScalarType>
LyapunovCoefficients<ScalarType> ComputeLyapunovCoefficients(LagrangePoint point, ScalarType mu) {
  LyapunovCoefficients<ScalarType> coeff;

  // L点位置を計算
  coeff.x_L = ComputeLagrangePointPosition(point, mu);

  // r1, r2の距離を計算
  ScalarType r1 = std::abs(coeff.x_L + mu);
  ScalarType r2 = std::abs(coeff.x_L - (1.0 - mu));

  // c2係数の計算
  coeff.c2 = (1.0 - mu) / (r1 * r1 * r1) + mu / (r2 * r2 * r2);

  // 特性方程式から固有値を計算
  // 面内運動の特性方程式: λ⁴ + (2 - c₂)λ² + (1 + 2c₂)(1 - c₂) = 0
  // これは λ² に関する2次方程式
  ScalarType a = 1.0;
  ScalarType b = 2.0 - coeff.c2;
  ScalarType c = (1.0 + 2.0 * coeff.c2) * (1.0 - coeff.c2);

  ScalarType discriminant = b * b - 4.0 * a * c;

  // サドル点（L1, L2, L3）では判別式 > 0
  if (discriminant < 0) {
    throw std::runtime_error("ComputeLyapunovCoefficients: Unexpected eigenvalue structure");
  }

  // λ² の2つの根
  ScalarType lambda2_1 = (-b + std::sqrt(discriminant)) / (2.0 * a);
  ScalarType lambda2_2 = (-b - std::sqrt(discriminant)) / (2.0 * a);

  // 一方は正（サドル）、一方は負（中心）
  if (lambda2_1 > 0 && lambda2_2 < 0) {
    coeff.lambda = std::sqrt(lambda2_1);     // サドル型固有値
    coeff.omega_xy = std::sqrt(-lambda2_2);  // 面内振動周波数
  } else if (lambda2_1 < 0 && lambda2_2 > 0) {
    coeff.lambda = std::sqrt(lambda2_2);
    coeff.omega_xy = std::sqrt(-lambda2_1);
  } else {
    throw std::runtime_error("ComputeLyapunovCoefficients: Unexpected eigenvalue signs");
  }

  // 面外振動周波数
  coeff.omega_z = std::sqrt(coeff.c2);

  // 振幅比係数 κ
  // κ = (ω² + 1 + 2c₂) / (2ω)
  coeff.kappa = (coeff.omega_xy * coeff.omega_xy + 1.0 + 2.0 * coeff.c2) / (2.0 * coeff.omega_xy);

  return coeff;
}

/**
 * @brief Lyapunov軌道の初期条件生成
 * @details 小振幅近似（0次近似）
 *          x(t) = x_L + A_x cos(ω t)
 *          y(t) = -κ A_x sin(ω t)
 *          z(t) = 0
 *
 *          t=0での初期条件:
 *          x₀ = x_L + A_x, y₀ = 0, z₀ = 0
 *          ẋ₀ = 0, ẏ₀ = -κ ω A_x, ż₀ = 0
 */
template <typename ScalarType>
State<ScalarType> GenerateLyapunovInitialGuess(LagrangePoint point, ScalarType mu,
                                               ScalarType amplitude) {
  // 係数を計算
  auto coeff = ComputeLyapunovCoefficients(point, mu);

  State<ScalarType> state;
  state.x = coeff.x_L + amplitude;
  state.y = 0.0;
  state.z = 0.0;
  state.vx = 0.0;
  state.vy = -coeff.kappa * coeff.omega_xy * amplitude;
  state.vz = 0.0;

  return state;
}

/**
 * @brief Lyapunov軌道の推定周期
 */
template <typename ScalarType>
ScalarType EstimateLyapunovPeriod(LagrangePoint point, ScalarType mu) {
  auto coeff = ComputeLyapunovCoefficients(point, mu);
  return 2.0 * std::numbers::pi_v<ScalarType> / coeff.omega_xy;
}

// ---------------------------------------------------------------------------
// Lyapunov軌道ジェネレータクラス（InitialOrbitGeneratorの具体実装）
// ---------------------------------------------------------------------------

template <typename ScalarType>
class LyapunovOrbitGenerator : public InitialOrbitGenerator<ScalarType> {
 public:
  LyapunovOrbitGenerator(LagrangePoint point, ScalarType mu)
      : point_(point), mu_(mu), coeff_(ComputeLyapunovCoefficients<ScalarType>(point, mu)) {}

  State<ScalarType> GenerateInitialGuess(ScalarType amplitude) override {
    return GenerateLyapunovInitialGuess(point_, mu_, amplitude);
  }

  ScalarType EstimatePeriod(ScalarType /*amplitude*/) override {
    // 線形近似では振幅に依存しない
    return 2.0 * std::numbers::pi_v<ScalarType> / coeff_.omega_xy;
  }

  OrbitFamilyType GetFamilyType() const override {
    switch (point_) {
      case LagrangePoint::L1:
        return OrbitFamilyType::LYAPUNOV_L1;
      case LagrangePoint::L2:
        return OrbitFamilyType::LYAPUNOV_L2;
      case LagrangePoint::L3:
        return OrbitFamilyType::LYAPUNOV_L3;
      default:
        return OrbitFamilyType::LYAPUNOV_L1;
    }
  }

  const LyapunovCoefficients<ScalarType>& GetCoefficients() const { return coeff_; }

 private:
  LagrangePoint point_;
  ScalarType mu_;
  LyapunovCoefficients<ScalarType> coeff_;
};

// ---------------------------------------------------------------------------
// ファクトリ関数（拡張性のため）
// ---------------------------------------------------------------------------

/**
 * @brief 軌道ジェネレータを作成するファクトリ関数
 * @param family_type 軌道族の種類
 * @param mu 質量比
 * @return ジェネレータへのunique_ptr
 */
template <typename ScalarType>
std::unique_ptr<InitialOrbitGenerator<ScalarType>> CreateOrbitGenerator(OrbitFamilyType family_type,
                                                                        ScalarType mu) {
  switch (family_type) {
    case OrbitFamilyType::LYAPUNOV_L1:
      return std::make_unique<LyapunovOrbitGenerator<ScalarType>>(LagrangePoint::L1, mu);
    case OrbitFamilyType::LYAPUNOV_L2:
      return std::make_unique<LyapunovOrbitGenerator<ScalarType>>(LagrangePoint::L2, mu);
    case OrbitFamilyType::LYAPUNOV_L3:
      return std::make_unique<LyapunovOrbitGenerator<ScalarType>>(LagrangePoint::L3, mu);

    // 将来の拡張用
    case OrbitFamilyType::HALO_L1_NORTH:
    case OrbitFamilyType::HALO_L1_SOUTH:
    case OrbitFamilyType::HALO_L2_NORTH:
    case OrbitFamilyType::HALO_L2_SOUTH:
      throw std::runtime_error("Halo orbit generators not implemented yet");

    case OrbitFamilyType::DRO:
      throw std::runtime_error("DRO generator not implemented yet");

    case OrbitFamilyType::RESONANT:
      throw std::runtime_error("Resonant orbit generator not implemented yet");

    default:
      throw std::invalid_argument("Unknown orbit family type");
  }
}

}  // namespace lyapunov

#endif  // LYAPUNOV_INITIAL_HPP
