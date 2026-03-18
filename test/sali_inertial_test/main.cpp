/**
 * @file main.cpp
 * @brief 慣性系 vs 回転座標系 SALI比較テスト
 *
 * GEO/LEO軌道のSALI振動が回転座標系の座標効果なのか、
 * CR3BPの物理的効果なのかを切り分けるテスト
 *
 * テスト1: 純粋二体問題のSALI（コントロール）
 * テスト2: CR3BP回転座標系のSALI（シンプレクティック）
 * テスト3: CR3BP回転座標系のSALI（RK Dopri5）
 * 補足: ヘシアンとコリオリ項のスケール確認
 *
 * @date 2026-02-10
 */

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

#include <utils.hpp>

#include "rtbp.hpp"

using namespace my_type;
using namespace crtbp;
using namespace param;

namespace fs = std::filesystem;

// ============================================================================
// 二体問題専用のポテンシャル関数（太陽項を除去）
// ============================================================================

/**
 * @brief 二体問題用の勾配とヘシアン計算（地球のmu項のみ）
 *
 * @details CR3BPの有効ポテンシャル U = -(1-mu)/r1 - mu/r2 から
 *          太陽項 -(1-mu)/r1 を除去し、地球項 -mu/r2 のみ残す。
 *          コリオリ項も除去するため、遠心力 (1/2)(x^2+y^2) も含まない。
 *
 * @param mu 質量比
 * @param x, y, z 第三体の位置
 * @param grad_U 出力: 勾配ベクトル
 * @param hessian_U 出力: ヘシアン行列
 */
template <typename ScalarType>
void CalculateGradientAndHessianU_TwoBody(const ScalarType mu, ScalarType x, ScalarType y,
                                           ScalarType z, ScalarType* grad_U,
                                           HessianMatrix<ScalarType>* hessian_U) {
  const ScalarType x2 = x - (1.0 - mu);
  const ScalarType r2_sq = x2 * x2 + y * y + z * z;

  constexpr ScalarType kMinDistSq = 1e-16;
  if (r2_sq < kMinDistSq) {
    throw std::runtime_error("Position too close to Earth in TwoBody potential.");
  }

  const ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));
  const ScalarType r2_inv5 = r2_inv3 / r2_sq;

  // 勾配: dU/dq = mu * (q - q_earth) / r2^3
  grad_U[0] = mu * x2 * r2_inv3;
  grad_U[1] = mu * y * r2_inv3;
  grad_U[2] = mu * z * r2_inv3;

  // ヘシアン
  ScalarType term;
  term = mu * (r2_inv3 - 3.0 * x2 * x2 * r2_inv5);
  hessian_U->hxx = term;
  term = mu * (r2_inv3 - 3.0 * y * y * r2_inv5);
  hessian_U->hyy = term;
  term = mu * (r2_inv3 - 3.0 * z * z * r2_inv5);
  hessian_U->hzz = term;
  hessian_U->hxy = mu * (-3.0 * x2 * y * r2_inv5);
  hessian_U->hxz = mu * (-3.0 * x2 * z * r2_inv5);
  hessian_U->hyz = mu * (-3.0 * y * z * r2_inv5);
}

/**
 * @brief 二体問題用の勾配のみ計算
 */
template <typename ScalarType>
void CalculateGradientU_TwoBody(const ScalarType mu, ScalarType x, ScalarType y, ScalarType z,
                                 ScalarType* grad_U) {
  const ScalarType x2 = x - (1.0 - mu);
  const ScalarType r2_sq = x2 * x2 + y * y + z * z;

  constexpr ScalarType kMinDistSq = 1e-16;
  if (r2_sq < kMinDistSq) {
    throw std::runtime_error("Position too close to Earth in TwoBody gradient.");
  }

  const ScalarType r2_inv3 = 1.0 / (r2_sq * std::sqrt(r2_sq));

  grad_U[0] = mu * x2 * r2_inv3;
  grad_U[1] = mu * y * r2_inv3;
  grad_U[2] = mu * z * r2_inv3;
}

// ============================================================================
// 二体問題用シンプレクティックSALI積分器
// コリオリ項（RotateStateSALI）を省略
// ============================================================================

/**
 * @brief 二体問題用Kick（運動量更新）ステップ
 */
template <typename ScalarType>
void UpdateMomentum_TwoBody(const ScalarType mu, CanonicalState<ScalarType>* state, ScalarType dt) {
  ScalarType grad_U[3];
  CalculateGradientU_TwoBody(mu, state->qx, state->qy, state->qz, grad_U);
  state->px -= dt * grad_U[0];
  state->py -= dt * grad_U[1];
  state->pz -= dt * grad_U[2];
}

/**
 * @brief 二体問題用KickステップSALI版
 */
template <typename ScalarType>
void UpdateMomentumSALI_TwoBody(const ScalarType mu, SaliState<ScalarType>* state, ScalarType dt) {
  ScalarType grad_U[3];
  HessianMatrix<ScalarType> hessian_U;
  CalculateGradientAndHessianU_TwoBody(mu, state->state.qx, state->state.qy, state->state.qz,
                                        grad_U, &hessian_U);
  // 主軌道
  state->state.px -= dt * grad_U[0];
  state->state.py -= dt * grad_U[1];
  state->state.pz -= dt * grad_U[2];

  // 偏差ベクトル1
  const ScalarType dq1x = state->w1.qx, dq1y = state->w1.qy, dq1z = state->w1.qz;
  state->w1.px -= dt * (hessian_U.hxx * dq1x + hessian_U.hxy * dq1y + hessian_U.hxz * dq1z);
  state->w1.py -= dt * (hessian_U.hxy * dq1x + hessian_U.hyy * dq1y + hessian_U.hyz * dq1z);
  state->w1.pz -= dt * (hessian_U.hxz * dq1x + hessian_U.hyz * dq1y + hessian_U.hzz * dq1z);

  // 偏差ベクトル2
  const ScalarType dq2x = state->w2.qx, dq2y = state->w2.qy, dq2z = state->w2.qz;
  state->w2.px -= dt * (hessian_U.hxx * dq2x + hessian_U.hxy * dq2y + hessian_U.hxz * dq2z);
  state->w2.py -= dt * (hessian_U.hxy * dq2x + hessian_U.hyy * dq2y + hessian_U.hyz * dq2z);
  state->w2.pz -= dt * (hessian_U.hxz * dq2x + hessian_U.hyz * dq2y + hessian_U.hzz * dq2z);
}

/**
 * @brief Position update for SALI (same as CR3BP: q += dt * p)
 */
template <typename ScalarType>
void UpdatePositionSALI_TwoBody(SaliState<ScalarType>* state, ScalarType dt) {
  // 主軌道
  state->state.qx += dt * state->state.px;
  state->state.qy += dt * state->state.py;
  state->state.qz += dt * state->state.pz;
  // w1
  state->w1.qx += dt * state->w1.px;
  state->w1.qy += dt * state->w1.py;
  state->w1.qz += dt * state->w1.pz;
  // w2
  state->w2.qx += dt * state->w2.px;
  state->w2.qy += dt * state->w2.py;
  state->w2.qz += dt * state->w2.pz;
}

/**
 * @brief 二体問題用 2次シンプレクティックSALIステップ
 * コリオリ項（RotateStateSALI）を省略した Leapfrog: K(h/2)*D(h)*K(h/2)
 */
template <typename ScalarType>
void SymplecticStepSALI_TwoBody(const ScalarType mu, SaliState<ScalarType>* state, ScalarType h) {
  // コリオリ項なし: 直接 Leapfrog
  UpdatePositionSALI_TwoBody(state, h / 2.0);
  UpdateMomentumSALI_TwoBody(mu, state, h);
  UpdatePositionSALI_TwoBody(state, h / 2.0);
}

/**
 * @brief 二体問題用 6次吉田法SALIステップ
 */
template <typename ScalarType>
void SymplecticStep6thOrderSALI_TwoBody(const ScalarType mu, SaliState<ScalarType>* state,
                                         ScalarType h) {
  constexpr ScalarType w1 = static_cast<ScalarType>(0.784513610477560);
  constexpr ScalarType w2 = static_cast<ScalarType>(0.235573213359357);
  constexpr ScalarType w3 = static_cast<ScalarType>(-1.17767998417887);
  constexpr ScalarType w4 = static_cast<ScalarType>(1.31518632068391);
  constexpr ScalarType weights[7] = {w1, w2, w3, w4, w3, w2, w1};

  for (auto c : weights) {
    SymplecticStepSALI_TwoBody(mu, state, c * h);
  }
}

// ============================================================================
// SALI時系列の計算と出力
// ============================================================================

struct SaliTimeSeriesConfig {
  std::string label;
  double x0, y0, z0;      // 初期位置（回転座標系、無次元）
  double vx0, vy0, vz0;   // 初期速度（回転座標系、無次元）
  double timestep;
  double total_time;
  int output_interval;     // 何ステップごとにCSV出力するか
};

/**
 * @brief テスト1: 二体問題SALI
 */
void RunTwoBodySALI(const SaliTimeSeriesConfig& cfg, double mu, const std::string& output_path) {
  std::ofstream ofs(output_path);
  ofs << "time,sali,w1_w2_dot,qx,qy,qz" << std::endl;
  ofs << std::scientific << std::setprecision(15);

  // 物理状態→正準状態に変換（ただし二体問題なので p = v）
  // 二体問題ではコリオリ項がないので正準運動量 = 速度
  SaliState<double> state;
  state.state = CanonicalState<double>{cfg.x0, cfg.y0, cfg.z0, cfg.vx0, cfg.vy0, cfg.vz0};
  state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  int num_steps = static_cast<int>(cfg.total_time / cfg.timestep);

  for (int step = 0; step <= num_steps; ++step) {
    if (step % cfg.output_interval == 0) {
      state.w1.Normalize();
      state.w2.Normalize();
      double norm_plus = (state.w1 + state.w2).Norm();
      double norm_minus = (state.w1 - state.w2).Norm();
      double sali = std::min(norm_plus, norm_minus);
      double w1w2_dot = state.w1.Dot(state.w2);
      double t = step * cfg.timestep;
      ofs << t << "," << sali << "," << w1w2_dot << "," << state.state.qx << "," << state.state.qy
          << "," << state.state.qz << std::endl;
    }
    if (step < num_steps) {
      SymplecticStep6thOrderSALI_TwoBody(mu, &state, cfg.timestep);
      state.w1.Normalize();
      state.w2.Normalize();
    }
  }
  ofs.close();
  std::cout << "  Output: " << output_path << std::endl;
}

/**
 * @brief テスト2: CR3BPシンプレクティックSALI
 */
void RunCR3BPSymplecticSALI(const SaliTimeSeriesConfig& cfg, double mu,
                             const std::string& output_path) {
  std::ofstream ofs(output_path);
  ofs << "time,sali,w1_w2_dot,qx,qy,qz" << std::endl;
  ofs << std::scientific << std::setprecision(15);

  // 物理状態→正準状態（CR3BPの正準変換: px = vx - y, py = vy + x）
  State<double> phys_state{cfg.x0, cfg.y0, cfg.z0, cfg.vx0, cfg.vy0, cfg.vz0};
  CanonicalState<double> canonical = ConvertToCanonical(phys_state);

  SaliState<double> state;
  state.state = canonical;
  state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  int num_steps = static_cast<int>(cfg.total_time / cfg.timestep);

  for (int step = 0; step <= num_steps; ++step) {
    if (step % cfg.output_interval == 0) {
      state.w1.Normalize();
      state.w2.Normalize();
      double norm_plus = (state.w1 + state.w2).Norm();
      double norm_minus = (state.w1 - state.w2).Norm();
      double sali = std::min(norm_plus, norm_minus);
      double w1w2_dot = state.w1.Dot(state.w2);
      double t = step * cfg.timestep;
      ofs << t << "," << sali << "," << w1w2_dot << "," << state.state.qx << "," << state.state.qy
          << "," << state.state.qz << std::endl;
    }
    if (step < num_steps) {
      SymplecticStep6thOrderSALI(mu, &state, cfg.timestep);
      state.w1.Normalize();
      state.w2.Normalize();
    }
  }
  ofs.close();
  std::cout << "  Output: " << output_path << std::endl;
}

/**
 * @brief テスト3: CR3BP RK (Dopri5) SALI
 *
 * Boost.Odeint の dopri5 を使用
 */
void RunCR3BPRKSALI(const SaliTimeSeriesConfig& cfg, double mu, const std::string& output_path,
                     const AstroConstants<double>& astro_params) {
  std::ofstream ofs(output_path);
  ofs << "time,sali,w1_w2_dot,qx,qy,qz" << std::endl;
  ofs << std::scientific << std::setprecision(15);

  State<double> phys_state{cfg.x0, cfg.y0, cfg.z0, cfg.vx0, cfg.vy0, cfg.vz0};
  CanonicalState<double> canonical = ConvertToCanonical(phys_state);

  SaliState<double> state;
  state.state = canonical;
  state.w1 = CanonicalState<double>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  state.w2 = CanonicalState<double>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  Dopri5SaliSystem<double> rk_system(astro_params);

  int num_steps = static_cast<int>(cfg.total_time / cfg.timestep);

  for (int step = 0; step <= num_steps; ++step) {
    if (step % cfg.output_interval == 0) {
      state.w1.Normalize();
      state.w2.Normalize();
      double norm_plus = (state.w1 + state.w2).Norm();
      double norm_minus = (state.w1 - state.w2).Norm();
      double sali = std::min(norm_plus, norm_minus);
      double w1w2_dot = state.w1.Dot(state.w2);
      double t = step * cfg.timestep;
      ofs << t << "," << sali << "," << w1w2_dot << "," << state.state.qx << "," << state.state.qy
          << "," << state.state.qz << std::endl;
    }
    if (step < num_steps) {
      // Simple RK4 for SaliState
      SaliState<double> k1, k2, k3, k4;
      double h = cfg.timestep;

      rk_system(state, k1, 0.0);
      SaliState<double> state_tmp = state + k1 * (h / 2.0);
      rk_system(state_tmp, k2, 0.0);
      state_tmp = state + k2 * (h / 2.0);
      rk_system(state_tmp, k3, 0.0);
      state_tmp = state + k3 * h;
      rk_system(state_tmp, k4, 0.0);

      state = state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);

      state.w1.Normalize();
      state.w2.Normalize();
    }
  }
  ofs.close();
  std::cout << "  Output: " << output_path << std::endl;
}

/**
 * @brief ヘシアンとコリオリ項のスケール比較
 */
void ReportHessianScale(double mu, double x, double y, double z, const std::string& label,
                         std::ostream& os) {
  double grad_U[3];
  HessianMatrix<double> hessian;
  CalculateGradientAndHessianU(mu, x, y, z, grad_U, &hessian);

  // ヘシアンのフロベニウスノルム
  double hessian_frobenius =
      std::sqrt(hessian.hxx * hessian.hxx + hessian.hyy * hessian.hyy + hessian.hzz * hessian.hzz +
                2.0 * (hessian.hxy * hessian.hxy + hessian.hxz * hessian.hxz +
                        hessian.hyz * hessian.hyz));

  // コリオリ行列のノルム: C = [[0, 2, 0], [-2, 0, 0], [0, 0, 0]]
  // フロベニウスノルム = sqrt(4+4) = 2*sqrt(2) ≈ 2.83
  const double kCoriolisFrobenius = 2.0 * std::sqrt(2.0);

  double r2 = calc_r2(x, y, z, mu);

  os << "=== " << label << " ===" << std::endl;
  os << std::scientific << std::setprecision(6);
  os << "  Position (x, y, z): " << x << ", " << y << ", " << z << std::endl;
  os << "  r2 (distance to Earth): " << r2 << std::endl;
  os << "  Hessian ||H_U||_F: " << hessian_frobenius << std::endl;
  os << "  Coriolis ||C||_F:  " << kCoriolisFrobenius << std::endl;
  os << "  Ratio ||H_U|| / ||C||: " << hessian_frobenius / kCoriolisFrobenius << std::endl;
  os << "  Hessian components:" << std::endl;
  os << "    Hxx=" << hessian.hxx << "  Hyy=" << hessian.hyy << "  Hzz=" << hessian.hzz
     << std::endl;
  os << "    Hxy=" << hessian.hxy << "  Hxz=" << hessian.hxz << "  Hyz=" << hessian.hyz
     << std::endl;
  os << std::endl;
}

// ============================================================================
// メイン
// ============================================================================

int main() {
  std::cout << "================================================================" << std::endl;
  std::cout << "  SALI Inertial vs Rotating Frame Comparison Test" << std::endl;
  std::cout << "================================================================" << std::endl;

  // 天文定数の読み込み
  const std::string kConfigFilePath = CONFIG_DIR;
  std::string astro_param_file = kConfigFilePath + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = utils::loadConstants<double>(astro_param_file);

  const double kGMSUN = astro_params.gm_sun;
  const double kGMEARTH = astro_params.gm_earth;
  const double kMU = kGMEARTH / (kGMEARTH + kGMSUN);

  std::cout << std::setprecision(10);
  std::cout << "  mu = " << kMU << std::endl;

  // 出力ディレクトリ
  std::string output_dir = std::string(OUTPUT_DIR) + "/sali_inertial_test";
  fs::create_directories(output_dir);

  // ============ 初期条件の設定 ============
  // GEO: 42164 km = 42164 / 149597870.7 AU ≈ 2.818e-4 (無次元)
  const double kGeoRadius_nd = 42164.0 / 149597870.7;
  // LEO (ISS): 6771 km ≈ 4.525e-5 (無次元)
  const double kLeoRadius_nd = 6771.0 / 149597870.7;

  // 地球中心の位置 (回転座標系)
  const double kEarthX = 1.0 - kMU;

  // ケプラー円軌道速度 v = sqrt(mu / r)  (二体問題の無次元系)
  const double kGeoVelocity = std::sqrt(kMU / kGeoRadius_nd);
  const double kLeoVelocity = std::sqrt(kMU / kLeoRadius_nd);

  std::cout << "\n  === Initial Conditions ===" << std::endl;
  std::cout << "  GEO: r2 = " << kGeoRadius_nd << ", v_kepler = " << kGeoVelocity << std::endl;
  std::cout << "  LEO: r2 = " << kLeoRadius_nd << ", v_kepler = " << kLeoVelocity << std::endl;

  // GEO: 地球の右側に配置、y方向に速度
  SaliTimeSeriesConfig geo_cfg;
  geo_cfg.label = "GEO";
  geo_cfg.x0 = kEarthX + kGeoRadius_nd;
  geo_cfg.y0 = 0.0;
  geo_cfg.z0 = 0.0;
  // 回転座標系での速度: v_inertial - omega×r
  // 慣性系でケプラー速度 v_kepler を y方向に持つ場合
  // 回転座標系では vx_rot = 0, vy_rot = v_kepler - 1*r (omega=1, r=distance from origin)
  // ただし、ここでは地球周りの局所的な円軌道を考える
  // 慣性系での速度: vy = v_kepler (地球重心周り)
  // 回転座標系: vy_rot = vy_inertial + omega*(x - x_barycenter)
  //           ≈ v_kepler - omega * (1-mu+r2) ≈ v_kepler - (1-mu) - r2
  // ただし簡単のため、回転座標系で地球周りの円軌道を直接設定:
  geo_cfg.vx0 = 0.0;
  geo_cfg.vy0 = kGeoVelocity;
  geo_cfg.vz0 = 0.0;
  // タイムステップと積分時間
  // GEOのケプラー周期: T = 2*pi*sqrt(r^3/mu) ≈ 2*pi*sqrt((2.8e-4)^3 / 3e-6) ≈ 0.017 (無次元)
  double geo_kepler_period = 2.0 * std::numbers::pi * std::sqrt(
      kGeoRadius_nd * kGeoRadius_nd * kGeoRadius_nd / kMU);
  std::cout << "  GEO Kepler period (nd): " << geo_kepler_period << std::endl;
  // 十分小さいタイムステップ: 周期の1/1000程度
  geo_cfg.timestep = geo_kepler_period / 1000.0;
  // 1年 = 2*pi の積分時間
  geo_cfg.total_time = 2.0 * std::numbers::pi;
  geo_cfg.output_interval = 100;  // 100ステップごとに出力

  // LEO設定
  SaliTimeSeriesConfig leo_cfg;
  leo_cfg.label = "LEO";
  leo_cfg.x0 = kEarthX + kLeoRadius_nd;
  leo_cfg.y0 = 0.0;
  leo_cfg.z0 = 0.0;
  leo_cfg.vx0 = 0.0;
  leo_cfg.vy0 = kLeoVelocity;
  leo_cfg.vz0 = 0.0;
  double leo_kepler_period = 2.0 * std::numbers::pi * std::sqrt(
      kLeoRadius_nd * kLeoRadius_nd * kLeoRadius_nd / kMU);
  std::cout << "  LEO Kepler period (nd): " << leo_kepler_period << std::endl;
  leo_cfg.timestep = leo_kepler_period / 1000.0;
  leo_cfg.total_time = 2.0 * std::numbers::pi;
  leo_cfg.output_interval = 1000;  // LEOは周期が短いので多めに間引き

  // ============ 補足テスト: ヘシアンスケール ============
  std::cout << "\n================================================================" << std::endl;
  std::cout << "  Hessian vs Coriolis Scale Report" << std::endl;
  std::cout << "================================================================" << std::endl;

  std::string hessian_report_path = output_dir + "/hessian_scale_report.txt";
  std::ofstream hessian_ofs(hessian_report_path);

  // GEO位置
  ReportHessianScale(kMU, kEarthX + kGeoRadius_nd, 0.0, 0.0, "GEO", std::cout);
  ReportHessianScale(kMU, kEarthX + kGeoRadius_nd, 0.0, 0.0, "GEO", hessian_ofs);

  // LEO位置
  ReportHessianScale(kMU, kEarthX + kLeoRadius_nd, 0.0, 0.0, "LEO (ISS)", std::cout);
  ReportHessianScale(kMU, kEarthX + kLeoRadius_nd, 0.0, 0.0, "LEO (ISS)", hessian_ofs);

  // L1近傍（研究対象の領域）
  double r_hill = std::pow(kMU / 3.0, 1.0 / 3.0);
  ReportHessianScale(kMU, kEarthX - r_hill, 0.0, 0.0, "L1 vicinity (Hill sphere)", std::cout);
  ReportHessianScale(kMU, kEarthX - r_hill, 0.0, 0.0, "L1 vicinity (Hill sphere)", hessian_ofs);

  hessian_ofs.close();
  std::cout << "  Report: " << hessian_report_path << std::endl;

  // ============ テスト実行 ============
  // GEOとLEO両方でテストを実行
  std::vector<SaliTimeSeriesConfig*> configs = {&geo_cfg, &leo_cfg};

  for (auto* cfg_ptr : configs) {
    auto& cfg = *cfg_ptr;
    std::string prefix = output_dir + "/" + cfg.label;

    std::cout << "\n================================================================" << std::endl;
    std::cout << "  Running tests for " << cfg.label << std::endl;
    std::cout << "================================================================" << std::endl;

    // テスト1: 二体問題
    std::cout << "\n  [Test 1] Two-body SALI (control)..." << std::endl;
    auto t1_start = std::chrono::steady_clock::now();
    RunTwoBodySALI(cfg, kMU, prefix + "_twobody.csv");
    auto t1_end = std::chrono::steady_clock::now();
    std::cout << "  Elapsed: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1_end - t1_start).count()
              << " ms" << std::endl;

    // テスト2: CR3BPシンプレクティック
    std::cout << "\n  [Test 2] CR3BP Symplectic SALI..." << std::endl;
    auto t2_start = std::chrono::steady_clock::now();
    RunCR3BPSymplecticSALI(cfg, kMU, prefix + "_crtbp_symplectic.csv");
    auto t2_end = std::chrono::steady_clock::now();
    std::cout << "  Elapsed: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2_end - t2_start).count()
              << " ms" << std::endl;

    // テスト3: CR3BP RK
    std::cout << "\n  [Test 3] CR3BP RK4 SALI..." << std::endl;
    auto t3_start = std::chrono::steady_clock::now();
    RunCR3BPRKSALI(cfg, kMU, prefix + "_crtbp_rk.csv", astro_params);
    auto t3_end = std::chrono::steady_clock::now();
    std::cout << "  Elapsed: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t3_end - t3_start).count()
              << " ms" << std::endl;
  }

  std::cout << "\n================================================================" << std::endl;
  std::cout << "  All tests completed." << std::endl;
  std::cout << "  Output directory: " << output_dir << std::endl;
  std::cout << "================================================================" << std::endl;

  return 0;
}
