/**
 * @file continuation.hpp
 * @brief 数値的継続法ライブラリ
 * @details 周期軌道族の系統的探索のための擬似弧長継続法・自然パラメータ継続法
 *          研究計画書 Step 1-2, 1-3 に対応
 * @date 2026-01-01
 */

#ifndef CONTINUATION_HPP
#define CONTINUATION_HPP

#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lyapunov_initial.hpp"
#include "periodic_orbit.hpp"
#include "periodic_orbit_impl.hpp"
#include "rtbp.hpp"

namespace continuation {

using namespace my_type;
using namespace periodic_orbit;
using namespace lyapunov;

// ---------------------------------------------------------------------------
// 列挙型定義
// ---------------------------------------------------------------------------

/**
 * @brief 継続法の種類
 */
enum class ContinuationMethod {
  NATURAL_PARAMETER,  ///< 自然パラメータ継続法（単純だが折り返し点に弱い）
  PSEUDO_ARCLENGTH,   ///< 擬似弧長継続法（折り返し点も追跡可能）
};

/**
 * @brief 分岐の種類
 */
enum class BifurcationType {
  PITCHFORK,        ///< ピッチフォーク分岐（固有値+1通過、Lyapunov→Halo）
  PERIOD_DOUBLING,  ///< 周期倍分岐（固有値-1通過）
  NEIMARK_SACKER,   ///< ネイマーク・サッカー分岐（複素固有値が単位円通過）
  UNKNOWN,
};

/**
 * @brief 継続の方向
 */
enum class ContinuationDirection {
  INCREASE_AMPLITUDE,  ///< 振幅を増加させる方向
  DECREASE_AMPLITUDE,  ///< 振幅を減少させる方向
  INCREASE_JACOBI,     ///< ヤコビ積分を増加させる方向
  DECREASE_JACOBI,     ///< ヤコビ積分を減少させる方向
};

// ---------------------------------------------------------------------------
// 構造体定義
// ---------------------------------------------------------------------------

/**
 * @brief 継続法の設定
 */
template <typename ScalarType>
struct ContinuationConfig {
  // 継続法の種類
  ContinuationMethod method = ContinuationMethod::PSEUDO_ARCLENGTH;

  // ステップ制御
  ScalarType initial_step_size = 0.001;  ///< 初期ステップサイズ
  ScalarType min_step_size = 1e-6;       ///< 最小ステップサイズ
  ScalarType max_step_size = 0.1;        ///< 最大ステップサイズ
  int max_steps = 2000;                  ///< 最大継続ステップ数

  // ヤコビ積分の範囲
  ScalarType jacobi_min = 2.85;  ///< ヤコビ積分下限
  ScalarType jacobi_max = 3.05;  ///< ヤコビ積分上限

  // Newton法の設定
  int newton_max_iterations = 50;
  ScalarType newton_tolerance = 1e-12;
  ScalarType max_integration_time = 500.0;
  ScalarType integration_timestep = 0.0001;

  // 分岐検出の設定
  bool detect_bifurcations = true;
  bool follow_branches = true;             ///< 分岐から派生族を追跡するか
  ScalarType eigenvalue_threshold = 1e-5;  ///< |λ - 1| < threshold で分岐点判定

  // ポアンカレ断面の設定
  int section_index = 1;           ///< 断面変数 (0:x, 1:y, 2:z)
  ScalarType section_value = 0.0;  ///< 断面の値

  // 開始軌道の設定
  OrbitFamilyType starting_family = OrbitFamilyType::LYAPUNOV_L1;
  ScalarType initial_amplitude = 0.001;

  // 並列化設定
  int num_threads = 0;  ///< 0 = 自動検出
};

/**
 * @brief 分岐点の情報
 */
template <typename ScalarType>
struct BifurcationPoint {
  int orbit_index;                               ///< 分岐が検出された軌道のインデックス
  ScalarType jacobi_integral;                    ///< 分岐点でのヤコビ積分
  BifurcationType type;                          ///< 分岐の種類
  std::complex<ScalarType> critical_eigenvalue;  ///< 分岐を引き起こした固有値
  std::array<ScalarType, 6> bifurcation_vector;  ///< 分岐方向ベクトル
  std::string description;                       ///< 分岐の説明
};

/**
 * @brief 周期軌道族
 */
template <typename ScalarType>
struct PeriodicOrbitFamily {
  std::string family_name;                                 ///< 族の名前 ("Lyapunov-L1"等)
  OrbitFamilyType family_type;                             ///< 族の種類
  std::vector<PeriodicOrbit<ScalarType>> orbits;           ///< 軌道のリスト
  std::vector<BifurcationPoint<ScalarType>> bifurcations;  ///< 検出された分岐点
  bool continuation_completed = false;                     ///< 継続が完了したか
  std::string termination_reason;                          ///< 終了理由
};

// ---------------------------------------------------------------------------
// 継続法クラス
// ---------------------------------------------------------------------------

/**
 * @brief 継続法エンジン
 */
template <typename ScalarType>
class OrbitContinuator {
 public:
  /**
   * @brief コンストラクタ
   * @param mu 質量比パラメータ
   * @param config 継続法の設定
   */
  OrbitContinuator(ScalarType mu, const ContinuationConfig<ScalarType>& config)
      : mu_(mu), config_(config) {
#ifdef _OPENMP
    if (config_.num_threads > 0) {
      omp_set_num_threads(config_.num_threads);
    }
#endif
  }

  /**
   * @brief 指定された起点から周期軌道族を継続
   * @param initial_orbit 継続の起点となる周期軌道
   * @return 周期軌道族
   */
  PeriodicOrbitFamily<ScalarType> ContinueFamily(const PeriodicOrbit<ScalarType>& initial_orbit) {
    PeriodicOrbitFamily<ScalarType> family;
    family.family_type = config_.starting_family;
    family.family_name = GetFamilyName(config_.starting_family);

    // 初期軌道を追加
    family.orbits.push_back(initial_orbit);

    switch (config_.method) {
      case ContinuationMethod::NATURAL_PARAMETER:
        ContinueNaturalParameter(&family);
        break;
      case ContinuationMethod::PSEUDO_ARCLENGTH:
        ContinuePseudoArclength(&family);
        break;
    }

    return family;
  }

  /**
   * @brief Lyapunov軌道族を自動生成して継続
   * @param point ラグランジュ点
   * @return 周期軌道族
   */
  PeriodicOrbitFamily<ScalarType> ContinueLyapunovFamily(LagrangePoint point) {
    // 初期条件を生成
    auto generator = std::make_unique<LyapunovOrbitGenerator<ScalarType>>(point, mu_);
    State<ScalarType> initial_guess = generator->GenerateInitialGuess(config_.initial_amplitude);

    std::cout << "=== Starting Lyapunov family continuation from "
              << GetFamilyName(generator->GetFamilyType()) << " ===\n";
    std::cout << "Initial amplitude: " << config_.initial_amplitude << "\n";
    std::cout << "Initial guess: (" << initial_guess.x << ", " << initial_guess.y << ", "
              << initial_guess.z << ", " << initial_guess.vx << ", " << initial_guess.vy << ", "
              << initial_guess.vz << ")\n";

    // Newton法で精密化
    NewtonConvergenceInfo<ScalarType> conv_info;
    PeriodicOrbit<ScalarType> initial_orbit;

    try {
      // Lyapunov軌道用の対称性ベース微分修正を使用
      initial_orbit = RefinePeriodicOrbitSymmetric(
          initial_guess, mu_,
          config_.newton_max_iterations,
          config_.newton_tolerance,
          config_.max_integration_time,
          config_.integration_timestep, &conv_info);

      if (!conv_info.converged) {
        throw std::runtime_error("Initial orbit refinement failed to converge");
      }

      std::cout << "Initial orbit refined successfully:\n";
      std::cout << "  Period: " << initial_orbit.period << "\n";
      std::cout << "  Jacobi constant: " << initial_orbit.jacobi_constant << "\n";

    } catch (const std::exception& e) {
      std::cerr << "Failed to refine initial orbit: " << e.what() << "\n";
      PeriodicOrbitFamily<ScalarType> empty_family;
      empty_family.family_type = generator->GetFamilyType();
      empty_family.family_name = GetFamilyName(generator->GetFamilyType());
      empty_family.termination_reason = std::string("Initial orbit refinement failed: ") + e.what();
      return empty_family;
    }

    // モノドロミー行列と安定性を計算
    initial_orbit.monodromy_matrix =
        ComputeMonodromyMatrix(initial_orbit, mu_, config_.integration_timestep);
    AnalyzeStability(&initial_orbit, 1.0);

    // 継続を実行
    PeriodicOrbitFamily<ScalarType> family = ContinueFamily(initial_orbit);
    family.family_type = generator->GetFamilyType();
    family.family_name = GetFamilyName(generator->GetFamilyType());

    return family;
  }

  /**
   * @brief 複数のラグランジュ点から並列に継続（OpenMP）
   * @param points ラグランジュ点のリスト
   * @return 複数の周期軌道族
   */
  std::vector<PeriodicOrbitFamily<ScalarType>> ContinueMultipleFamilies(
      const std::vector<LagrangePoint>& points) {
    std::vector<PeriodicOrbitFamily<ScalarType>> families(points.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
      families[i] = ContinueLyapunovFamily(points[i]);
    }

    return families;
  }

  /**
   * @brief ピッチフォーク分岐点からHalo族を生成して継続
   * @param bifurcation_orbit 分岐が検出されたLyapunov軌道
   * @param north_branch trueならNorth Halo, falseならSouth Halo
   * @return Halo族
   * 
   * @details Lyapunov軌道はxy平面内（z=0）なので、モノドロミー行列の固有ベクトルも
   *          z成分がほぼゼロになる。Halo軌道を生成するには、z方向に明示的に
   *          摂動を与える必要がある。摂動量はLyapunov軌道のx振幅に比例させる。
   */
  PeriodicOrbitFamily<ScalarType> ContinueHaloFromBifurcation(
      const PeriodicOrbit<ScalarType>& bifurcation_orbit, bool north_branch = true) {
    PeriodicOrbitFamily<ScalarType> family;

    // L2点の位置を推定
    ScalarType x_L2 = 1.0 + std::cbrt(mu_ / 3.0);
    
    // Lyapunov軌道のx振幅
    ScalarType x_amplitude = std::abs(bifurcation_orbit.initial_state.x - x_L2);
    
    // z方向の摂動量: x振幅に比例（Halo軌道の典型的なスケール）
    // 小さすぎるとNewtonが収束しにくい、大きすぎると別の軌道に飛ぶ
    ScalarType z_perturbation = x_amplitude * 0.1;  // x振幅の10%
    
    // North/Southの符号
    ScalarType z_sign = north_branch ? 1.0 : -1.0;

    // 分岐点から z方向に摂動した初期条件を生成
    State<ScalarType> initial_guess = bifurcation_orbit.initial_state;
    initial_guess.z = z_sign * z_perturbation;
    // vz は対称性から初期条件で0（y=0断面でのHalo対称性）

    std::cout << "=== Starting Halo family continuation from bifurcation point ===" << std::endl;
    std::cout << "  Branch: " << (north_branch ? "North" : "South") << std::endl;
    std::cout << "  Lyapunov x-amplitude: " << x_amplitude << std::endl;
    std::cout << "  Initial z perturbation: " << initial_guess.z << std::endl;
    std::cout << "  Initial state: (" << initial_guess.x << ", " << initial_guess.y 
              << ", " << initial_guess.z << ", " << initial_guess.vx 
              << ", " << initial_guess.vy << ", " << initial_guess.vz << ")" << std::endl;

    // Newton法で周期軌道として収束させる
    // Halo軌道はz対称性を使える: y=0で開始し、次のy=0交差まで積分
    NewtonConvergenceInfo<ScalarType> conv_info;
    PeriodicOrbit<ScalarType> initial_orbit;

    try {
      // Halo軌道用の3変数対称Newton法を使用
      // 変数: (x₀, z₀, vy₀), 目標: vx(T/2)=0, vz(T/2)=0
      initial_orbit = RefinePeriodicOrbitHalo(initial_guess, mu_,
                                               config_.newton_max_iterations,
                                               config_.newton_tolerance,
                                               config_.max_integration_time,
                                               config_.integration_timestep, &conv_info);

      if (!conv_info.converged) {
        throw std::runtime_error("Halo orbit refinement failed to converge");
      }

      std::cout << "  Halo orbit refined successfully!" << std::endl;
      std::cout << "  Period: " << initial_orbit.period << std::endl;
      std::cout << "  Final z: " << initial_orbit.initial_state.z << std::endl;
      std::cout << "  Jacobi: " << initial_orbit.jacobi_constant << std::endl;

    } catch (const std::exception& e) {
      std::cerr << "Failed to refine Halo orbit: " << e.what() << std::endl;
      family.family_type = north_branch ? OrbitFamilyType::HALO_L2_NORTH
                                        : OrbitFamilyType::HALO_L2_SOUTH;
      family.family_name = GetFamilyName(family.family_type);
      family.termination_reason = std::string("Initial Halo orbit refinement failed: ") + e.what();
      return family;
    }

    // モノドロミー行列と安定性を計算
    initial_orbit.monodromy_matrix =
        ComputeMonodromyMatrix(initial_orbit, mu_, config_.integration_timestep);
    AnalyzeStability(&initial_orbit, 1.0);

    // 継続を実行（Halo継続モードを有効化）
    is_halo_continuation_ = true;
    expected_z_sign_ = north_branch ? 1.0 : -1.0;  // z符号を記録
    family = ContinueFamily(initial_orbit);
    is_halo_continuation_ = false;  // 終了後にリセット
    expected_z_sign_ = 0.0;

    // 族の種類を設定
    family.family_type = north_branch ? OrbitFamilyType::HALO_L2_NORTH
                                      : OrbitFamilyType::HALO_L2_SOUTH;
    family.family_name = GetFamilyName(family.family_type);

    return family;
  }


 private:
  ScalarType mu_;
  ContinuationConfig<ScalarType> config_;
  ScalarType current_step_size_;
  bool is_halo_continuation_ = false;  ///< Halo継続モードフラグ
  ScalarType expected_z_sign_ = 0.0;   ///< 期待されるz符号 (+1: North, -1: South, 0: auto)

  /**
   * @brief 自然パラメータ継続法
   */
  void ContinueNaturalParameter(PeriodicOrbitFamily<ScalarType>* family) {
    current_step_size_ = config_.initial_step_size;
    int step_count = 0;

    while (step_count < config_.max_steps) {
      const auto& prev_orbit = family->orbits.back();

      // ヤコビ積分の範囲チェック
      if (prev_orbit.jacobi_constant < config_.jacobi_min ||
          prev_orbit.jacobi_constant > config_.jacobi_max) {
        family->termination_reason = "Jacobi integral out of range";
        family->continuation_completed = true;
        break;
      }

      // 予測ステップ: 振幅を増加
      State<ScalarType> predicted_state = PredictNextState(prev_orbit);

      // Newton法で補正
      try {
        NewtonConvergenceInfo<ScalarType> conv_info;
        PeriodicOrbit<ScalarType> next_orbit;
        if (is_halo_continuation_) {
          // Halo軌道用の3変数対称Newton法を使用
          next_orbit = RefinePeriodicOrbitHalo(
              predicted_state, mu_,
              config_.newton_max_iterations, config_.newton_tolerance,
              config_.max_integration_time, config_.integration_timestep, &conv_info,
              expected_z_sign_);
        } else {
          // Lyapunov軌道用の2変数対称Newton法を使用
          next_orbit = RefinePeriodicOrbitSymmetric(
              predicted_state, mu_,
              config_.newton_max_iterations, config_.newton_tolerance,
              config_.max_integration_time, config_.integration_timestep, &conv_info);
        }

        if (!conv_info.converged) {
          // ステップサイズを縮小して再試行
          current_step_size_ *= 0.5;
          if (current_step_size_ < config_.min_step_size) {
            family->termination_reason = "Step size became too small";
            family->continuation_completed = true;
            break;
          }
          continue;
        }

        // 安定性解析
        next_orbit.monodromy_matrix =
            ComputeMonodromyMatrix(next_orbit, mu_, config_.integration_timestep);
        AnalyzeStability(&next_orbit, 1.0);

        // 分岐検出
        if (config_.detect_bifurcations) {
          auto bifurcation = DetectBifurcation(prev_orbit, next_orbit, family->orbits.size());
          if (bifurcation.type != BifurcationType::UNKNOWN) {
            family->bifurcations.push_back(bifurcation);
            std::cout << "  Bifurcation detected at orbit " << family->orbits.size() << ": "
                      << bifurcation.description << "\n";
          }
        }

        family->orbits.push_back(next_orbit);
        step_count++;

        // 進捗表示
        if (step_count % 50 == 0) {
          std::cout << "  Step " << step_count << ": C = " << next_orbit.jacobi_constant
                    << ", T = " << next_orbit.period
                    << ", stable = " << (next_orbit.is_stable ? "yes" : "no") << "\n";
        }

        // ステップサイズの適応的調整
        if (conv_info.iterations < 5) {
          current_step_size_ = std::min(current_step_size_ * 1.2, config_.max_step_size);
        }

      } catch (const std::exception& e) {
        current_step_size_ *= 0.5;
        if (current_step_size_ < config_.min_step_size) {
          family->termination_reason = std::string("Continuation failed: ") + e.what();
          family->continuation_completed = true;
          break;
        }
      }
    }

    if (step_count >= config_.max_steps) {
      family->termination_reason = "Maximum steps reached";
      family->continuation_completed = true;
    }
  }

  /**
   * @brief 擬似弧長継続法
   */
  void ContinuePseudoArclength(PeriodicOrbitFamily<ScalarType>* family) {
    current_step_size_ = config_.initial_step_size;
    int step_count = 0;

    // 接線方向を推定するために少なくとも2つの軌道が必要
    if (family->orbits.size() < 2) {
      // まず自然パラメータ法で1ステップ進める
      ContinueNaturalParameterOneStep(family);
      if (family->orbits.size() < 2) {
        family->termination_reason = "Failed to get second orbit for pseudo-arclength";
        return;
      }
    }

    while (step_count < config_.max_steps) {
      const auto& prev_orbit = family->orbits.back();
      const auto& prev_prev_orbit = family->orbits[family->orbits.size() - 2];

      // ヤコビ積分の範囲チェック
      if (prev_orbit.jacobi_constant < config_.jacobi_min ||
          prev_orbit.jacobi_constant > config_.jacobi_max) {
        family->termination_reason = "Jacobi integral out of range";
        family->continuation_completed = true;
        break;
      }

      // 接線方向を計算（前の2軌道の差分）
      std::array<ScalarType, 5> tangent = ComputeTangentVector(prev_prev_orbit, prev_orbit);

      // 予測ステップ
      State<ScalarType> predicted_state = PredictAlongTangent(prev_orbit, tangent);

      // Newton法で補正（擬似弧長条件付き）
      try {
        NewtonConvergenceInfo<ScalarType> conv_info;
        PeriodicOrbit<ScalarType> next_orbit;
        if (is_halo_continuation_) {
          // Halo軌道用の3変数対称Newton法を使用
          next_orbit = RefinePeriodicOrbitHalo(
              predicted_state, mu_,
              config_.newton_max_iterations, config_.newton_tolerance,
              config_.max_integration_time, config_.integration_timestep, &conv_info,
              expected_z_sign_);
        } else {
          // Lyapunov軌道用の2変数対称Newton法を使用
          next_orbit = RefinePeriodicOrbitSymmetric(
              predicted_state, mu_,
              config_.newton_max_iterations, config_.newton_tolerance,
              config_.max_integration_time, config_.integration_timestep, &conv_info);
        }

        if (!conv_info.converged) {
          current_step_size_ *= 0.5;
          if (current_step_size_ < config_.min_step_size) {
            family->termination_reason = "Step size became too small";
            family->continuation_completed = true;
            break;
          }
          continue;
        }

        // 安定性解析
        next_orbit.monodromy_matrix =
            ComputeMonodromyMatrix(next_orbit, mu_, config_.integration_timestep);
        AnalyzeStability(&next_orbit, 1.0);

        // 分岐検出
        if (config_.detect_bifurcations) {
          auto bifurcation = DetectBifurcation(prev_orbit, next_orbit, family->orbits.size());
          if (bifurcation.type != BifurcationType::UNKNOWN) {
            family->bifurcations.push_back(bifurcation);
            std::cout << "  Bifurcation detected at orbit " << family->orbits.size() << ": "
                      << bifurcation.description << "\n";
          }
        }

        family->orbits.push_back(next_orbit);
        step_count++;

        // 進捗表示
        if (step_count % 50 == 0) {
          std::cout << "  Step " << step_count << ": C = " << next_orbit.jacobi_constant
                    << ", T = " << next_orbit.period
                    << ", stable = " << (next_orbit.is_stable ? "yes" : "no") << "\n";
        }

        // ステップサイズの適応的調整
        if (conv_info.iterations < 5) {
          current_step_size_ = std::min(current_step_size_ * 1.2, config_.max_step_size);
        }

      } catch (const std::exception& e) {
        current_step_size_ *= 0.5;
        if (current_step_size_ < config_.min_step_size) {
          family->termination_reason = std::string("Continuation failed: ") + e.what();
          family->continuation_completed = true;
          break;
        }
      }
    }

    if (step_count >= config_.max_steps) {
      family->termination_reason = "Maximum steps reached";
      family->continuation_completed = true;
    }
  }

  /**
   * @brief 自然パラメータ法で1ステップだけ進める
   */
  void ContinueNaturalParameterOneStep(PeriodicOrbitFamily<ScalarType>* family) {
    const auto& prev_orbit = family->orbits.back();
    State<ScalarType> predicted_state = PredictNextState(prev_orbit);

    try {
      NewtonConvergenceInfo<ScalarType> conv_info;
      PeriodicOrbit<ScalarType> next_orbit;
      if (is_halo_continuation_) {
        // Halo軌道用の3変数対称Newton法を使用
        next_orbit = RefinePeriodicOrbitHalo(
            predicted_state, mu_,
            config_.newton_max_iterations, config_.newton_tolerance,
            config_.max_integration_time, config_.integration_timestep, &conv_info,
            expected_z_sign_);
      } else {
        // Lyapunov軌道用の2変数対称Newton法を使用
        next_orbit = RefinePeriodicOrbitSymmetric(
            predicted_state, mu_,
            config_.newton_max_iterations, config_.newton_tolerance,
            config_.max_integration_time, config_.integration_timestep, &conv_info);
      }

      if (conv_info.converged) {
        next_orbit.monodromy_matrix =
            ComputeMonodromyMatrix(next_orbit, mu_, config_.integration_timestep);
        AnalyzeStability(&next_orbit, 1.0);
        family->orbits.push_back(next_orbit);
      }
    } catch (...) {
      // 失敗した場合は何もしない
    }
  }

  /**
   * @brief 次の状態を予測（振幅増加方向）
   * @details Lyapunov軌道では振幅(x - x_L)に比例してvyもスケールする
   *          ステップサイズは振幅の相対増加率として解釈
   *          L1 (x < 1) と L2 (x > 1) の両方に対応
   */
  State<ScalarType> PredictNextState(const PeriodicOrbit<ScalarType>& orbit) {
    // Halo継続モードの場合はHalo用予測を使用
    if (is_halo_continuation_) {
      return PredictNextStateHalo(orbit);
    }
    
    State<ScalarType> predicted = orbit.initial_state;

    // L1/L2の自動判定: x < 1 なら L1、x > 1 なら L2
    LagrangePoint point;
    if (orbit.initial_state.x < 1.0) {
      point = LagrangePoint::L1;
    } else {
      point = LagrangePoint::L2;
    }
    
    // Lyapunov係数を計算
    auto coeff = ComputeLyapunovCoefficients(point, mu_);
    
    // 現在の振幅 = |x - x_L|
    ScalarType current_amplitude = std::abs(orbit.initial_state.x - coeff.x_L);
    
    // ステップサイズを相対増加率として使用（例：0.001 = 0.1%増加）
    // 最小値を0.001 (0.1%)に設定して小振幅軌道に対応
    ScalarType relative_step = std::max(0.001, std::min(0.5, current_step_size_ * 100.0));
    ScalarType scale_factor = 1.0 + relative_step;
    
    // 新しい振幅
    ScalarType new_amplitude = current_amplitude * scale_factor;
    
    // Lyapunov軌道の初期条件生成ロジックと同じ方法で予測
    // x₀ = x_L + A_x（L1の場合も + で良い。ComputeLyapunovCoefficientsがx_Lを正しく計算）
    // vy₀ = -κ ω A_x
    predicted.x = coeff.x_L + new_amplitude;
    predicted.vy = -coeff.kappa * coeff.omega_xy * new_amplitude;
    
    return predicted;
  }


  
  /**
   * @brief Halo軌道用の次の状態を予測（3D振幅増加方向）
   * @details Halo軌道では x, z, vy を同じ比率でスケールする
   */
  State<ScalarType> PredictNextStateHalo(const PeriodicOrbit<ScalarType>& orbit) {
    State<ScalarType> predicted = orbit.initial_state;

    // L2点の位置を推定
    ScalarType x_L2 = 1.0 + std::cbrt(mu_ / 3.0);
    
    // ステップサイズを相対増加率として使用
    ScalarType relative_step = std::max(0.05, std::min(0.5, current_step_size_ * 100.0));
    ScalarType scale_factor = 1.0 + relative_step;
    
    // x, z, vy を同じ比率でスケール（Halo軌道の形状を保持）
    predicted.x = x_L2 + (orbit.initial_state.x - x_L2) * scale_factor;
    predicted.z = orbit.initial_state.z * scale_factor;
    predicted.vy = orbit.initial_state.vy * scale_factor;
    
    return predicted;
  }

  /**
   * @brief 接線方向ベクトルを計算（5次元: Halo対応）
   */
  std::array<ScalarType, 5> ComputeTangentVector(const PeriodicOrbit<ScalarType>& orbit1,
                                                 const PeriodicOrbit<ScalarType>& orbit2) {
    // (x0, z0, vy0, T, C) 空間での接線（Halo軌道対応）
    std::array<ScalarType, 5> tangent;
    tangent[0] = orbit2.initial_state.x - orbit1.initial_state.x;
    tangent[1] = orbit2.initial_state.z - orbit1.initial_state.z;  // z成分を追加
    tangent[2] = orbit2.initial_state.vy - orbit1.initial_state.vy;
    tangent[3] = orbit2.period - orbit1.period;
    tangent[4] = orbit2.jacobi_constant - orbit1.jacobi_constant;

    // 正規化
    ScalarType norm = 0.0;
    for (const auto& t : tangent) norm += t * t;
    norm = std::sqrt(norm);
    if (norm > 1e-14) {
      for (auto& t : tangent) t /= norm;
    }

    return tangent;
  }

  /**
   * @brief 接線方向に沿って予測（5次元: Halo対応）
   */
  State<ScalarType> PredictAlongTangent(const PeriodicOrbit<ScalarType>& orbit,
                                        const std::array<ScalarType, 5>& tangent) {
    State<ScalarType> predicted = orbit.initial_state;
    predicted.x += current_step_size_ * tangent[0];
    predicted.z += current_step_size_ * tangent[1];  // z成分を追加
    predicted.vy += current_step_size_ * tangent[2];
    return predicted;
  }

  /**
   * @brief 分岐点を検出
   */
  BifurcationPoint<ScalarType> DetectBifurcation(const PeriodicOrbit<ScalarType>& prev_orbit,
                                                 const PeriodicOrbit<ScalarType>& curr_orbit,
                                                 int orbit_index) {
    BifurcationPoint<ScalarType> bif;
    bif.orbit_index = static_cast<int>(orbit_index);
    bif.jacobi_integral = curr_orbit.jacobi_constant;
    bif.type = BifurcationType::UNKNOWN;

    if (!prev_orbit.stability_computed || !curr_orbit.stability_computed) {
      return bif;
    }

    // 固有値の変化を調べる
    for (size_t i = 0; i < prev_orbit.eigenvalues.size() && i < curr_orbit.eigenvalues.size();
         ++i) {
      const auto& prev_ev = prev_orbit.eigenvalues[i];
      const auto& curr_ev = curr_orbit.eigenvalues[i];

      // +1 通過（ピッチフォーク分岐）
      ScalarType prev_dist_to_1 = std::abs(prev_ev - std::complex<ScalarType>(1.0, 0.0));
      ScalarType curr_dist_to_1 = std::abs(curr_ev - std::complex<ScalarType>(1.0, 0.0));

      if ((prev_dist_to_1 > config_.eigenvalue_threshold &&
           curr_dist_to_1 < config_.eigenvalue_threshold) ||
          (prev_dist_to_1 < config_.eigenvalue_threshold &&
           curr_dist_to_1 > config_.eigenvalue_threshold)) {
        // 符号変化を確認
        ScalarType prev_real_minus_1 = prev_ev.real() - 1.0;
        ScalarType curr_real_minus_1 = curr_ev.real() - 1.0;

        if (prev_real_minus_1 * curr_real_minus_1 < 0) {
          bif.type = BifurcationType::PITCHFORK;
          bif.critical_eigenvalue = curr_ev;
          bif.description = "Pitchfork bifurcation (eigenvalue crossed +1)";
          return bif;
        }
      }

      // -1 通過（周期倍分岐）
      ScalarType prev_dist_to_m1 = std::abs(prev_ev - std::complex<ScalarType>(-1.0, 0.0));
      ScalarType curr_dist_to_m1 = std::abs(curr_ev - std::complex<ScalarType>(-1.0, 0.0));

      if ((prev_dist_to_m1 > config_.eigenvalue_threshold &&
           curr_dist_to_m1 < config_.eigenvalue_threshold) ||
          (prev_dist_to_m1 < config_.eigenvalue_threshold &&
           curr_dist_to_m1 > config_.eigenvalue_threshold)) {
        ScalarType prev_real_plus_1 = prev_ev.real() + 1.0;
        ScalarType curr_real_plus_1 = curr_ev.real() + 1.0;

        if (prev_real_plus_1 * curr_real_plus_1 < 0) {
          bif.type = BifurcationType::PERIOD_DOUBLING;
          bif.critical_eigenvalue = curr_ev;
          bif.description = "Period-doubling bifurcation (eigenvalue crossed -1)";
          return bif;
        }
      }

      // 単位円通過（ネイマーク・サッカー分岐）
      if (std::abs(prev_ev.imag()) > config_.eigenvalue_threshold &&
          std::abs(curr_ev.imag()) > config_.eigenvalue_threshold) {
        ScalarType prev_abs = std::abs(prev_ev);
        ScalarType curr_abs = std::abs(curr_ev);

        if ((prev_abs - 1.0) * (curr_abs - 1.0) < 0) {
          bif.type = BifurcationType::NEIMARK_SACKER;
          bif.critical_eigenvalue = curr_ev;
          bif.description = "Neimark-Sacker bifurcation (complex eigenvalue crossed unit circle)";
          return bif;
        }
      }
    }

    return bif;
  }

  /**
   * @brief 族の名前を取得
   */
  static std::string GetFamilyName(OrbitFamilyType type) {
    switch (type) {
      case OrbitFamilyType::LYAPUNOV_L1:
        return "Lyapunov-L1";
      case OrbitFamilyType::LYAPUNOV_L2:
        return "Lyapunov-L2";
      case OrbitFamilyType::LYAPUNOV_L3:
        return "Lyapunov-L3";
      case OrbitFamilyType::HALO_L1_NORTH:
        return "Halo-L1-North";
      case OrbitFamilyType::HALO_L1_SOUTH:
        return "Halo-L1-South";
      case OrbitFamilyType::HALO_L2_NORTH:
        return "Halo-L2-North";
      case OrbitFamilyType::HALO_L2_SOUTH:
        return "Halo-L2-South";
      case OrbitFamilyType::DRO:
        return "DRO";
      case OrbitFamilyType::RESONANT:
        return "Resonant";
      default:
        return "Unknown";
    }
  }
};

// ---------------------------------------------------------------------------
// 出力関数
// ---------------------------------------------------------------------------

/**
 * @brief 周期軌道族をCSVファイルに出力
 */
template <typename ScalarType>
void ExportFamilyToCsv(const PeriodicOrbitFamily<ScalarType>& family,
                       const std::string& output_path) {
  std::ofstream file(output_path);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open output file: " + output_path);
  }

  file << std::setprecision(16);
  file << "# Periodic Orbit Family: " << family.family_name << "\n";
  file << "# Columns: index,x0,y0,z0,vx0,vy0,vz0,period,jacobi,stable,stability_index\n";

  for (size_t i = 0; i < family.orbits.size(); ++i) {
    const auto& orbit = family.orbits[i];
    file << i << "," << orbit.initial_state.x << "," << orbit.initial_state.y << ","
         << orbit.initial_state.z << "," << orbit.initial_state.vx << "," << orbit.initial_state.vy
         << "," << orbit.initial_state.vz << "," << orbit.period << "," << orbit.jacobi_constant
         << "," << (orbit.is_stable ? 1 : 0) << "," << orbit.stability_index << "\n";
  }

  file.close();
}

/**
 * @brief 分岐点情報をCSVファイルに出力
 */
template <typename ScalarType>
void ExportBifurcationsToCsv(const PeriodicOrbitFamily<ScalarType>& family,
                             const std::string& output_path) {
  std::ofstream file(output_path);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open output file: " + output_path);
  }

  file << std::setprecision(16);
  file << "# Bifurcation Points for Family: " << family.family_name << "\n";
  file << "# Columns: orbit_index,jacobi,type,eigenvalue_real,eigenvalue_imag,description\n";

  for (const auto& bif : family.bifurcations) {
    std::string type_str;
    switch (bif.type) {
      case BifurcationType::PITCHFORK:
        type_str = "pitchfork";
        break;
      case BifurcationType::PERIOD_DOUBLING:
        type_str = "period_doubling";
        break;
      case BifurcationType::NEIMARK_SACKER:
        type_str = "neimark_sacker";
        break;
      default:
        type_str = "unknown";
    }

    file << bif.orbit_index << "," << bif.jacobi_integral << "," << type_str << ","
         << bif.critical_eigenvalue.real() << "," << bif.critical_eigenvalue.imag() << ","
         << "\"" << bif.description << "\"\n";
  }

  file.close();
}

/**
 * @brief GnuPlotスクリプトを生成
 */
template <typename ScalarType>
void GenerateFamilyGnuplotScript(const PeriodicOrbitFamily<ScalarType>& family,
                                 const std::string& output_dir, const std::string& data_filename) {
  std::string script_path = output_dir + "/plot_family.gp";
  std::ofstream file(script_path);

  file << "# Gnuplot script for " << family.family_name << "\n\n";

  file << "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n";
  file << "set output '" << output_dir << "/family_plot.png'\n\n";

  file << "set multiplot layout 2,2 title '" << family.family_name << " Family'\n\n";

  // Plot 1: Jacobi vs Period
  file << "set xlabel 'Jacobi Integral C'\n";
  file << "set ylabel 'Period T'\n";
  file << "plot '" << data_filename << "' using 9:8 with lines title 'C vs T', \\\n";
  file << "     '' using 9:8:(($10==1)?1:2) with points pt 7 ps 0.5 palette notitle\n\n";

  // Plot 2: Jacobi vs x0
  file << "set xlabel 'Jacobi Integral C'\n";
  file << "set ylabel 'x_0'\n";
  file << "plot '" << data_filename << "' using 9:2 with lines title 'C vs x_0'\n\n";

  // Plot 3: x0 vs vy0
  file << "set xlabel 'x_0'\n";
  file << "set ylabel 'v_{y,0}'\n";
  file << "plot '" << data_filename << "' using 2:6 with lines title 'x_0 vs v_{y,0}'\n\n";

  // Plot 4: Stability index
  file << "set xlabel 'Orbit Index'\n";
  file << "set ylabel 'Stability Index'\n";
  file << "set logscale y\n";
  file << "plot '" << data_filename << "' using 1:11 with lines title 'Stability Index'\n";

  file << "\nunset multiplot\n";
  file.close();
}

}  // namespace continuation

#endif  // CONTINUATION_HPP
