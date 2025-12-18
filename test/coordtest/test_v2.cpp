/**
 * @file test_v2.cpp
 * @brief ConvertInertial2RotatingV2 の数学的正しさを検証するテスト
 *
 * テストケース:
 * 1. 理想円軌道（e=0, i=0）
 * 2. 実際の地球軌道パラメータ（e=0.0167, i≈0）- 近日点
 * 3. 実際の地球軌道パラメータ - 遠日点
 * 4. 実際の地球軌道パラメータ - 中間位置（真近点角 = 90°）
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <rtbp.hpp>

using namespace my_type;
using namespace crtbp;
using namespace param;

/**
 * @brief 軌道要素から状態量を計算
 *
 * @param a 軌道長半径 [AU]
 * @param e 離心率
 * @param i 軌道傾斜角 [rad]
 * @param Omega 昇交点赤経 [rad]
 * @param omega 近点引数 [rad]
 * @param f 真近点角 [rad]
 * @param gm 重力定数 (G*M) [AU^3/T*^2] where T* = sqrt(AU^3/GM_sun)
 * @param state 出力状態量
 */
void OrbitalElements2State(double a, double e, double i, double Omega, double omega, double f,
                           double gm, State<double>* state) {
  // 軌道面内座標
  const double r = a * (1.0 - e * e) / (1.0 + e * std::cos(f));

  // 軌道面内の位置・速度
  const double x_orb = r * std::cos(f);
  const double y_orb = r * std::sin(f);

  // vis-viva equation から速度の大きさを計算
  // v^2 = GM * (2/r - 1/a)
  const double v_sq = gm * (2.0 / r - 1.0 / a);
  const double v = std::sqrt(v_sq);

  // 動径方向と接線方向の速度成分
  const double h = std::sqrt(gm * a * (1.0 - e * e));  // 比角運動量
  const double vr = gm * e * std::sin(f) / h;          // 動径方向速度
  const double vt = h / r;                             // 接線方向速度

  const double vx_orb = vr * std::cos(f) - vt * std::sin(f);
  const double vy_orb = vr * std::sin(f) + vt * std::cos(f);

  // 3次元回転（軌道面 -> 慣性系）
  // R = R_z(-Omega) * R_x(-i) * R_z(-omega)
  const double cos_O = std::cos(Omega);
  const double sin_O = std::sin(Omega);
  const double cos_i = std::cos(i);
  const double sin_i = std::sin(i);
  const double cos_w = std::cos(omega);
  const double sin_w = std::sin(omega);

  // 回転行列の成分
  const double r11 = cos_O * cos_w - sin_O * sin_w * cos_i;
  const double r12 = -cos_O * sin_w - sin_O * cos_w * cos_i;
  const double r21 = sin_O * cos_w + cos_O * sin_w * cos_i;
  const double r22 = -sin_O * sin_w + cos_O * cos_w * cos_i;
  const double r31 = sin_w * sin_i;
  const double r32 = cos_w * sin_i;

  // 位置変換
  state->x = r11 * x_orb + r12 * y_orb;
  state->y = r21 * x_orb + r22 * y_orb;
  state->z = r31 * x_orb + r32 * y_orb;

  // 速度変換
  state->vx = r11 * vx_orb + r12 * vy_orb;
  state->vy = r21 * vx_orb + r22 * vy_orb;
  state->vz = r31 * vx_orb + r32 * vy_orb;
}

/**
 * @brief テストケースを実行
 *
 * @note ConvertInertial2RotatingV2 の入力仕様:
 *       - 小惑星速度: [AU/T*] (無次元時間系)
 *       - 地球速度: [AU/Day] (JPL DE ephemeris)
 */
void RunTestCase(const std::string& name, const State<double>& earth_nd,
                 const State<double>& asteroid, const AstroConstants<double>& params, double mu) {
  // 変換係数: AU/T* → AU/Day
  // v[AU/Day] = v[AU/T*] / (365.25 / 2π)
  constexpr double k_ND_TO_AU_DAY = 2.0 * std::numbers::pi / 365.25;  // ≈ 0.0172

  // 地球速度を AU/Day に変換 (関数の入力仕様に合わせる)
  State<double> earth_input = earth_nd;
  earth_input.vx = earth_nd.vx * k_ND_TO_AU_DAY;
  earth_input.vy = earth_nd.vy * k_ND_TO_AU_DAY;
  earth_input.vz = earth_nd.vz * k_ND_TO_AU_DAY;

  std::cout << "\n----------------------------------------------------------------" << std::endl;
  std::cout << "  Test Case: " << name << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;

  // 入力表示
  double r_earth =
      std::sqrt(earth_nd.x * earth_nd.x + earth_nd.y * earth_nd.y + earth_nd.z * earth_nd.z);
  double v_earth_nd =
      std::sqrt(earth_nd.vx * earth_nd.vx + earth_nd.vy * earth_nd.vy + earth_nd.vz * earth_nd.vz);
  double v_earth_day = v_earth_nd * k_ND_TO_AU_DAY;
  std::cout << "  Earth |r| = " << r_earth << " AU" << std::endl;
  std::cout << "  Earth |v| = " << v_earth_nd << " AU/T* = " << v_earth_day << " AU/Day"
            << std::endl;
  std::cout << "  Earth pos: (" << earth_nd.x << ", " << earth_nd.y << ", " << earth_nd.z << ")"
            << std::endl;
  std::cout << "  Earth vel (ND): (" << earth_nd.vx << ", " << earth_nd.vy << ", " << earth_nd.vz
            << ")" << std::endl;

  // 地球の変換
  // ast_state (第1引数) には AU/T* の速度を持つ earth_nd を渡す
  // p2_state (第2引数) には AU/Day の速度を持つ earth_input を渡す
  State<double> earth_rot = ConvertInertial2RotatingV2(earth_nd, earth_input, params);

  std::cout << "\n  [Earth in Rotating Frame]" << std::endl;
  std::cout << "    Position: (" << earth_rot.x << ", " << earth_rot.y << ", " << earth_rot.z << ")"
            << std::endl;
  std::cout << "    Expected x: " << (1.0 - mu) << std::endl;

  double pos_err = std::sqrt(std::pow(earth_rot.x - (1.0 - mu), 2) + std::pow(earth_rot.y, 2) +
                             std::pow(earth_rot.z, 2));
  std::cout << "    Position error: " << pos_err << std::endl;

  double vel_mag = std::sqrt(earth_rot.vx * earth_rot.vx + earth_rot.vy * earth_rot.vy +
                             earth_rot.vz * earth_rot.vz);
  std::cout << "    Velocity: (" << earth_rot.vx << ", " << earth_rot.vy << ", " << earth_rot.vz
            << ")" << std::endl;
  std::cout << "    |v| (should be 0 for circular): " << vel_mag << std::endl;

  // 小惑星の変換
  State<double> ast_rot = ConvertInertial2RotatingV2(asteroid, earth_input, params);
  double jacobi = calc_jacobi_integral(ast_rot, mu);

  std::cout << "\n  [Asteroid in Rotating Frame]" << std::endl;
  std::cout << "    Position: (" << ast_rot.x << ", " << ast_rot.y << ", " << ast_rot.z << ")"
            << std::endl;
  std::cout << "    Velocity: (" << ast_rot.vx << ", " << ast_rot.vy << ", " << ast_rot.vz << ")"
            << std::endl;
  std::cout << "    Jacobi integral: " << jacobi << std::endl;

  // 判定
  bool pos_ok = pos_err < 1e-10;
  bool vel_ok_approx = true;  // 楕円軌道では速度が0にならない

  std::cout << "\n  [Result] Position test: " << (pos_ok ? "PASS" : "FAIL") << std::endl;
}

int main() {
  std::cout << std::setprecision(15);
  std::cout << "================================================================" << std::endl;
  std::cout << "   ConvertInertial2RotatingV2 Mathematical Verification Test" << std::endl;
  std::cout << "   (with realistic Earth orbital parameters)" << std::endl;
  std::cout << "================================================================\n" << std::endl;

  // 天文定数
  AstroConstants<double> params;
  params.au = 1.495978707e11;        // [m]
  params.gm_sun = 1.32712440018e20;  // [m^3/s^2]
  params.gm_earth = 3.986004418e14;  // [m^3/s^2]

  const double gm_total = params.gm_sun + params.gm_earth;
  const double mu = params.gm_earth / gm_total;

  // 無次元系での太陽重力定数
  // T* = sqrt(AU^3 / GM_sun) のとき、GM_sun = 4*pi^2 [AU^3/T*^2]
  const double gm_nd = 4.0 * std::numbers::pi * std::numbers::pi;  // ≈ 39.478

  std::cout << "[Parameters]" << std::endl;
  std::cout << "  mu (mass ratio) = " << mu << std::endl;
  std::cout << "  GM_sun (non-dim) = " << gm_nd << " [AU^3/T*^2]" << std::endl;

  // 地球の軌道要素（実際の値）
  const double a_earth = 1.00000261;  // 軌道長半径 [AU]
  const double e_earth = 0.01671123;  // 離心率
  const double i_earth = 0.00005;     // 軌道傾斜角 [rad] (黄道面基準でほぼ0)
  const double Omega_earth = 0.0;     // 昇交点赤経 [rad] (黄道面基準で0)
  const double omega_earth = 1.7968;  // 近点引数 [rad] ≈ 102.9°

  std::cout << "\n[Earth Orbital Elements]" << std::endl;
  std::cout << "  Semi-major axis a = " << a_earth << " AU" << std::endl;
  std::cout << "  Eccentricity e = " << e_earth << std::endl;
  std::cout << "  Inclination i = " << i_earth << " rad" << std::endl;
  std::cout << "  Arg. of perihelion = " << omega_earth * 180.0 / std::numbers::pi << " deg"
            << std::endl;

  // ========================================
  // Test Case 1: 理想円軌道（比較用）
  // ========================================
  {
    State<double> earth, asteroid;
    OrbitalElements2State(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, gm_nd, &earth);

    // 小惑星: 地球より少し外側
    OrbitalElements2State(1.02, 0.05, 0.01, 0.0, 0.0, 0.1, gm_nd, &asteroid);

    RunTestCase("Ideal Circular Orbit (e=0, i=0, f=0)", earth, asteroid, params, mu);
  }

  // ========================================
  // Test Case 2: 実際の地球軌道 - 近日点 (f=0)
  // ========================================
  {
    State<double> earth, asteroid;
    OrbitalElements2State(a_earth, e_earth, i_earth, Omega_earth, omega_earth, 0.0, gm_nd, &earth);

    // 小惑星: 地球から少しずれた位置
    double f_ast = 0.05;
    OrbitalElements2State(1.02, 0.05, 0.01, 0.0, 0.0, f_ast, gm_nd, &asteroid);

    RunTestCase("Realistic Earth at Perihelion (f=0)", earth, asteroid, params, mu);
  }

  // ========================================
  // Test Case 3: 実際の地球軌道 - 遠日点 (f=π)
  // ========================================
  {
    State<double> earth, asteroid;
    OrbitalElements2State(a_earth, e_earth, i_earth, Omega_earth, omega_earth, std::numbers::pi,
                          gm_nd, &earth);

    OrbitalElements2State(1.02, 0.05, 0.01, 0.0, 0.0, std::numbers::pi + 0.05, gm_nd, &asteroid);

    RunTestCase("Realistic Earth at Aphelion (f=pi)", earth, asteroid, params, mu);
  }

  // ========================================
  // Test Case 4: 実際の地球軌道 - 中間位置 (f=π/2)
  // ========================================
  {
    State<double> earth, asteroid;
    OrbitalElements2State(a_earth, e_earth, i_earth, Omega_earth, omega_earth,
                          std::numbers::pi / 2.0, gm_nd, &earth);

    OrbitalElements2State(1.02, 0.05, 0.01, 0.0, 0.0, std::numbers::pi / 2.0 + 0.05, gm_nd,
                          &asteroid);

    RunTestCase("Realistic Earth at f=90deg", earth, asteroid, params, mu);
  }

  // ========================================
  // Test Case 5: 傾斜角ありの地球軌道（仮想）
  // ========================================
  {
    State<double> earth, asteroid;
    const double i_test = 0.05;  // 約2.9度の傾斜
    OrbitalElements2State(a_earth, e_earth, i_test, 0.5, omega_earth, 0.3, gm_nd, &earth);

    OrbitalElements2State(1.02, 0.05, 0.02, 0.0, 0.0, 0.35, gm_nd, &asteroid);

    RunTestCase("Hypothetical Inclined Earth (i=2.9deg)", earth, asteroid, params, mu);
  }

  // ========================================
  // Test Case 6: 実データサンプル (OBT_2_Earth.txt 空行後)
  // ========================================
  {
    std::cout << "\n----------------------------------------------------------------" << std::endl;
    std::cout << "  Test Case: Real Data Sample (from OBT_2_Earth.txt after blank line)"
              << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;

    // データ: JD = 2468660.061541989
    // 小惑星: x y z vx vy vz (位置: AU, 速度: AU/T*)
    State<double> asteroid;
    asteroid.x = 6.550096686453792e-01;
    asteroid.y = 7.189035128391313e-01;
    asteroid.z = 2.339597200619310e-02;
    asteroid.vx = -3.869069813848611e-01;
    asteroid.vy = 1.100467346833496e+00;
    asteroid.vz = 1.248060297393993e-01;

    // 地球: x y z vx vy vz (位置: AU, 速度: AU/Day)
    State<double> earth_day;
    earth_day.x = 6.670157198760734e-01;
    earth_day.y = 7.320278924419075e-01;
    earth_day.z = -8.736685309734815e-05;
    earth_day.vx = -1.299452285447134e-02;
    earth_day.vy = 1.056617814344062e-02;
    earth_day.vz = 4.579643795894608e-03;

    // 地球速度を AU/T* に変換 (表示用)
    constexpr double k_AU_DAY_TO_ND = 365.25 / (2.0 * std::numbers::pi);
    State<double> earth_nd = earth_day;
    earth_nd.vx = earth_day.vx * k_AU_DAY_TO_ND;
    earth_nd.vy = earth_day.vy * k_AU_DAY_TO_ND;
    earth_nd.vz = earth_day.vz * k_AU_DAY_TO_ND;

    double r_earth = std::sqrt(earth_day.x * earth_day.x + earth_day.y * earth_day.y +
                               earth_day.z * earth_day.z);
    double v_earth_day = std::sqrt(earth_day.vx * earth_day.vx + earth_day.vy * earth_day.vy +
                                   earth_day.vz * earth_day.vz);
    std::cout << "  Earth |r| = " << r_earth << " AU" << std::endl;
    std::cout << "  Earth |v| = " << v_earth_day << " AU/Day" << std::endl;
    std::cout << "  Earth pos: (" << earth_day.x << ", " << earth_day.y << ", " << earth_day.z
              << ")" << std::endl;
    std::cout << "  Earth vel (AU/Day): (" << earth_day.vx << ", " << earth_day.vy << ", "
              << earth_day.vz << ")" << std::endl;

    double r_ast =
        std::sqrt(asteroid.x * asteroid.x + asteroid.y * asteroid.y + asteroid.z * asteroid.z);
    std::cout << "  Asteroid |r| = " << r_ast << " AU" << std::endl;
    std::cout << "  Asteroid pos: (" << asteroid.x << ", " << asteroid.y << ", " << asteroid.z
              << ")" << std::endl;
    std::cout << "  Asteroid vel (AU/T*): (" << asteroid.vx << ", " << asteroid.vy << ", "
              << asteroid.vz << ")" << std::endl;

    // 座標変換実行
    State<double> ast_rot = ConvertInertial2RotatingV2(asteroid, earth_day, params);

    // 地球も変換（検証用）
    State<double> earth_rot = ConvertInertial2RotatingV2(earth_nd, earth_day, params);

    std::cout << "\n  [Earth in Rotating Frame]" << std::endl;
    std::cout << "    Position: (" << earth_rot.x << ", " << earth_rot.y << ", " << earth_rot.z
              << ")" << std::endl;
    double earth_pos_err = std::sqrt(std::pow(earth_rot.x - (1.0 - mu), 2) +
                                     std::pow(earth_rot.y, 2) + std::pow(earth_rot.z, 2));
    std::cout << "    Position error: " << earth_pos_err << std::endl;
    double earth_vel = std::sqrt(earth_rot.vx * earth_rot.vx + earth_rot.vy * earth_rot.vy +
                                 earth_rot.vz * earth_rot.vz);
    std::cout << "    Velocity: (" << earth_rot.vx << ", " << earth_rot.vy << ", " << earth_rot.vz
              << ")" << std::endl;
    std::cout << "    |v| (should be ~0): " << earth_vel << std::endl;

    double jacobi = calc_jacobi_integral(ast_rot, mu);
    std::cout << "\n  [Asteroid in Rotating Frame]" << std::endl;
    std::cout << "    Position: (" << ast_rot.x << ", " << ast_rot.y << ", " << ast_rot.z << ")"
              << std::endl;
    std::cout << "    Velocity: (" << ast_rot.vx << ", " << ast_rot.vy << ", " << ast_rot.vz << ")"
              << std::endl;
    std::cout << "    Jacobi integral: " << jacobi << std::endl;

    bool pos_ok = earth_pos_err < 1e-10;
    std::cout << "\n  [Result] Position test: " << (pos_ok ? "PASS" : "FAIL") << std::endl;
  }

  std::cout << "\n================================================================" << std::endl;
  std::cout << "   All tests completed" << std::endl;
  std::cout << "================================================================" << std::endl;

  return 0;
}
