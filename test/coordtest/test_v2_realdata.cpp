/**
 * @file test_v2_realdata.cpp
 * @brief OBT_2_Earth.txt の実データを使った座標変換の健全性テスト
 *
 * テスト内容:
 * - 実際の天体暦データを読み込み
 * - ConvertInertial2RotatingV2 で座標変換を実行
 * - 地球が (1-μ, 0, 0) 付近に変換されることを確認
 * - ヤコビ積分が妥当な値であることを確認
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <rtbp.hpp>
#include <sstream>
#include <string>
#include <vector>

using namespace my_type;
using namespace crtbp;
using namespace param;

/**
 * @brief OBT_2_Earth.txt の1行をパースしてデータを取得
 *
 * フォーマット (13列):
 * JD  ast_x ast_y ast_z  ast_vx ast_vy ast_vz  earth_x earth_y earth_z  earth_vx earth_vy earth_vz
 *
 * @note 小惑星速度: AU/T* (無次元時間系)
 * @note 地球速度: AU/Day (JPL DE ephemeris)
 */
struct EphemerisData {
  double jd;            // ユリウス日
  State<double> ast;    // 小惑星状態量 (位置: AU, 速度: AU/T*)
  State<double> earth;  // 地球状態量 (位置: AU, 速度: AU/Day)
};

bool ParseLine(const std::string& line, EphemerisData* data) {
  std::istringstream iss(line);
  double values[13];
  for (int i = 0; i < 13; ++i) {
    if (!(iss >> values[i])) {
      return false;
    }
  }
  data->jd = values[0];
  data->ast = State<double>(values[1], values[2], values[3], values[4], values[5], values[6]);
  data->earth = State<double>(values[7], values[8], values[9], values[10], values[11], values[12]);
  return true;
}

/**
 * @brief データファイルから複数行を読み込む
 */
std::vector<EphemerisData> LoadEphemeris(const std::string& filename, int max_lines = 100) {
  std::vector<EphemerisData> data;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open file: " << filename << std::endl;
    return data;
  }

  std::string line;
  int count = 0;
  while (std::getline(file, line) && count < max_lines) {
    EphemerisData ephemeris;
    if (ParseLine(line, &ephemeris)) {
      data.push_back(ephemeris);
      count++;
    }
  }
  return data;
}

int main(int argc, char* argv[]) {
  std::cout << std::setprecision(10);
  std::cout << "================================================================" << std::endl;
  std::cout << "   Real Data Coordinate Transformation Sanity Test" << std::endl;
  std::cout << "================================================================\n" << std::endl;

  // データファイルパス
  std::string datafile = "test/data/eph_cpp_result/OBT_RESULT_PHAs_GR4_test/OBT_2_Earth.txt";
  if (argc > 1) {
    datafile = argv[1];
  }

  std::cout << "[Loading ephemeris data from: " << datafile << "]" << std::endl;

  // 天文定数
  AstroConstants<double> params;
  params.au = 1.495978707e11;        // [m]
  params.gm_sun = 1.32712440018e20;  // [m^3/s^2]
  params.gm_earth = 3.986004418e14;  // [m^3/s^2]

  const double gm_total = params.gm_sun + params.gm_earth;
  const double mu = params.gm_earth / gm_total;

  std::cout << "  mu = " << mu << std::endl;

  // データ読み込み
  std::vector<EphemerisData> ephemeris = LoadEphemeris(datafile, 20);
  if (ephemeris.empty()) {
    std::cerr << "Error: No data loaded." << std::endl;
    return 1;
  }

  std::cout << "  Loaded " << ephemeris.size() << " data points.\n" << std::endl;

  // 統計用
  double max_earth_pos_err = 0.0;
  double max_earth_vel = 0.0;
  double min_jacobi = 1e10;
  double max_jacobi = -1e10;
  int test_count = 0;
  int pass_count = 0;

  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "  Testing coordinate transformations..." << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;

  for (size_t i = 0; i < ephemeris.size(); ++i) {
    const EphemerisData& eph = ephemeris[i];

    // 座標変換実行
    State<double> ast_rot = ConvertInertial2RotatingV2(eph.ast, eph.earth, params);

    // 地球も変換（検証用）
    // 地球速度を AU/T* に変換してから渡す
    constexpr double k_AU_DAY_TO_ND = 365.25 / (2.0 * std::numbers::pi);
    State<double> earth_nd = eph.earth;
    earth_nd.vx = eph.earth.vx * k_AU_DAY_TO_ND;
    earth_nd.vy = eph.earth.vy * k_AU_DAY_TO_ND;
    earth_nd.vz = eph.earth.vz * k_AU_DAY_TO_ND;

    State<double> earth_rot = ConvertInertial2RotatingV2(earth_nd, eph.earth, params);

    // 地球位置の誤差 (期待値: (1-mu, 0, 0))
    double earth_pos_err = std::sqrt(std::pow(earth_rot.x - (1.0 - mu), 2) +
                                     std::pow(earth_rot.y, 2) + std::pow(earth_rot.z, 2));
    double earth_vel = std::sqrt(earth_rot.vx * earth_rot.vx + earth_rot.vy * earth_rot.vy +
                                 earth_rot.vz * earth_rot.vz);

    // ヤコビ積分
    double jacobi = calc_jacobi_integral(ast_rot, mu);

    // 統計更新
    max_earth_pos_err = std::max(max_earth_pos_err, earth_pos_err);
    max_earth_vel = std::max(max_earth_vel, earth_vel);
    min_jacobi = std::min(min_jacobi, jacobi);
    max_jacobi = std::max(max_jacobi, jacobi);
    test_count++;

    // 健全性チェック
    bool earth_pos_ok = earth_pos_err < 1e-10;
    bool jacobi_ok = (jacobi > 2.5 && jacobi < 3.5);  // 地球近傍の妥当な範囲

    if (earth_pos_ok && jacobi_ok) {
      pass_count++;
    }

    // 最初と最後、および代表的なデータを出力
    if (i < 3 || i == ephemeris.size() - 1 || i == ephemeris.size() / 2) {
      std::cout << "\n  [Sample " << i << "] JD = " << eph.jd << std::endl;
      std::cout << "    Input (J2000):" << std::endl;
      std::cout << "      Earth |r| = "
                << std::sqrt(eph.earth.x * eph.earth.x + eph.earth.y * eph.earth.y +
                             eph.earth.z * eph.earth.z)
                << " AU" << std::endl;
      std::cout << "          state =" << eph.earth.x << " " << eph.earth.y << " " << eph.earth.z
                << " " << eph.earth.vx << " " << eph.earth.vy << " " << eph.earth.vz << std::endl;
      std::cout << "      Asteroid |r| = "
                << std::sqrt(eph.ast.x * eph.ast.x + eph.ast.y * eph.ast.y + eph.ast.z * eph.ast.z)
                << " AU" << std::endl;
      std::cout << "          state =" << eph.ast.x << " " << eph.ast.y << " " << eph.ast.z << " "
                << eph.ast.vx << " " << eph.ast.vy << " " << eph.ast.vz << std::endl;
      std::cout << "    Output (Rotating):" << std::endl;
      std::cout << "      Earth pos: (" << earth_rot.x << ", " << earth_rot.y << ", " << earth_rot.z
                << ")" << std::endl;
      std::cout << "      Earth pos error: " << earth_pos_err << std::endl;
      std::cout << "      Earth |v|: " << earth_vel << std::endl;
      std::cout << "      Earth   v: " << earth_rot.vx << " " << earth_rot.vy << " " << earth_rot.vz
                << std::endl;
      std::cout << "      Asteroid pos: (" << ast_rot.x << ", " << ast_rot.y << ", " << ast_rot.z
                << ")" << std::endl;
      std::cout << "      Asteroid v: " << ast_rot.vx << " " << ast_rot.vy << " " << ast_rot.vz
                << std::endl;
      std::cout << "      Jacobi integral: " << jacobi << std::endl;
      std::cout << "      Status: " << (earth_pos_ok && jacobi_ok ? "PASS" : "FAIL") << std::endl;
    }
  }

  // サマリー
  std::cout << "\n================================================================" << std::endl;
  std::cout << "   Test Summary" << std::endl;
  std::cout << "================================================================" << std::endl;
  std::cout << "  Tests run: " << test_count << std::endl;
  std::cout << "  Tests passed: " << pass_count << std::endl;
  std::cout << "  Pass rate: " << (100.0 * pass_count / test_count) << "%" << std::endl;
  std::cout << "\n  Statistics:" << std::endl;
  std::cout << "    Max Earth position error: " << max_earth_pos_err << std::endl;
  std::cout << "    Max Earth velocity in rotating frame: " << max_earth_vel << std::endl;
  std::cout << "    Jacobi integral range: [" << min_jacobi << ", " << max_jacobi << "]"
            << std::endl;

  // 全体の健全性判定
  bool overall_pass = (pass_count == test_count) && (max_earth_pos_err < 1e-8) &&
                      (min_jacobi > 2.0 && max_jacobi < 4.0);

  std::cout << "\n  Overall Result: " << (overall_pass ? "PASS" : "FAIL") << std::endl;
  std::cout << "================================================================" << std::endl;

  return overall_pass ? 0 : 1;
}
