/**
 * @file main.cpp
 * @brief ヤコビ積分計算アプリケーション
 *
 * OBT_\d+Earth.txt ファイルの2番目のセクション（2連続空行で区切り）の
 * 最初の行を読み込み、CR3BP座標変換後にヤコビ積分を計算してCSV出力する。
 *
 * @note TrajectorySaliの実装を参考に作成
 */

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <rtbp.hpp>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector>

namespace fs = std::filesystem;

using namespace my_type;
using namespace crtbp;
using namespace param;
using namespace utils;

/**
 * @brief 軌道データ1行分
 */
struct OrbitDataRow {
  double time_j2000;       ///< J2000時刻
  State<double> asteroid;  ///< 小惑星状態量
  State<double> earth;     ///< 地球状態量
};

/**
 * @brief 2連続空行で区切られた2番目のセクションの最初の行を読み込む
 * @param[in] filepath ファイルパス
 * @param[out] row 読み込んだデータ
 * @return 成功時true
 */
bool LoadSecondSectionFirstLine(const std::string& filepath, OrbitDataRow* row) {
  std::ifstream ifs(filepath);
  if (!ifs) {
    std::cerr << "<> !err! Cannot open orbit data file: " << filepath << std::endl;
    return false;
  }

  std::string line;
  int consecutive_empty = 0;
  bool in_second_section = false;

  while (std::getline(ifs, line)) {
    // 空行判定
    bool is_empty = line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos;

    if (is_empty) {
      consecutive_empty++;
      if (consecutive_empty >= 2) {
        in_second_section = true;
      }
    } else {
      if (in_second_section) {
        // 2番目のセクションの最初の有効行
        std::stringstream ss(line);
        ss >> row->time_j2000;
        ss >> row->asteroid.x >> row->asteroid.y >> row->asteroid.z;
        ss >> row->asteroid.vx >> row->asteroid.vy >> row->asteroid.vz;
        ss >> row->earth.x >> row->earth.y >> row->earth.z;
        ss >> row->earth.vx >> row->earth.vy >> row->earth.vz;

        if (!ss.fail()) {
          return true;
        } else {
          std::cerr << "<> !err! Failed to parse line: " << line << std::endl;
          return false;
        }
      }
      consecutive_empty = 0;
    }
  }

  std::cerr << "<> !err! Second section not found in: " << filepath << std::endl;
  return false;
}

/**
 * @brief DATA_DIRからEarth軌道ファイルを発見
 * @param[in] data_dir データディレクトリ
 * @return ファイルパスのリスト
 */
std::vector<std::string> DiscoverOrbitFiles(const std::string& data_dir) {
  std::vector<std::string> files;
  if (!fs::exists(data_dir)) {
    std::cerr << "<> !err! DATA_DIR does not exist: " << data_dir << std::endl;
    return files;
  }

  // OBT_\d+Earth.txt または OBT_\d+_Earth.txt の両方にマッチ
  const std::regex pattern(R"(OBT_\d+_?Earth\.txt)");
  for (const auto& entry : fs::directory_iterator(data_dir)) {
    if (!entry.is_regular_file()) continue;
    std::string filename = entry.path().filename().string();
    if (std::regex_match(filename, pattern)) {
      files.push_back(entry.path().string());
    }
  }

  std::sort(files.begin(), files.end());
  return files;
}

int main(int argc, char* argv[]) {
  // ヘッダー出力
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>            Jacobi Integral Calculator ver1.1" << std::endl;
  std::cout << "<>            (Output sorted by Jacobi integral)" << std::endl;
  std::cout << "<>----------------------------------------------------------------\n\n"
            << std::endl;

  // コマンドライン引数チェック
  if (argc < 2) {
    std::cerr << "<> Usage: " << argv[0] << " <data_directory>" << std::endl;
    std::cerr << "<>        Processes OBT_*Earth.txt files in the specified directory."
              << std::endl;
    return -1;
  }

  std::string data_dir = argv[1];
  std::cout << "<>    Data directory: " << data_dir << std::endl;

  // 天文定数読み込み (CONFIG_DIRマクロを使用)
  const std::string kAstroParamFile = std::string(CONFIG_DIR) + "/astro_param/astro_param.txt";
  AstroConstants<double> astro_params = loadConstants<double>(kAstroParamFile);
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  std::cout << "<>    mu parameter: " << std::setprecision(15) << kMU << std::endl;
  std::cout << "<>" << std::endl;

  // 軌道データファイル発見
  std::vector<std::string> orbit_files = DiscoverOrbitFiles(data_dir);
  if (orbit_files.empty()) {
    std::cerr << "<> !err! No orbit files found in: " << data_dir << std::endl;
    return -1;
  }

  std::cout << "<>    Found " << orbit_files.size() << " orbit files" << std::endl;
  std::cout << "<>" << std::endl;

  // 結果を格納する構造体
  struct JacobiResult {
    std::string filename;
    double jacobi;
    double time_j2000;
    State<double> state_rot;
  };

  std::vector<JacobiResult> results;
  results.reserve(orbit_files.size());

  int success_count = 0;
  int fail_count = 0;

  // 各ファイルを処理
  for (const auto& orbit_file : orbit_files) {
    std::string filename = fs::path(orbit_file).filename().string();

    OrbitDataRow row;
    if (!LoadSecondSectionFirstLine(orbit_file, &row)) {
      std::cerr << "<> !warn! Skipping: " << filename << std::endl;
      fail_count++;
      continue;
    }

    // J2000 → CR3BP回転座標変換
    State<double> state_rot = ConvertInertial2RotatingV2(row.asteroid, row.earth, astro_params);

    // ヤコビ積分計算
    double jacobi = calc_jacobi_integral(state_rot, kMU);

    // 結果を保存
    results.push_back({filename, jacobi, row.time_j2000, state_rot});

    success_count++;

    // 進捗表示（100ファイルごと）
    if (success_count % 100 == 0) {
      std::cout << "\r<>    Processed: " << success_count << " files" << std::flush;
    }
  }

  std::cout << "\r<>    Processed: " << success_count << " files" << std::endl;

  // ヤコビ積分でソート（昇順）
  std::sort(results.begin(), results.end(),
            [](const JacobiResult& a, const JacobiResult& b) { return a.jacobi < b.jacobi; });

  std::cout << "<>    Sorting by Jacobi integral..." << std::endl;

  // 出力CSVファイル
  std::string output_csv = data_dir + "/jacobi_integral_results.csv";
  std::ofstream ofs(output_csv);
  if (!ofs) {
    std::cerr << "<> !err! Cannot create output file: " << output_csv << std::endl;
    return -1;
  }

  ofs << std::setprecision(15) << std::fixed;
  ofs << "# Jacobi Integral Calculation Results (Sorted by Jacobi Integral)\n";
  ofs << "# Data directory: " << data_dir << "\n";
  ofs << "# mu = " << kMU << "\n";
  ofs << "filename,jacobi_integral,time_j2000,x_rot,y_rot,z_rot,vx_rot,vy_rot,vz_rot\n";

  // ソートされた結果を出力
  for (const auto& result : results) {
    ofs << result.filename << "," << result.jacobi << "," << result.time_j2000 << ","
        << result.state_rot.x << "," << result.state_rot.y << "," << result.state_rot.z << ","
        << result.state_rot.vx << "," << result.state_rot.vy << "," << result.state_rot.vz << "\n";
  }

  ofs.close();

  std::cout << "<>" << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;
  std::cout << "<>    Processing complete!" << std::endl;
  std::cout << "<>    Success: " << success_count << " files" << std::endl;
  std::cout << "<>    Failed:  " << fail_count << " files" << std::endl;
  std::cout << "<>    Output:  " << output_csv << std::endl;
  std::cout << "<>----------------------------------------------------------------" << std::endl;

  return 0;
}
