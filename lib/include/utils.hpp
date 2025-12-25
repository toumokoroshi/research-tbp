#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <rtbp.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

namespace utils {
namespace fs = std::filesystem;
using namespace my_type;

// ---------------------------------------------------------------------------
// 共通の列挙型定義
// ---------------------------------------------------------------------------

/**
 * @brief 数値積分法の種類
 */
enum class IntegratorType {
  kSymplectic4th,  ///< 4次シンプレクティック法
  kSymplectic6th,  ///< 6次シンプレクティック法 (デフォルト)
  kRungeKutta4th,  ///< 4次ルンゲクッタ法
  kDormandPrince   ///< ドルマンプリンス法 (DOPRI5/RK45)
};

/**
 * @brief カオス指標の種類
 */
enum class ChaosIndexType {
  NONE,  ///< カオス指標を計算しない
  SALI,  ///< SALI (K=2)
  GALI,  ///< GALI (K可変)
  LLE    ///< 最大リヤプノフ指数
};

#ifndef __KEYWAIT__
#define __KEYWAIT__                                           \
  {                                                           \
    std::cout << "Press any key to continue..." << std::endl; \
    char i;                                                   \
    std::cin >> i;                                            \
  }
#endif

// 入力キー待ち
template <typename T>
inline void WaitForKey(const std::string& message) {
  std::cout << message << std::endl;
  T input;
  std::cin >> input;
}

// エンターキー入力まで待機する関数
inline void WaitForEnter(const std::string& message = "<> Press Enter to continue...") {
  std::cout << message << std::endl;
  while (std::cin.get() != '\n') continue;
}

/**
 * @brief ファイルの各行の開始位置をインデックス化する
 * @param filename ファイルパス
 * @return 各行の開始位置を格納したベクター
 * @throws std::runtime_error ファイルを開けない場合
 */
inline std::vector<std::streampos> indexFileLines(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("ファイルを開けませんでした: " + filename);
  }

  std::vector<std::streampos> linePositions;
  std::string line;

  // 最初の行の開始位置を記録
  linePositions.push_back(file.tellg());

  // 各行の開始位置を記録
  while (std::getline(file, line)) {
    linePositions.push_back(file.tellg());
  }

  return linePositions;
}

/**
 * @brief インデックスを使用して特定の行を高速に読み込む
 * @param filename ファイルパス
 * @param linePositions indexFileLinesで取得したインデックス
 * @param targetLine 読み込む行番号 (1-indexed)
 * @return 指定行の内容
 * @throws std::out_of_range 行番号が範囲外の場合
 */
inline std::string readLineAtIndex(const std::string& filename,
                                   const std::vector<std::streampos>& linePositions,
                                   int targetLine) {
  if (targetLine < 1 || targetLine >= static_cast<int>(linePositions.size())) {
    throw std::out_of_range("指定した行が範囲外です: " + std::to_string(targetLine));
  }

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("ファイルを開けませんでした: " + filename);
  }

  // 指定した行の位置にシーク
  file.seekg(linePositions[targetLine - 1]);

  std::string line;
  std::getline(file, line);
  return line;
}

/**
 * @brief 文字列の前後にある空白文字(スペース, タブなど)を削除する
 * @param s 対象の文字列
 * @return トリム後の文字列
 */
inline std::string trim(const std::string& s) {
  const std::string WHITESPACE = " \t\n\r\f\v";
  size_t first = s.find_first_not_of(WHITESPACE);
  if (std::string::npos == first) {
    return "";
  }
  size_t last = s.find_last_not_of(WHITESPACE);
  return s.substr(first, (last - first + 1));
}

/**
 * @brief キー=値形式の設定ファイルをパースするクラス
 *
 * 設定ファイルの各行を「KEY = VALUE」形式でパースし、
 * 型安全な値取得メソッドを提供する。
 * 同一キーが複数回出現する場合、すべての値を保持する。
 */
class ConfigParser {
 public:
  /**
   * @brief コンストラクタ
   * @param filepath 設定ファイルのパス
   * @throws std::runtime_error ファイルを開けない場合
   */
  explicit ConfigParser(const std::string& filepath) {
    std::ifstream ifs(filepath);
    if (!ifs) {
      throw std::runtime_error("設定ファイルを開けませんでした: " + filepath);
    }
    std::string line;
    while (std::getline(ifs, line)) {
      // 空白を除去
      line = trim(line);
      // コメント行や空行をスキップ
      if (line.empty() || line[0] == '#' || line[0] == ';') {
        continue;
      }
      // '=' で分割
      size_t pos = line.find('=');
      if (pos != std::string::npos) {
        std::string key = trim(line.substr(0, pos));
        std::string value = trim(line.substr(pos + 1));
        config_multimap_[key].push_back(value);
      }
    }
  }

  /**
   * @brief double値を取得（最初の値）
   * @param key キー名
   * @param default_val キーが存在しない場合のデフォルト値
   * @return パースしたdouble値
   */
  double GetDouble(const std::string& key, double default_val = 0.0) const {
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end() && !it->second.empty()) {
      try {
        return std::stod(it->second[0]);
      } catch (...) {
        return default_val;
      }
    }
    return default_val;
  }

  /**
   * @brief int値を取得（最初の値）
   * @param key キー名
   * @param default_val キーが存在しない場合のデフォルト値
   * @return パースしたint値
   */
  int GetInt(const std::string& key, int default_val = 0) const {
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end() && !it->second.empty()) {
      try {
        return std::stoi(it->second[0]);
      } catch (...) {
        return default_val;
      }
    }
    return default_val;
  }

  /**
   * @brief 文字列値を取得（最初の値）
   * @param key キー名
   * @param default_val キーが存在しない場合のデフォルト値
   * @return 文字列値
   */
  std::string GetString(const std::string& key, const std::string& default_val = "") const {
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end() && !it->second.empty()) {
      return it->second[0];
    }
    return default_val;
  }

  /**
   * @brief bool値を取得（最初の値）
   * @param key キー名
   * @param default_val キーが存在しない場合のデフォルト値
   * @return パースしたbool値 (1, true, on -> true)
   */
  bool GetBool(const std::string& key, bool default_val = false) const {
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end() && !it->second.empty()) {
      std::string val = it->second[0];
      // 小文字に変換して比較
      std::transform(val.begin(), val.end(), val.begin(), ::tolower);
      if (val == "1" || val == "true" || val == "on" || val == "yes") {
        return true;
      } else if (val == "0" || val == "false" || val == "off" || val == "no") {
        return false;
      }
    }
    return default_val;
  }

  /**
   * @brief スペース区切りの3つの値をState3dとして取得（最初の値）
   * @param key キー名
   * @param default_val キーが存在しない場合のデフォルト値
   * @return State3d値
   */
  template <typename T>
  State3d<T> GetState3d(const std::string& key, State3d<T> default_val = {}) const {
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end() && !it->second.empty()) {
      std::stringstream ss(it->second[0]);
      T x, y, z;
      if (ss >> x >> y >> z) {
        return State3d<T>{x, y, z};
      }
    }
    return default_val;
  }

  // ========== 複数値取得メソッド ==========

  /**
   * @brief 同一キーのすべての文字列値を取得
   * @param key キー名
   * @return 文字列のベクター（キーが存在しない場合は空）
   */
  std::vector<std::string> GetMultiString(const std::string& key) const {
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end()) {
      return it->second;
    }
    return {};
  }

  /**
   * @brief 同一キーのすべてのdouble値を取得
   * @param key キー名
   * @return doubleのベクター（パース失敗した値はスキップ）
   */
  std::vector<double> GetMultiDouble(const std::string& key) const {
    std::vector<double> result;
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end()) {
      for (const auto& val : it->second) {
        try {
          result.push_back(std::stod(val));
        } catch (...) {
          // パース失敗はスキップ
        }
      }
    }
    return result;
  }

  /**
   * @brief 同一キーのすべてのカンマ区切り6値をState<double>として取得
   * @param key キー名
   * @return State<double>のベクター（パース失敗した値はスキップ）
   * @note COORD = x, y, z, vx, vy, vz 形式のパースに使用
   */
  std::vector<State<double>> GetMultiState6d(const std::string& key) const {
    std::vector<State<double>> result;
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end()) {
      for (const auto& val : it->second) {
        // カンマ区切りでパース
        std::vector<double> values;
        std::stringstream ss(val);
        std::string token;
        while (std::getline(ss, token, ',')) {
          std::string trimmed = trim(token);
          if (!trimmed.empty()) {
            try {
              values.push_back(std::stod(trimmed));
            } catch (...) {
              break;
            }
          }
        }
        if (values.size() >= 6) {
          result.push_back(
              State<double>{values[0], values[1], values[2], values[3], values[4], values[5]});
        }
      }
    }
    return result;
  }

  /**
   * @brief キーが存在するか確認
   * @param key キー名
   * @return 存在する場合true
   */
  bool HasKey(const std::string& key) const {
    return config_multimap_.find(key) != config_multimap_.end();
  }

  /**
   * @brief 指定キーの値の数を取得
   * @param key キー名
   * @return 値の数（キーが存在しない場合は0）
   */
  size_t GetValueCount(const std::string& key) const {
    auto it = config_multimap_.find(key);
    if (it != config_multimap_.end()) {
      return it->second.size();
    }
    return 0;
  }

  /**
   * @brief すべてのキーを取得（重複なし）
   * @return キーのベクター
   */
  std::vector<std::string> GetKeys() const {
    std::vector<std::string> keys;
    for (const auto& pair : config_multimap_) {
      keys.push_back(pair.first);
    }
    return keys;
  }

 private:
  std::map<std::string, std::vector<std::string>> config_multimap_;
};

/**
 * @brief TOML形式の設定ファイルをパースするクラス
 *
 * toml++ライブラリを使用してTOML形式の設定ファイルを読み込み、
 * 階層的なキーアクセスと配列取得をサポートする。
 *
 * 使用例:
 *   TomlConfigParser config("config.toml");
 *   double dt = config.GetDouble("simulation.calc_timestep", 0.001);
 *   auto coords = config.GetCoordsArray("coords");
 */
class TomlConfigParser {
 public:
  /**
   * @brief コンストラクタ
   * @param filepath TOML設定ファイルのパス
   * @throws std::runtime_error ファイルを開けない/パースエラーの場合
   */
  explicit TomlConfigParser(const std::string& filepath) {
    try {
      config_ = toml::parse_file(filepath);
    } catch (const toml::parse_error& err) {
      throw std::runtime_error("TOMLパースエラー: " + std::string(err.description()));
    }
  }

  /**
   * @brief double値を取得（ドット区切りのパスをサポート）
   * @param path キーパス (例: "simulation.calc_timestep")
   * @param default_val キーが存在しない場合のデフォルト値
   * @return パースしたdouble値
   */
  double GetDouble(const std::string& path, double default_val = 0.0) const {
    auto node = navigateToNode(path);
    if (node && node->is_number()) {
      return node->value_or(default_val);
    }
    return default_val;
  }

  /**
   * @brief int値を取得
   * @param path キーパス
   * @param default_val キーが存在しない場合のデフォルト値
   * @return パースしたint値
   */
  int GetInt(const std::string& path, int default_val = 0) const {
    auto node = navigateToNode(path);
    if (node && node->is_integer()) {
      return static_cast<int>(node->value_or(static_cast<int64_t>(default_val)));
    }
    return default_val;
  }

  /**
   * @brief 文字列値を取得
   * @param path キーパス
   * @param default_val キーが存在しない場合のデフォルト値
   * @return 文字列値
   */
  std::string GetString(const std::string& path, const std::string& default_val = "") const {
    auto node = navigateToNode(path);
    if (node && node->is_string()) {
      return std::string(node->value_or(default_val));
    }
    return default_val;
  }

  /**
   * @brief bool値を取得
   * @param path キーパス
   * @param default_val キーが存在しない場合のデフォルト値
   * @return bool値
   */
  bool GetBool(const std::string& path, bool default_val = false) const {
    auto node = navigateToNode(path);
    if (node && node->is_boolean()) {
      return node->value_or(default_val);
    }
    return default_val;
  }

  /**
   * @brief double配列を取得
   * @param path キーパス
   * @return doubleのベクター
   */
  std::vector<double> GetDoubleArray(const std::string& path) const {
    std::vector<double> result;
    auto node = navigateToNode(path);
    if (node && node->is_array()) {
      for (const auto& elem : *node->as_array()) {
        if (elem.is_number()) {
          result.push_back(elem.value_or(0.0));
        }
      }
    }
    return result;
  }

  /**
   * @brief 座標配列を取得 (TOML配列のテーブルから)
   * @param path キーパス (例: "coords")
   * @return State<double>のベクター
   *
   * 期待するTOML形式:
   *   [[coords]]
   *   position = [1.0, 0.0, 0.0]
   *   velocity = [0.1, 0.0, 0.0]
   */
  std::vector<State<double>> GetCoordsArray(const std::string& path) const {
    std::vector<State<double>> result;
    auto node = navigateToNode(path);
    if (node && node->is_array()) {
      for (const auto& elem : *node->as_array()) {
        if (elem.is_table()) {
          const auto& tbl = *elem.as_table();
          auto pos = tbl["position"].as_array();
          auto vel = tbl["velocity"].as_array();
          if (pos && pos->size() >= 3 && vel && vel->size() >= 3) {
            State<double> state;
            state.x = (*pos)[0].value_or(0.0);
            state.y = (*pos)[1].value_or(0.0);
            state.z = (*pos)[2].value_or(0.0);
            state.vx = (*vel)[0].value_or(0.0);
            state.vy = (*vel)[1].value_or(0.0);
            state.vz = (*vel)[2].value_or(0.0);
            result.push_back(state);
          }
        }
      }
    }
    return result;
  }

  /**
   * @brief State3d配列を取得
   * @param path キーパス
   * @return State3d<double>のベクター
   */
  template <typename T>
  State3d<T> GetState3d(const std::string& path, State3d<T> default_val = {}) const {
    auto arr = GetDoubleArray(path);
    if (arr.size() >= 3) {
      return State3d<T>{static_cast<T>(arr[0]), static_cast<T>(arr[1]), static_cast<T>(arr[2])};
    }
    return default_val;
  }

  /**
   * @brief キーが存在するか確認
   * @param path キーパス
   * @return 存在する場合true
   */
  bool HasKey(const std::string& path) const { return navigateToNode(path) != nullptr; }

  /**
   * @brief 生のtoml::tableへのアクセス（上級者向け）
   * @return toml::tableへの参照
   */
  const toml::table& GetRawTable() const { return config_; }

 private:
  /**
   * @brief ドット区切りのパスをナビゲートしてノードを取得
   * @param path キーパス (例: "simulation.calc_timestep")
   * @return ノードへのポインタ（存在しない場合はnullptr）
   */
  const toml::node* navigateToNode(const std::string& path) const {
    std::vector<std::string> keys;
    std::stringstream ss(path);
    std::string token;
    while (std::getline(ss, token, '.')) {
      keys.push_back(token);
    }

    const toml::node* current = &config_;
    for (const auto& key : keys) {
      if (current->is_table()) {
        current = current->as_table()->get(key);
        if (!current) return nullptr;
      } else {
        return nullptr;
      }
    }
    return current;
  }

  toml::table config_;
};

/**
 * @brief シミュレーション出力ファイルを書き込むためのクラス
 *
 * ヘッダーコメント、カラム名、データ行を一貫した形式で出力する。
 */
class SimulationOutputWriter {
 public:
  /**
   * @brief コンストラクタ
   * @param filepath 出力ファイルのパス
   * @param precision 数値の精度 (デフォルト: 15)
   */
  explicit SimulationOutputWriter(const std::string& filepath, int precision = 15)
      : filepath_(filepath), precision_(precision), header_written_(false) {
    ofs_.open(filepath);
    if (!ofs_) {
      throw std::runtime_error("出力ファイルを開けませんでした: " + filepath);
    }
    ofs_ << std::setprecision(precision_) << std::fixed;
  }

  /**
   * @brief ヘッダーコメントを追加 (文字列値)
   * @param key キー名
   * @param value 値
   */
  void AddHeaderComment(const std::string& key, const std::string& value) {
    ofs_ << "# " << key << "=" << value << "\n";
  }

  /**
   * @brief ヘッダーコメントを追加 (数値)
   * @param key キー名
   * @param value 数値
   */
  void AddHeaderComment(const std::string& key, double value) {
    ofs_ << "# " << key << "=" << value << "\n";
  }

  /**
   * @brief ヘッダーコメントを追加 (整数)
   * @param key キー名
   * @param value 整数
   */
  void AddHeaderComment(const std::string& key, int value) {
    ofs_ << "# " << key << "=" << value << "\n";
  }

  /**
   * @brief カラム名を設定 (CSV形式)
   * @param column_names カラム名のベクター
   */
  void SetColumns(const std::vector<std::string>& column_names) {
    for (size_t i = 0; i < column_names.size(); ++i) {
      ofs_ << column_names[i];
      if (i < column_names.size() - 1) {
        ofs_ << ",";
      }
    }
    ofs_ << "\n";
    header_written_ = true;
  }

  /**
   * @brief データ行を書き込む (可変引数テンプレート)
   * @param values 書き込む値
   */
  template <typename... Args>
  void WriteRow(Args&&... values) {
    WriteRowImpl(std::forward<Args>(values)...);
    ofs_ << "\n";
  }

  /**
   * @brief スペース区切りでデータ行を書き込む
   * @param values 書き込む値
   */
  template <typename... Args>
  void WriteRowSpace(Args&&... values) {
    WriteRowSpaceImpl(std::forward<Args>(values)...);
    ofs_ << "\n";
  }

  /**
   * @brief ファイルを閉じる
   */
  void Close() {
    if (ofs_.is_open()) {
      ofs_.close();
    }
  }

  /**
   * @brief 内部ストリームへの参照を取得（複雑な出力用）
   * @return 出力ストリームへの参照
   */
  std::ofstream& GetStream() { return ofs_; }

  /**
   * @brief 空行を書き込む（gnuplotインデックス区切り用など）
   */
  void WriteBlankLine() { ofs_ << "\n"; }

  /**
   * @brief コメント行を書き込む
   * @param comment コメント内容
   */
  void WriteComment(const std::string& comment) { ofs_ << "# " << comment << "\n"; }

  /**
   * @brief ファイルパスを取得
   * @return 出力ファイルのパス
   */
  std::string GetFilePath() const { return filepath_; }

  /**
   * @brief デストラクタ
   */
  ~SimulationOutputWriter() { Close(); }

 private:
  // カンマ区切り用の再帰実装
  template <typename T>
  void WriteRowImpl(T&& value) {
    ofs_ << std::forward<T>(value);
  }

  template <typename T, typename... Rest>
  void WriteRowImpl(T&& first, Rest&&... rest) {
    ofs_ << std::forward<T>(first) << ",";
    WriteRowImpl(std::forward<Rest>(rest)...);
  }

  // スペース区切り用の再帰実装
  template <typename T>
  void WriteRowSpaceImpl(T&& value) {
    ofs_ << std::forward<T>(value);
  }

  template <typename T, typename... Rest>
  void WriteRowSpaceImpl(T&& first, Rest&&... rest) {
    ofs_ << std::forward<T>(first) << " ";
    WriteRowSpaceImpl(std::forward<Rest>(rest)...);
  }

  std::ofstream ofs_;
  std::string filepath_;
  int precision_;
  bool header_written_;
};

template <typename ScalarType>
inline std::vector<State3d<ScalarType>> createSphereMesh(const ScalarType ROI_radius,
                                                         const int divisions,
                                                         const State3d<ScalarType>& center) {
  std::vector<State3d<ScalarType>> meshPoints;
  if (divisions <= 0) return meshPoints;

  ScalarType step = (2.0 * ROI_radius) / (divisions - 1);
  ScalarType radiusSquared = ROI_radius * ROI_radius;
  const ScalarType cx = center.x;
  const ScalarType cy = center.y;
  const ScalarType cz = center.z;
  for (int i = 0; i < divisions; ++i) {
    for (int j = 0; j < divisions; ++j) {
      for (int k = 0; k < divisions; ++k) {
        ScalarType x = cx - ROI_radius + i * step;
        ScalarType y = cy - ROI_radius + j * step;
        ScalarType z = cz - ROI_radius + k * step;

        ScalarType distanceSquared =
            std::pow(x - cx, 2) + std::pow(y - cy, 2) + std::pow(z - cz, 2);
        if (distanceSquared <= radiusSquared) {
          // 初回だけプリント
          if (meshPoints.empty()) {
            std::cout << " <> mesh point : " << x << ", " << y << ", " << z << std::endl;
          }
          meshPoints.push_back({x, y, z});
        }
      }
    }
  }

  std::sort(meshPoints.begin(), meshPoints.end(),
            [](const State3d<ScalarType>& a, const State3d<ScalarType>& b) {
              if (a.z != b.z) return a.z < b.z;
              if (a.y != b.y) return a.y < b.y;
              return a.x < b.x;
            });

  return meshPoints;
}
template <typename ScalarType>
inline std::vector<State3d<ScalarType>> CreateCircleMesh(const ScalarType ROI_radius,
                                                         const int divisions,
                                                         const State3d<ScalarType>& center) {
  std::vector<State3d<ScalarType>> meshPoints;
  if (divisions <= 0) return meshPoints;

  ScalarType step = (2.0 * ROI_radius) / (divisions - 1);
  ScalarType radiusSquared = ROI_radius * ROI_radius;

  const ScalarType cx = center.x;
  const ScalarType cy = center.y;
  const ScalarType cz = center.z;
  for (int i = 0; i < divisions; ++i) {
    for (int j = 0; j < divisions; ++j) {
      ScalarType x = cx - ROI_radius + i * step;
      ScalarType y = cy - ROI_radius + j * step;
      ScalarType z = cz;

      ScalarType distanceSquared = std::pow(x - cx, 2) + std::pow(y - cy, 2) + std::pow(z - cz, 2);
      if (distanceSquared <= radiusSquared) {
        // 初回だけプリント
        if (meshPoints.empty()) {
          std::cout << " <> mesh point : " << x << ", " << y << ", " << z << std::endl;
        }
        meshPoints.push_back({x, y, z});
      }
    }
  }

  std::sort(meshPoints.begin(), meshPoints.end(),
            [](const State3d<ScalarType>& a, const State3d<ScalarType>& b) {
              if (a.y != b.y) return a.y < b.y;
              return a.x < b.x;
            });

  return meshPoints;
}

template <typename ScalarType>
inline std::vector<State3d<ScalarType>> createCubeMesh(const State3d<ScalarType>& center,
                                                       const ScalarType spacing,
                                                       const ScalarType halfExtentInSteps = 1) {
  std::vector<State3d<ScalarType>> meshPoints;
  if (spacing <= ScalarType(0) || halfExtentInSteps < 0) return meshPoints;

  const int gridWidth = halfExtentInSteps * 2 + 1;
  meshPoints.reserve(static_cast<std::size_t>(gridWidth) * gridWidth * gridWidth);

  for (int dx = -halfExtentInSteps; dx <= halfExtentInSteps; ++dx) {
    for (int dy = -halfExtentInSteps; dy <= halfExtentInSteps; ++dy) {
      for (int dz = -halfExtentInSteps; dz <= halfExtentInSteps; ++dz) {
        meshPoints.push_back({center.x + static_cast<ScalarType>(dx) * spacing,
                              center.y + static_cast<ScalarType>(dy) * spacing,
                              center.z + static_cast<ScalarType>(dz) * spacing});
      }
    }
  }

  std::sort(meshPoints.begin(), meshPoints.end(),
            [](const State3d<ScalarType>& a, const State3d<ScalarType>& b) {
              if (a.z != b.z) return a.z < b.z;
              if (a.y != b.y) return a.y < b.y;
              return a.x < b.x;
            });

  return meshPoints;
}

/**
 * @brief Create a dimensionless Cartesian mesh around @p center.
 * @param center    Mesh center in dimensionless coordinates.
 * @param halfWidth Positive half-width per axis (extent from center to each boundary).
 * @param divisions Number of samples per axis (>= 1).
 */
template <typename ScalarType>
inline std::vector<State3d<ScalarType>> createDimensionlessCartesianMesh(
    const State3d<ScalarType>& center, const State3d<ScalarType>& halfWidth,
    const State3d<int>& divisions) {
  if (halfWidth.x < ScalarType(0) || halfWidth.y < ScalarType(0) || halfWidth.z < ScalarType(0)) {
    throw std::invalid_argument(
        "createDimensionlessCartesianMesh: halfWidth must be non-negative.");
  }
  if (divisions.x <= 0 || divisions.y <= 0 || divisions.z <= 0) {
    throw std::invalid_argument("createDimensionlessCartesianMesh: divisions must be positive.");
  }

  const State3d<ScalarType> minCorner{center.x - halfWidth.x, center.y - halfWidth.y,
                                      center.z - halfWidth.z};
  const auto computeStep = [](ScalarType width, int div) -> ScalarType {
    return (div > 1) ? (width * ScalarType(2) / static_cast<ScalarType>(div - 1)) : ScalarType(0);
  };
  const ScalarType stepX = computeStep(halfWidth.x, divisions.x);
  const ScalarType stepY = computeStep(halfWidth.y, divisions.y);
  const ScalarType stepZ = computeStep(halfWidth.z, divisions.z);

  const std::size_t totalPoints = static_cast<std::size_t>(divisions.x) *
                                  static_cast<std::size_t>(divisions.y) *
                                  static_cast<std::size_t>(divisions.z);
  std::vector<State3d<ScalarType>> meshPoints;
  meshPoints.reserve(totalPoints);

  for (int ix = 0; ix < divisions.x; ++ix) {
    for (int iy = 0; iy < divisions.y; ++iy) {
      for (int iz = 0; iz < divisions.z; ++iz) {
        meshPoints.push_back({minCorner.x + static_cast<ScalarType>(ix) * stepX,
                              minCorner.y + static_cast<ScalarType>(iy) * stepY,
                              minCorner.z + static_cast<ScalarType>(iz) * stepZ});
      }
    }
  }

  return meshPoints;
}

/**
 * @brief Convert a mesh from dimensionless units into simulation units.
 * @param mesh        Dimensionless mesh (e.g., output of createDimensionlessCartesianMesh()).
 * @param scale       Per-axis scaling factors to reach the target unit system.
 * @param offset      Optional translation applied after scaling (defaults to origin).
 * @return            Mesh represented in simulation units.
 */
template <typename ScalarType>
inline std::vector<State3d<ScalarType>> convertMeshToSimulationUnits(
    const std::vector<State3d<ScalarType>>& mesh, const State3d<ScalarType>& scale,
    const State3d<ScalarType>& offset = State3d<ScalarType>{0, 0, 0}) {
  std::vector<State3d<ScalarType>> converted;
  converted.reserve(mesh.size());

  for (const auto& point : mesh) {
    converted.push_back(
        {point.x * scale.x + offset.x, point.y * scale.y + offset.y, point.z * scale.z + offset.z});
  }

  return converted;
}

// template <typename T, typename ScalarType>
// inline std::vector<T> create_cube_mesh(const ScalarType ROI_length, const int divisions,
// }
// template <typename ScalarType>
// inline std::vector<Coord3D<ScalarType>> create_square_mesh(
//     const ScalarType ROI_length, const int divisions,
//     const Coord3D<ScalarType>& center = Coord3D<ScalarType>(0, 0, 0), const double z = 0) {
//   std::vector<Coord3D<ScalarType>> meshPoints;
//   if (divisions <= 0) return meshPoints;

//   double step = (2 * ROI_length) / (divisions - 1);

//   for (int i = 0; i < divisions; ++i) {
//     for (int j = 0; j < divisions; ++j) {
//       for (int k = 0; k < divisions; ++k) {
//         ScalarType x = center[0] - ROI_length + i * step;
//         ScalarType y = center[1] - ROI_length + j * step;
//         ScalarType z = center[2] - ROI_length + j * step;
//         meshPoints.push_back({x, y, z});
//       }
//     }
//   }

//   std::sort(meshPoints.begin(), meshPoints.end(),
//             [](const Coord3D<ScalarType>& a, const Coord3D<ScalarType>& b) {
//               if (a[2] != b[2]) return a[2] < b[2];
//               if (a[1] != b[1]) return a[1] < b[1];
//               return a[0] < b[0];
//             });

//   return meshPoints;
// }

inline void displayProgressBar(double progress, int barWidth = 40) {
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      std::cout << "=";
    else if (i == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

// スレッドセーフな進捗表示関数
inline void displayProgressBarThreadSafe(double progress, int barWidth = 40) {
#pragma omp critical(progress_display)
  {
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos)
        std::cout << "=";
      else if (i == pos)
        std::cout << ">";
      else
        std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }
}

inline std::string getcurrent_date() {
  // 現在の時刻を取得
  std::time_t now = std::time(nullptr);
  std::tm* ltm = std::localtime(&now);

  // 日付時刻をフォーマット
  std::ostringstream oss;
  oss << std::put_time(ltm, "%y_%m%d_%H%M");
  return oss.str();
}

/**
 * @brief 設定ファイルを読み込み、キーと値のマップを返す
 * @param filename 読み込むファイル名
 * @return キー(string)と値(long double)の std::map
 * @throws std::runtime_error ファイルが開けない場合
 */
template <typename ScalarType>
inline param::AstroConstants<ScalarType> loadConstants(const std::string& filename) {
  std::cout << "<>-------------------------------" << std::endl;

  std::map<std::string, double> constants;
  std::ifstream file(filename);
  param::AstroConstants<ScalarType> astroConstants;

  if (!file.is_open()) {
    // ファイルが開けなかった場合、例外を投げる
    throw std::runtime_error("エラー: ファイルを開けません: " + filename);
  }
  std::cout << "----" << std::endl;

  std::string line;
  int line_number = 0;

  while (std::getline(file, line)) {
    line_number++;

    size_t comment_pos = line.find("//");
    if (comment_pos != std::string::npos) {
      line = line.substr(0, comment_pos);
    }

    size_t eq_pos = line.find('=');
    if (eq_pos == std::string::npos) {
      continue;
    }

    std::string key = trim(line.substr(0, eq_pos));
    std::string val_str = trim(line.substr(eq_pos + 1));

    if (key.empty() || val_str.empty()) {
      continue;
    }

    try {
      // マップに直接キーと値を格納
      constants[key] = std::stold(val_str);
    } catch (const std::exception& e) {
      std::cerr << "警告 (L" << line_number << "): 値のパースに失敗しました: '" << val_str << "'"
                << std::endl;
    }
  }

  file.close();

  astroConstants.au = constants.at("au");
  astroConstants.gm_sun = constants.at("gm_sun");
  astroConstants.gm_earth = constants.at("gm_earth");
  astroConstants.G = constants.at("G");
  return astroConstants;  // 読み込んだ構造体を返す
}

/**
 * @brief プレフィックスに一致する設定ファイルをディレクトリから列挙。
 */
inline std::vector<std::string> DiscoverConfigFiles(const std::string& directory,
                                                    const std::string& prefix) {
  std::vector<std::string> files;
  if (!fs::exists(directory)) {
    std::cerr << "Config directory does not exist: " << directory << "\n";
    return files;
  }
  const std::regex pattern("^" + prefix + "(?:_\\d+)?\\.txt$");
  try {
    for (const auto& entry : fs::directory_iterator(directory)) {
      if (!entry.is_regular_file()) continue;
      const auto name = entry.path().filename().string();
      if (std::regex_match(name, pattern)) {
        files.push_back(fs::absolute(entry.path()).string());
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "Failed to read config directory: " << e.what() << "\n";
    return files;
  }

  auto sorter = [](const std::string& a, const std::string& b) {
    const std::string stem_a = fs::path(a).stem().string();
    const std::string stem_b = fs::path(b).stem().string();
    auto number_from_stem = [](const std::string& stem) -> int {
      const auto pos = stem.find_last_of('_');
      if (pos == std::string::npos) return 0;
      return std::stoi(stem.substr(pos + 1));
    };
    return number_from_stem(stem_a) < number_from_stem(stem_b);
  };
  std::sort(files.begin(), files.end(), sorter);
  return files;
}

// ---------------------------------------------------------------------------
// TOML設定ファイル検出用の共通構造体・関数
// ---------------------------------------------------------------------------

/**
 * @brief 設定ファイル検出オプション
 */
struct ConfigDiscoveryOptions {
  bool exclude_sample = true;  ///< _sample付きを除外するか
  bool continuous_mode =
      false;  ///< 連番ファイルのみ対象（true: prefix_1.toml等、false: prefix.toml）
};

/**
 * @brief TOML設定ファイルをディレクトリから列挙
 * @param directory 検索ディレクトリ
 * @param prefix ファイル名プレフィックス（例: "3DSALIconfig"）
 * @param options 検出オプション
 * @return 検出したファイルパスのベクター（数値順にソート済み）
 */
inline std::vector<std::string> DiscoverConfigFilesToml(
    const std::string& directory, const std::string& prefix,
    const ConfigDiscoveryOptions& options = {}) {
  std::vector<std::string> files;
  if (!fs::exists(directory)) {
    std::cerr << "<> Config directory does not exist: " << directory << std::endl;
    return files;
  }

  // パターン: continuous_mode=true -> prefix_数字.toml、false -> prefix.toml
  std::string pattern_str;
  if (options.continuous_mode) {
    pattern_str = "^" + prefix + "_\\d+\\.toml$";
  } else {
    pattern_str = "^" + prefix + "\\.toml$";
  }
  const std::regex pattern(pattern_str);

  // _sample除外用パターン
  const std::regex sample_pattern(".*_sample.*", std::regex_constants::icase);

  try {
    for (const auto& entry : fs::directory_iterator(directory)) {
      if (!entry.is_regular_file()) continue;
      const auto filename = entry.path().filename().string();

      // _sample除外チェック
      if (options.exclude_sample && std::regex_match(filename, sample_pattern)) {
        continue;
      }

      if (std::regex_match(filename, pattern)) {
        files.push_back(fs::absolute(entry.path()).string());
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "<> Failed to read config directory: " << e.what() << std::endl;
    return files;
  }

  // 数値順にソート
  auto sorter = [](const std::string& a, const std::string& b) {
    const std::string stem_a = fs::path(a).stem().string();
    const std::string stem_b = fs::path(b).stem().string();
    auto number_from_stem = [](const std::string& stem) -> int {
      const auto pos = stem.find_last_of('_');
      if (pos == std::string::npos) return 0;
      try {
        return std::stoi(stem.substr(pos + 1));
      } catch (...) {
        return 0;
      }
    };
    return number_from_stem(stem_a) < number_from_stem(stem_b);
  };
  std::sort(files.begin(), files.end(), sorter);
  return files;
}

// ---------------------------------------------------------------------------
// コマンドライン引数パース用の共通構造体・関数
// ---------------------------------------------------------------------------

/**
 * @brief コマンドライン引数の解析結果を格納する構造体
 */
struct CommonArgs {
  bool is_continuous = false;   ///< 連続シミュレーションモード
  bool skip_wait = false;       ///< WaitForEnterをスキップ
  std::string output_tag = "";  ///< 出力フォルダに付与するタグ
};

/**
 * @brief 共通のコマンドライン引数をパースする
 * @param argc 引数の数
 * @param argv 引数の配列
 * @param app_name アプリケーション名（ヘルプ表示用、空の場合argv[0]を使用）
 * @return パースされた引数
 *
 * サポートするオプション:
 *   -c, --continuous  連続シミュレーションモード
 *   -n, --no-wait     ユーザー確認の待機をスキップ
 *   -t, --tag <TAG>   出力フォルダに付与するタグ
 *   -h, --help        ヘルプを表示して終了
 */
inline CommonArgs ParseCommonArgs(int argc, char* argv[], const std::string& app_name = "") {
  CommonArgs args;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--continuous" || arg == "-c") {
      args.is_continuous = true;
    } else if (arg == "--no-wait" || arg == "-n") {
      args.skip_wait = true;
    } else if ((arg == "--tag" || arg == "-t") && i + 1 < argc) {
      args.output_tag = argv[++i];
    } else if (arg == "--help" || arg == "-h") {
      std::string name = app_name.empty() ? argv[0] : app_name;
      std::cout
          << "Usage: " << name << " [options]\n"
          << "Options:\n"
          << "  -c, --continuous  連続シミュレーションモード（複数configファイルを順次処理）\n"
          << "  -n, --no-wait     ユーザー確認のための待機をスキップ\n"
          << "  -t, --tag <TAG>   出力フォルダに付与するタグ\n"
          << "  -h, --help        このヘルプを表示\n"
          << std::endl;
      std::exit(0);
    }
  }
  return args;
}

}  // namespace utils
#endif  // UTILS_HPP
