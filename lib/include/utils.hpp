#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifndef __KEYWAIT__
#define __KEYWAIT__                                                            \
  {                                                                            \
    std::cout << "Press any key to continue..." << std::endl;                  \
    char i;                                                                    \
    std::cin >> i;                                                             \
  }
#endif

// 入力キー待ち
template <typename T> inline void WaitForKey(const std::string &message) {
  std::cout << message << std::endl;
  T input;
  std::cin >> input;
}

// エンターキー入力まで待機する関数
inline void
WaitForEnter(const std::string &message = "Press Enter to continue...") {
  std::cout << message << std::endl;
  while (std::cin.get() != '\n')
    continue;
}

struct Point3D {
  double x, y, z;
  Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

template <typename T>
inline std::vector<T> createSphereMesh(const double ROI_radius,
                                       const int divisions, const T &center) {
  std::vector<T> meshPoints;
  if (divisions <= 0)
    return meshPoints;

  double step = (2.0 * ROI_radius) / (divisions - 1);
  double radiusSquared = ROI_radius * ROI_radius;

  for (int i = 0; i < divisions; ++i) {
    for (int j = 0; j < divisions; ++j) {
      for (int k = 0; k < divisions; ++k) {
        double cx, cy, cz;
        if constexpr (std::is_same<T, Point3D>::value) {
          cx = center.x;
          cy = center.y;
          cz = center.z;
        } else if constexpr (std::is_same<T, std::array<double, 3>>::value) {
          cx = center[0];
          cy = center[1];
          cz = center[2];
        }
        double x = cx - ROI_radius + i * step;
        double y = cy - ROI_radius + j * step;
        double z = cz - ROI_radius + k * step;

        double distanceSquared =
            std::pow(x - cx, 2) + std::pow(y - cy, 2) + std::pow(z - cz, 2);
        if (distanceSquared <= radiusSquared) {
          // 初回だけプリント
          if (meshPoints.empty()) {
            std::cout << " <> mesh point : " << x << ", " << y << ", " << z
                      << std::endl;
          }
          meshPoints.emplace_back(x, y, z);
        }
      }
    }
  }

  std::sort(meshPoints.begin(), meshPoints.end(), [](const T &a, const T &b) {
    if constexpr (std::is_same<T, Point3D>::value) {
      if (a.z != b.z)
        return a.z < b.z;
      if (a.y != b.y)
        return a.y < b.y;
      return a.x < b.x;
    } else if constexpr (std::is_same<T, std::array<double, 3>>::value) {
      if (a[2] != b[2])
        return a[2] < b[2];
      if (a[1] != b[1])
        return a[1] < b[1];
      return a[0] < b[0];
    }
  });

  return meshPoints;
}

template <typename T>
inline std::vector<T> create_cube_mesh(const double ROI_length,
                                       const int divisions, const T &center) {
  std::vector<T> meshPoints;
  if (divisions <= 0)
    return meshPoints;

  double step = (ROI_length) / (divisions - 1);

  for (int i = 0; i < divisions; ++i) {
    for (int j = 0; j < divisions; ++j) {
      for (int k = 0; k < divisions; ++k) {
        if constexpr (std::is_same<T, Point3D>::value) {
          double x = center.x - ROI_length + i * step;
          double y = center.y - ROI_length + j * step;
          double z = center.z - ROI_length + k * step;
          meshPoints.emplace_back(x, y, z);
        } else if constexpr (std::is_same<T, std::array<double, 3>>::value) {
          double x = center[0] - ROI_length + i * step;
          double y = center[1] - ROI_length + j * step;
          double z = center[2] - ROI_length + k * step;
          meshPoints.emplace_back(x, y, z);
        }
      }
    }
  }

  std::sort(meshPoints.begin(), meshPoints.end(), [](const T &a, const T &b) {
    if constexpr (std::is_same<T, Point3D>::value) {
      if (a.z != b.z)
        return a.z < b.z;
      if (a.y != b.y)
        return a.y < b.y;
      return a.x < b.x;
    } else if constexpr (std::is_same<T, std::array<double, 3>>::value) {
      if (a[2] != b[2])
        return a[2] < b[2];
      if (a[1] != b[1])
        return a[1] < b[1];
      return a[0] < b[0];
    }
  });

  return meshPoints;
}

inline std::vector<Point3D>
create_square_mesh(const double ROI_length, const int divisions,
                   const Point3D &center = Point3D(0, 0, 0),
                   const double z = 0) {
  std::vector<Point3D> meshPoints;
  if (divisions <= 0)
    return meshPoints;

  double step = (2 * ROI_length) / (divisions - 1);

  for (int i = 0; i < divisions; ++i) {
    for (int j = 0; j < divisions; ++j) {
      for (int k = 0; k < divisions; ++k) {
        double x = center.x - ROI_length + i * step;
        double y = center.y - ROI_length + j * step;
        meshPoints.emplace_back(x, y, z);
      }
    }
  }

  std::sort(meshPoints.begin(), meshPoints.end(),
            [](const Point3D &a, const Point3D &b) {
              if (a.z != b.z)
                return a.z < b.z;
              if (a.y != b.y)
                return a.y < b.y;
              return a.x < b.x;
            });

  return meshPoints;
}

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

inline std::string getcurrent_date() {
  // 現在の時刻を取得
  std::time_t now = std::time(nullptr);
  std::tm *ltm = std::localtime(&now);

  // 日付時刻をフォーマット
  std::ostringstream oss;
  oss << std::put_time(ltm, "%y_%m%d_%H%M");
  return oss.str();
}

/**
 * @brief 構造体の概要説明
 *
 * 詳細説明
 */
struct AstroConstants {
  double au;       ///< 天文単位 (m)
  double gm_sun;   ///< 太陽の重力定数 (m^3/s^2)
  double gm_earth; ///< 地球の重力定数 (m^3/s^2)
  double mu;       ///< 質量比
  double G;        ///< 万有引力定数 (kg^-1 m^3 s^-2)
};
/**
 * @brief 文字列の前後にある空白文字(スペース, タブなど)を削除する
 * @param s 対象の文字列
 * @return トリム後の文字列
 */
std::string trim(const std::string &s) {
  const std::string WHITESPACE = " \t\n\r\f\v";
  size_t first = s.find_first_not_of(WHITESPACE);
  if (std::string::npos == first) {
    return "";
  }
  size_t last = s.find_last_not_of(WHITESPACE);
  return s.substr(first, (last - first + 1));
}
/**
 * @brief 設定ファイルを読み込み、キーと値のマップを返す
 * @param filename 読み込むファイル名
 * @return キー(string)と値(long double)の std::map
 * @throws std::runtime_error ファイルが開けない場合
 */
inline AstroConstants loadConstants(const std::string &filename) {

  std::map<std::string, double> constants;
  std::ifstream file(filename);
  AstroConstants astroConstants;

  if (!file.is_open()) {
    // ファイルが開けなかった場合、例外を投げる
    throw std::runtime_error("エラー: ファイルを開けません: " + filename);
  }

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
    } catch (const std::exception &e) {
      std::cerr << "警告 (L" << line_number << "): 値のパースに失敗しました: '"
                << val_str << "'" << std::endl;
    }
  }

  file.close();

  astroConstants.au = constants.at("au") astroConstants.gm_sun =
      constants.at("gm_sun");
  astroConstants.gm_earth = constants.at("gm_earth");
  astroConstants.G = constants.at("G");
  astroConstants.mu = astroConstants.gm_earth /
                      (astroConstants.gm_sun +
                       astroConstants.gm_earth); // mu parameter of Earth-Sun
  return astroConstants;                         // 読み込んだ構造体を返す
}
#endif // UTILS_HPP