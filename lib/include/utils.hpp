#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <rtbp.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace utils {
using namespace my_type;
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
}  // namespace utils
#endif  // UTILS_HPP
