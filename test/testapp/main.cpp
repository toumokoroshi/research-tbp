#include <fstream>
#include <iostream>
#include <rtbp.hpp>
#include <utils.hpp>
#include <vector>

int main(void) {
  using namespace my_type;
  using namespace utils;

  const State3d<double> center{1.0, 0.0, 0.0};
  const State3d<double> half_width{0.01, 0.01, 0.01};
  const State3d<int> divisions{10, 10, 10};
  std::vector<State3d<double>> mesh =
      createDimensionlessCartesianMesh(center, half_width, divisions);
  std::ofstream ofs{"mesh.txt"};
  for (const auto& x : mesh) {
    ofs << x.x << " " << x.y << " " << x.z << "\n";
  }
  return 0;
}
