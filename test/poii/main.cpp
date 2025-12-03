#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <crtbp.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector3d.hpp>
#include <vector>

int main(void) {
  using namespace crtbp;
  using namespace utils;
  using namespace my_type;
  const std::string kConfigDir = std::string(CONFIG_DIR);

  AstroConstants<double> astro_params =
      loadConstants<double>(kConfigDir + "/astro_param/astro_param.txt");
  const double kMU = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);
  //   8.968485078526299e-01 4.333538985868564e-01 -5.249679263683027e-05
  //   -7.768896407219534e-03 1.415004964214564e-02 6.132486209437444e-03
  // 8.831891841530757e-01 4.672704590252861e-01 5.495494782733364e-03
  // -4.839153999380735e-02 8.985770771166410e-01 -8.173884161897890e-02

  State<double> earth_state_helio = {8.968485078526299e-01,  4.333538985868564e-01,
                                     -5.249679263683027e-05, -7.768896407219534e-03,
                                     1.415004964214564e-02,  6.132486209437444e-03};
  State<double> asteroid_state_helio = {8.831891841530757e-01, 4.672704590252861e-01,
                                        5.495494782733364e-03, -4.839153999380735e-02,
                                        8.985770771166410e-01, -8.173884161897890e-02};
  State<double> earth_state_nd;
  State<double> asteroid_state_nd;

  // ---座標変換---
  earth_state_nd =
      crtbp::ConvertInertial2Rotating____(earth_state_helio, earth_state_helio, astro_params);
  asteroid_state_nd =
      crtbp::ConvertInertial2Rotating____(asteroid_state_helio, earth_state_helio, astro_params);

  std::cout << " ------- earth ------" << std::endl;
  std::cout << "state(before) : " << earth_state_helio.x << " " << earth_state_helio.y << " "
            << earth_state_helio.z << " " << earth_state_helio.vx << " " << earth_state_helio.vy
            << " " << earth_state_helio.vz << std::endl;
  std::cout << "state(after) : " << earth_state_nd.x << " " << earth_state_nd.y << " "
            << earth_state_nd.z << " " << earth_state_nd.vx << " " << earth_state_nd.vy << " "
            << earth_state_nd.vz << std::endl;
  std::cout << "velo(before:au/day) "
            << std::sqrt(earth_state_helio.vx * earth_state_helio.vx +
                         earth_state_helio.vy * earth_state_helio.vy +
                         earth_state_helio.vz * earth_state_helio.vz)
            << std::endl;
  std::cout << "velo(after:au/TU) "
            << std::sqrt(earth_state_nd.vx * earth_state_nd.vx +
                         earth_state_nd.vy * earth_state_nd.vy +
                         earth_state_nd.vz * earth_state_nd.vz)
            << std::endl;
  std::cout << "velo(after:au/day) "
            << std::sqrt(earth_state_nd.vx * earth_state_nd.vx +
                         earth_state_nd.vy * earth_state_nd.vy +
                         earth_state_nd.vz * earth_state_nd.vz) /
                   365.25
            << std::endl;

  return 0;
}