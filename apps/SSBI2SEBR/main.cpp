
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <rtbp.hpp>
#include <utils.hpp>
#include <vector3d.hpp>

// SSB(solar system barycenter) to SEBR(sun-earth barycenter rotating frame)
// conversion tool
int main() {
  using namespace rtbp::crtbp;
  std::string config_base_path = CONFIG_DIR;
  std::string astro_param_file =
      config_base_path + "/astro_param/astro_param.txt";
  AstroConstants astro_params = loadConstants(astro_param_file);
  const double kAU = astro_params.au; // astronomical unit in meters
  const double kGMSUN =
      astro_params.gm_sun; // heliocentric gravitational constant m3 s-2
  const double kGMEARTH =
      astro_params.gm_earth; // geocentric gravitational constant m3 s-2
  const double kMU = astro_params.mu; // mu parameter of Earth-Sun

  // $$SOE
  // 2460984.458333333, A.D. 2025-Nov-04
  // 23:00:00.0000,  7.304196237871223E-01,  6.098274565642665E-01,  2.646650323945715E-01,
  // -1.224114608078465E-02,  1.204646234284942E-02,  5.261882690463215E-03,
  // 2460999.458333333, A.D. 2025-Nov-19
  // 23:00:00.0000,  5.257984240277522E-01,  7.577655799250410E-01,  3.283858264463194E-01,
  // -1.433600204809335E-02,  8.127139041478151E-03,  3.503870113807898E-03,
  // $$EOE
  Vector3d target_pos_SSB_perigee{7.304196237871223E-01, 6.098274565642665E-01,
                                  2.646650323945715E-01};
  Vector3d target_vel_SSB_perigee{-1.224114608078465E-02, 1.204646234284942E-02,
                                  5.261882690463215E-03};
  Vector3d target_pos_SSB_apogee{5.257984240277522E-01, 7.577655799250410E-01,
                                 3.283858264463194E-01};
  Vector3d target_vel_SSB_apogee{-1.433600204809335E-02, 8.127139041478151E-03,
                                 3.503870113807898E-03};

  // $$SOE
  // 2460984.458333333, A.D. 2025-Nov-04
  // 23:00:00.0000,  7.284388574310334E-01,  6.086682022418486E-01,  2.639921635473907E-01,
  // -1.187458565021005E-02,  1.158622746231203E-02,  5.021849812165025E-03,
  // 2460999.458333333, A.D. 2025-Nov-19
  // 23:00:00.0000,  5.273914666583941E-01,  7.596884055466107E-01,  3.294606499098446E-01,
  // -1.478997985135310E-02,  8.422006642370750E-03,  3.650959645440828E-03,
  // $$EOE
  Vector3d p2_pos_SSB_perigee{7.284388574310334E-01, 6.086682022418486E-01,
                              2.639921635473907E-01};
  Vector3d p2_vel_SSB_perigee{-1.187458565021005E-02, 1.158622746231203E-02,
                              5.021849812165025E-03};
  Vector3d p2_pos_SSB_apogee{5.273914666583941E-01, 7.596884055466107E-01,
                             3.294606499098446E-01};
  Vector3d p2_vel_SSB_apogee{-1.478997985135310E-02, 8.422006642370750E-03,
                             3.650959645440828E-03};

  std::array<Vector3d, 2> target_state_SSB_perigee{target_pos_SSB_perigee,
                                                   target_vel_SSB_perigee};
  std::array<Vector3d, 2> target_state_SSB_apogee{target_pos_SSB_apogee,
                                                  target_vel_SSB_apogee};
  std::array<Vector3d, 2> p2_state_SSB_perigee{p2_pos_SSB_perigee,
                                               p2_vel_SSB_perigee};
  std::array<Vector3d, 2> p2_state_SSB_apogee{p2_pos_SSB_apogee,
                                              p2_vel_SSB_apogee};

  std::array<Vector3d, 2> target_state_SEBR_apogee =
      FrameTransformation<Vector3d, 2>(target_state_SSB_apogee,
                                       p2_state_SSB_apogee, astro_params);
  std::array<Vector3d, 2> target_state_SEBR_perigee =
      FrameTransformation<Vector3d, 2>(target_state_SSB_perigee,
                                       p2_state_SSB_perigee, astro_params);
  std::cout << std::setprecision(15);
  std::cout << "target pos SEBR at apogee = " << target_state_SEBR_apogee[0].x()
            << " " << target_state_SEBR_apogee[0].y() << " "
            << target_state_SEBR_apogee[0].z() << std::endl;
  std::cout << "target vel SEBR at apogee = " << target_state_SEBR_apogee[1].x()
            << " " << target_state_SEBR_apogee[1].y() << " "
            << target_state_SEBR_apogee[1].z() << std::endl;
  std::cout << "target pos SEBR at perigee = "
            << target_state_SEBR_perigee[0].x() << " "
            << target_state_SEBR_perigee[0].y() << " "
            << target_state_SEBR_perigee[0].z() << std::endl;
  std::cout << "target vel SEBR at perigee = "
            << target_state_SEBR_perigee[1].x() << " "
            << target_state_SEBR_perigee[1].y() << " "
            << target_state_SEBR_perigee[1].z() << std::endl;

  //   calc jacobi constant
  double x_apogee = target_state_SEBR_apogee[0].x();
  double y_apogee = target_state_SEBR_apogee[0].y();
  double z_apogee = target_state_SEBR_apogee[0].z();
  double vx_apogee = target_state_SEBR_apogee[1].x();
  double vy_apogee = target_state_SEBR_apogee[1].y();
  double vz_apogee = target_state_SEBR_apogee[1].z();
  double x_perigee = target_state_SEBR_perigee[0].x();
  double y_perigee = target_state_SEBR_perigee[0].y();
  double z_perigee = target_state_SEBR_perigee[0].z();
  double vx_perigee = target_state_SEBR_perigee[1].x();
  double vy_perigee = target_state_SEBR_perigee[1].y();
  double vz_perigee = target_state_SEBR_perigee[1].z();
  double r1_apogee = std::sqrt((x_apogee + kMU) * (x_apogee + kMU) +
                               y_apogee * y_apogee + z_apogee * z_apogee);
  double r1_perigee = std::sqrt((x_perigee + kMU) * (x_perigee + kMU) +
                                y_perigee * y_perigee + z_perigee * z_perigee);
  double r2_apogee = std::sqrt((x_apogee - 1 + kMU) * (x_apogee - 1 + kMU) +
                               y_apogee * y_apogee + z_apogee * z_apogee);
  double r2_perigee = std::sqrt((x_perigee - 1 + kMU) * (x_perigee - 1 + kMU) +
                                y_perigee * y_perigee + z_perigee * z_perigee);

  double jacobi_constant_apogee =
      x_apogee * x_apogee + y_apogee * y_apogee + 2 * (1 - kMU) / r1_apogee +
      2 * kMU / r2_apogee -
      (vx_apogee * vx_apogee + vy_apogee * vy_apogee + vz_apogee * vz_apogee);
  double jacobi_constant_perigee =
      x_perigee * x_perigee + y_perigee * y_perigee +
      2 * (1 - kMU) / r1_perigee + 2 * kMU / r2_perigee -
      (vx_perigee * vx_perigee + vy_perigee * vy_perigee +
       vz_perigee * vz_perigee);
  //    kmに直してプリント
  std::cout << "target pos SEBR at apogee (km) = "
            << target_state_SEBR_apogee[0].x() * kAU / 1000.0 << " "
            << target_state_SEBR_apogee[0].y() * kAU / 1000.0 << " "
            << target_state_SEBR_apogee[0].z() * kAU / 1000.0 << " "
            << std::endl;
  std::cout << "Jacobi constant at apogee = " << jacobi_constant_apogee
            << std::endl;
  std::cout << "Jacobi constant at perigee = " << jacobi_constant_perigee
            << std::endl;
  return 0;
}
