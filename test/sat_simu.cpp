// Copyright [Year] [Your Name or Organization]
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Original MATLAB script description:
// % 250305
// % effects on relative orbit from perturbations (J2 and drag)
// % relative orbit in RT plane test
// % multiple chasers in one same closed ellipse relative orbit
//
// % reference: autonomous formation flying in low earth orbit

#include <cmath>
#include <vector>
#include <iostream> // For potential output (not used for plotting)
#include <numbers> // For std::numbers::pi (C++20)
#include <numeric> // For std::atan2
#include <iomanip> // For potential output formatting
#include <fstream> // For file output (std::ofstream)

int main() {
  // Natural parameters
  const double kMu = 3.986e14;         // Gravitational parameter of Earth [m^3/s^2]
  const double kEarthRadius = 6371.0 * 1000.0; // Earth radius [m]
  const double kJ2 = 1.08263e-3;       // J2 perturbation coefficient

  // Set target orbital elements
  const double kHeightKm = 500.0;                 // Height [km]
  const double target_a = kEarthRadius + kHeightKm * 1000.0; // Semi major axis [m]
  const double target_e = 0.0;                    // Eccentricity
  const double target_i_deg = 97.5;               // Inclination [deg]
  const double target_omega_big_deg = 10.0;       // RAAN [deg]
  const double target_omega_small_deg = 0.0;      // Argument of perigee [deg]
  const double target_m_deg = 0.0;                // Mean anomaly [deg]

  // Convert degrees to radians
  const double target_i_rad = target_i_deg * std::numbers::pi / 180.0;
  const double target_omega_big_rad = target_omega_big_deg * std::numbers::pi / 180.0;
  const double target_omega_small_rad = target_omega_small_deg * std::numbers::pi / 180.0;
  const double target_m_rad = target_m_deg * std::numbers::pi / 180.0;

  // Calculate initial true longitude and mean motion
  const double u0_rad = target_omega_small_rad + target_m_rad; // Initial true longitude [rad]
  const double n_rad_per_sec = std::sqrt(kMu / std::pow(target_a, 3)); // Mean motion [rad/s]

  // Simulation time
  const double t_end_s = 10000.0; // End time [s]
  const double u_end_rad = u0_rad + n_rad_per_sec * t_end_s;

  // Generate u history (true longitude)
  std::vector<double> u_history_rad;
  double current_u = u0_rad;
  // Add a small epsilon for safety in floating point comparison
  const double u_epsilon = n_rad_per_sec * 1e-9;
  while (current_u <= u_end_rad + u_epsilon) {
    u_history_rad.push_back(current_u);
    current_u += n_rad_per_sec; // Corresponds to MATLAB's u = u0:n:u_end
                                // Note: step size is n [rad/s], not dt [s].
                                // This assumes steps in 'u' are equal to 'n'.
                                // A time-based loop might be more conventional:
                                // double dt = 1.0; // example time step
                                // for (double t = 0; t <= t_end_s; t += dt) {
                                //    u_history_rad.push_back(u0_rad + n_rad_per_sec * t);
                                // }
                                // Sticking to the original MATLAB logic for direct translation.
  }
  const size_t num_points = u_history_rad.size();

  // Target ballistic coefficient (defined but not used in provided calculations)
  const double b_target = 2.2 * 0.01 / 4.0;

  // --- Chaser 1 ---
  // Chaser 1 Absolute Orbital Elements (AOEs)
  const double chaser1_a = kEarthRadius + kHeightKm * 1000.0;
  const double chaser1_e = 1.0e-5;
  const double chaser1_i_deg = 97.8;
  const double chaser1_omega_big_deg = 10.0;
  const double chaser1_omega_small_deg = 10.0;
  const double chaser1_m_deg = 0.0;

  // Convert degrees to radians for Chaser 1
  const double chaser1_i_rad = chaser1_i_deg * std::numbers::pi / 180.0;
  const double chaser1_omega_big_rad = chaser1_omega_big_deg * std::numbers::pi / 180.0;
  const double chaser1_omega_small_rad = chaser1_omega_small_deg * std::numbers::pi / 180.0;
  // const double chaser1_m_rad = chaser1_m_deg * std::numbers::pi / 180.0; // Not used directly

  // Initial delta orbital elements for Chaser 1
  const double d_a_1_ini = chaser1_a - target_a; // Should be 0 based on values
  const double d_e_1_ini = chaser1_e - target_e;
  const double d_omega_big_1_ini = chaser1_omega_big_rad - target_omega_big_rad; // Should be 0

  // Initial Relative Orbital Elements (ROEs) for Chaser 1
  const double d_lam_1_ini = 0.0; // delta mean longitude [rad] (temporarily set)
  const double d_ex_1_ini = d_e_1_ini * std::cos(chaser1_omega_small_rad) -
                            target_e * std::cos(target_omega_small_rad);
  const double d_ey_1_ini = d_e_1_ini * std::sin(chaser1_omega_small_rad) -
                            target_e * std::sin(target_omega_small_rad);
  const double phi1_ini_rad = std::atan2(d_ey_1_ini, d_ex_1_ini);
  const double d_ix_1_ini = chaser1_i_rad - target_i_rad; // delta inclination x [rad]
  const double d_iy_1_ini = d_omega_big_1_ini * std::sin(target_i_rad); // delta inclination y [rad]

  // Drag parameters for Chaser 1 (defined but not used in calculations)
  const double b_chaser1 = 2.2 * 0.02 / 4.0;
  const double delta_b1 = b_chaser1 - b_target;

  // J2 effect parameters
  const double ita = 1.0; // sqrt(1 - e^2), for target e=0, ita=1
  const double gamma = kJ2 / 2.0 * std::pow(kEarthRadius / target_a, 2) / std::pow(ita, 4);
  // Rate of change of phi due to J2 (relative argument of perigee rate)
  const double phi1_dot_rad_per_u =
      (3.0 / 2.0) * gamma * (5.0 * std::pow(std::cos(target_i_rad), 2) - 1.0);

  // Propagate ROEs for Chaser 1 and calculate RTN coordinates
  std::vector<double> d_a_1(num_points);
  std::vector<double> d_lam_1(num_points);
  std::vector<double> d_ex_1(num_points);
  std::vector<double> d_ey_1(num_points);
  std::vector<double> d_ix_1(num_points);
  std::vector<double> d_iy_1(num_points);

  std::vector<double> r1(num_points);
  std::vector<double> t1(num_points);
  std::vector<double> n1(num_points);

  std::vector<double> d_e1_mag(num_points);
  std::vector<double> d_i1_mag(num_points);
  std::vector<double> phi1_angle(num_points);
  std::vector<double> theta1_angle(num_points);


  for (size_t k = 0; k < num_points; ++k) {
    double u_k = u_history_rad[k];
    double delta_u = u_k - u0_rad; // Difference in true longitude from start

    // Propagate ROEs according to J2 secular effects
    d_a_1[k] = d_a_1_ini; // d_a is constant in this model
    d_lam_1[k] = d_lam_1_ini -
                 (21.0 / 2.0) * (gamma * std::sin(2.0 * target_i_rad) * d_ix_1_ini +
                                 (1.0 / 7.0) * d_a_1_ini) * delta_u;
    double phi1_k = phi1_ini_rad + phi1_dot_rad_per_u * delta_u;
    d_ex_1[k] = d_e_1_ini * std::cos(phi1_k);
    d_ey_1[k] = d_e_1_ini * std::sin(phi1_k);
    d_ix_1[k] = d_ix_1_ini; // d_ix is constant in this model
    d_iy_1[k] = d_iy_1_ini +
                 3.0 * gamma * std::pow(std::sin(target_i_rad), 2) * d_ix_1_ini * delta_u;

    // Calculate relative trajectory in RTN frame
    r1[k] = (-d_ex_1[k] * std::cos(u_k) - d_ey_1[k] * std::sin(u_k)) * target_a;
    t1[k] = (d_lam_1[k] + 2.0 * d_ex_1[k] * std::sin(u_k) -
             2.0 * d_ey_1[k] * std::cos(u_k)) * target_a;
    n1[k] = (d_ix_1[k] * std::sin(u_k) - d_iy_1[k] * std::cos(u_k)) * target_a;

    // Calculate derived ROE magnitudes and angles
    d_e1_mag[k] = std::sqrt(std::pow(d_ex_1[k], 2) + std::pow(d_ey_1[k], 2));
    d_i1_mag[k] = std::sqrt(std::pow(d_ix_1[k], 2) + std::pow(d_iy_1[k], 2));
    phi1_angle[k] = std::atan2(d_ey_1[k], d_ex_1[k]);
    theta1_angle[k] = std::atan2(d_iy_1[k], d_ix_1[k]);
  }

  // --- Output / Visualization Placeholder ---
  // The MATLAB code plots the results. Here, we just print the final state
  // as an example. For plotting, results need to be written to a file
  // or passed to a plotting library.


  std::ofstream output_file("relative_orbit_output.csv");
  if(!output_file){
    std::cerr << "Error opening output file!" << std::endl;
    return 1;
  }

  std::cout << std::fixed << std::setprecision(6); // Set output precision

    output_file << "Time(s)," // Time in seconds
                << "R1(m),T1(m),N1(m)\n";

    for (size_t k = 0; k < num_points; ++k) {
        // Output each point to the file
        output_file << k * n_rad_per_sec << "," // Time in seconds
                    // << d_a_1[k] << ","
                    // << d_lam_1[k] << ","
                    // << d_ex_1[k] << ","
                    // << d_ey_1[k] << ","
                    // << d_ix_1[k] << ","
                    // << d_iy_1[k] << ","
                    << r1[k] << ","
                    << t1[k] << ","
                    << n1[k] << ","
                    // << d_e1_mag[k] << ","
                    // << d_i1_mag[k] << ","
                    // << phi1_angle[k] << ","
                    // << theta1_angle[k] << ","
                    // // Chaser 2 data follows:
                    // << d_a_2[k] << ","
                    // << d_lam_2[k] << ","
                    // << d_ex_2[k] << ","
                    // << d_ey_2[k] << ","
                    // << d_ix_2[k] << ","
                    // << d_iy_2[k] << ","
                    // << r2[k] << ","
                    // << t2[k] << ","
                    // << n2[k] << ","
                    // << d_e2_mag[k] << ","
                    // << d_i2_mag[k] << ","
                    // << phi2_angle[k]  << ","
                    // << theta2_angle[k]
                // End of line for this point.
                // Note: No newline at the end of the last line to avoid empty line in CSV.
                // output_file.flush(); // Optional: flush after each write for large files.
                // std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Optional: slow down output for readability.
                // std::cout.flush(); // Optional: flush stdout after each write for real-time output.
                // std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Optional: slow down output for readability.
                ;

  if (num_points > 0) {
    size_t last_idx = num_points - 1;
    std::cout << "Simulation complete. Number of points: " << num_points << std::endl;
    std::cout << "Final state (k = " << last_idx << ", u = " << u_history_rad[last_idx] << " rad):" << std::endl;

    std::cout << "\nChaser 1:" << std::endl;
    std::cout << "  RTN: R=" << r1[last_idx] << " m, T=" << t1[last_idx] << " m, N=" << n1[last_idx] << " m" << std::endl;


    std::cout << "No simulation points generated." << std::endl;
  }

  return 0;
}
}
