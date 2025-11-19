#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <rtbp.hpp>
#include <string>
#include <utils.hpp>
#include <vector>

struct IntegratorResult {
  std::vector<double> t;
  std::vector<double> energy;
  double mean = 0.0;
};

int main() {
  using namespace param;
  using namespace crtbp;
  using namespace my_type;
  namespace fs = std::filesystem;

  const fs::path config_dir = CONFIG_DIR;
  const fs::path output_dir = fs::path(OUTPUT_DIR) / "error_compare";
  fs::create_directories(output_dir);

  const fs::path astro_param_file = config_dir / "astro_param" / "astro_param.txt";
  const AstroConstants<double> astro_params =
      utils::loadConstants<double>(astro_param_file.string());
  const double mu = astro_params.gm_earth / (astro_params.gm_earth + astro_params.gm_sun);

  // Off-equilibrium initial condition so energy drift becomes visible.
  const State<double> initial_state{0.8, 0.0, 0.0, 0.0, 0.3, 0.0};
  const double initial_jacobi = calc_jacobi_integral(initial_state, mu);
  EquationOfMotion<double> eom(astro_params);

  // step and horizon
  const double dt = 0.02;
  const int steps = 20000;  // total simulated time = 400

  auto run_integrator = [&](auto&& stepper, State<double> state) -> IntegratorResult {
    IntegratorResult r;
    r.t.reserve(static_cast<std::size_t>(steps) + 1);
    r.energy.reserve(static_cast<std::size_t>(steps) + 1);
    double time = 0.0;
    r.t.push_back(time);
    r.energy.push_back(calc_jacobi_integral(state, mu));
    for (int i = 0; i < steps; ++i) {
      state = stepper(state, time, dt);
      time += dt;
      r.t.push_back(time);
      r.energy.push_back(calc_jacobi_integral(state, mu));
    }
    double sum = 0.0;
    for (double v : r.energy) sum += v;
    r.mean = sum / static_cast<double>(r.energy.size());
    return r;
  };

  auto rk4_step = [&](const State<double>& s, double /*t*/, double h) {
    return RungeKutta4Step(eom, s, 0.0, h);
  };
  auto y4_step = [&](const State<double>& s, double /*t*/, double h) {
    return SymplecticStep4thOrder(mu, s, h);
  };
  auto bm4_step = [&](const State<double>& s, double /*t*/, double h) {
    return SymplecticStep4thOrderBM(mu, s, h);
  };
  auto y6_step = [&](const State<double>& s, double /*t*/, double h) {
    return SymplecticStep6thOrder(mu, s, h);
  };

  auto res_rk4 = run_integrator(rk4_step, initial_state);
  auto res_y4 = run_integrator(y4_step, initial_state);
  auto res_bm4 = run_integrator(bm4_step, initial_state);
  auto res_y6 = run_integrator(y6_step, initial_state);

  const fs::path csv_path = output_dir / "energy_error.csv";
  std::ofstream csv(csv_path);
  if (!csv) {
    std::cerr << "Failed to open output file: " << csv_path << std::endl;
    return 1;
  }

  csv << std::fixed << std::setprecision(15);
  csv << "time,rk4_error,yoshida4_error,bm4_error,yoshida6_error,"
         "rk4_shadow,bm4_shadow,yoshida6_shadow\n";
  for (std::size_t i = 0; i < res_rk4.t.size(); ++i) {
    const double t = res_rk4.t[i];
    const double er_rk4 = res_rk4.energy[i] - initial_jacobi;
    const double er_y4 = res_y4.energy[i] - initial_jacobi;
    const double er_bm4 = res_bm4.energy[i] - initial_jacobi;
    const double er_y6 = res_y6.energy[i] - initial_jacobi;
    const double sh_rk4 = res_rk4.energy[i] - res_rk4.mean;
    const double sh_bm4 = res_bm4.energy[i] - res_bm4.mean;
    const double sh_y6 = res_y6.energy[i] - res_y6.mean;
    csv << t << "," << er_rk4 << "," << er_y4 << "," << er_bm4 << "," << er_y6 << "," << sh_rk4
        << "," << sh_bm4 << "," << sh_y6 << "\n";
  }
  csv.close();

  const fs::path png_path = output_dir / "energy_error.png";
  const fs::path gp_path = output_dir / "energy_error.gp";

  std::ofstream gp(gp_path);
  if (gp) {
    const std::string csv_file = csv_path.generic_string();
    const std::string png_file = png_path.generic_string();

    gp << "set datafile separator ','\n";
    gp << "set term pngcairo size 1400,900\n";
    gp << "set output '" << png_file << "'\n";
    gp << "set multiplot layout 2,1 title 'Jacobi error (RK4, Yoshida4, BM4, Yoshida6)'\n";
    gp << "set grid\n";
    gp << "set key top left\n";
    gp << "set ylabel 'delta C_J'\n";
    gp << "plot \\\n";
    gp << "  '" << csv_file
       << "' using 1:2 with lines lw 2 lc rgb '#d62728' title 'RK4', \\\n";
    gp << "  '" << csv_file
       << "' using 1:3 with lines lw 2 lc rgb '#1f77b4' title 'Yoshida4', \\\n";
    gp << "  '" << csv_file
       << "' using 1:4 with lines lw 2 lc rgb '#2ca02c' title 'BM4', \\\n";
    gp << "  '" << csv_file
       << "' using 1:5 with lines lw 2 lc rgb '#9467bd' title 'Yoshida6'\n";
    gp << "set ylabel 'delta (C_J - mean)'\n";
    gp << "plot \\\n";
    gp << "  '" << csv_file
       << "' using 1:6 with lines lw 2 lc rgb '#ff7f0e' title 'RK4 shadow', \\\n";
    gp << "  '" << csv_file
       << "' using 1:7 with lines lw 2 lc rgb '#17becf' title 'BM4 shadow', \\\n";
    gp << "  '" << csv_file
       << "' using 1:8 with lines lw 2 lc rgb '#8c564b' title 'Yoshida6 shadow'\n";
    gp << "unset multiplot\n";
  }
  gp.close();

  const std::string gnuplot_cmd = "gnuplot \"" + gp_path.generic_string() + "\"";
  std::system(gnuplot_cmd.c_str());

  std::cout << "Data written to: " << csv_path << std::endl;
  std::cout << "Plot script: " << gp_path << "\nImage (if gnuplot present): " << png_path
            << std::endl;

  return 0;
}
