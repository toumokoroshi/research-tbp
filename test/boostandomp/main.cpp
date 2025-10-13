#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>

// OpenMPヘッダー
#include <omp.h>

// 状態ベクトルを定義 (x, v)
using state_type = std::vector<double>;
constexpr double M_PI = 3.14159265358979323846;
// 単振動の運動方程式を定義
void harmonic_oscillator(const state_type &x, state_type &dxdt, double t) {
  const double omega = 2.0 * M_PI; // 角振動数 ω = 2π (周期1秒)
  dxdt[0] = x[1];                  // dx/dt = v
  dxdt[1] = -omega * omega * x[0]; // dv/dt = -ω^2 * x
}

// 結果を保存するためのオブザーバー
struct observer {
  std::vector<double> &m_t;
  std::vector<state_type> &m_x;

  observer(std::vector<double> &t, std::vector<state_type> &x)
      : m_t(t), m_x(x) {}

  void operator()(const state_type &x, double t) {
    m_t.push_back(t);
    m_x.push_back(x);
  }
};

// シミュレーションを1回実行する関数
void run_simulation(int run_index) {
  // 初期条件: x(0) = 1.0, v(0) = 0.0
  state_type x = {1.0, 0.0};

  // 結果を格納するコンテナ
  std::vector<double> t_vec;
  std::vector<state_type> x_vec;

  // Dormand-Prince 5(4)法ソルバー
  using namespace boost::numeric::odeint;
  runge_kutta_dopri5<state_type> stepper;

  // 0秒から100秒まで、0.01秒ステップで計算
  integrate_adaptive(stepper, harmonic_oscillator, x, 0.0, 100.0, 0.01,
                     observer(t_vec, x_vec));

  // 最初の実行結果のみファイルに保存
  if (run_index == 0) {
    std::ofstream ofs("result.dat");
    for (size_t i = 0; i < t_vec.size(); ++i) {
      ofs << t_vec[i] << "\t" << x_vec[i][0] << "\t" << x_vec[i][1]
          << std::endl;
    }
  }
}

int main() {
  const int N = 1000; // シミュレーションの繰り返し回数

  // --- OpenMPなし (逐次実行) ---
  std::cout << "Running " << N << " simulations sequentially..." << std::endl;
  auto start_seq = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < N; ++i) {
    run_simulation(i);
  }

  auto end_seq = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seq = end_seq - start_seq;
  std::cout << "Sequential execution time: " << elapsed_seq.count() << " s"
            << std::endl;
  std::cout << "Average time per simulation: " << elapsed_seq.count() / N
            << " s" << std::endl;
  std::cout << "------------------------------------------" << std::endl;

  // --- OpenMPあり (並列実行) ---
  std::cout << "Running " << N << " simulations in parallel with OpenMP..."
            << std::endl;

  auto start_omp = std::chrono::high_resolution_clock::now();

// まず #pragma omp parallel で並列領域を開始する
#pragma omp parallel
  {
// マスターブランチ（通常はスレッド0）が一度だけ実行する
#pragma omp master
    {
      std::cout << "[Info] Parallel region is running with "
                << omp_get_num_threads() << " threads." << std::endl;
    }

// 次に #pragma omp for でループ処理を各スレッドに分散させる
#pragma omp for
    for (int i = 0; i < N; ++i) {
      run_simulation(i);
    }
  } // ここで並列領域が終了する

  auto end_omp = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_omp = end_omp - start_omp;
  std::cout << "Parallel execution time (OpenMP): " << elapsed_omp.count()
            << " s" << std::endl;
  std::cout << "Average time per simulation: " << elapsed_omp.count() / N
            << " s" << std::endl;
  std::cout << "------------------------------------------" << std::endl;

  double speedup = elapsed_seq.count() / elapsed_omp.count();
  std::cout << "Speedup: " << speedup << "x" << std::endl;
  std::cout << "\nSimulation data for one run has been saved to 'result.dat'."
            << std::endl;
  std::cout << "Use 'plot.plt' with gnuplot to visualize the result."
            << std::endl;

  return 0;
}