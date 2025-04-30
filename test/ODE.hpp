/**
 * @file ODE.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-01-30
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef ODE_HPP
#define ODE_HPP

#include <vector>

using state_type = std::vector<double>;

// ルンゲクッタ法を実装したソルバークラス
#include <vector>
#include <functional>
#include <iostream>

using state_type = std::vector<double>;

class ODESolver {
public:
    // ソルバーの種類を定義
    enum class SolverType {
        EULER,
        RK4
    };

    ODESolver(double start_time, double end_time, double dt)
        : t_start(start_time), t_end(end_time), dt(dt) {}

    // メインの解法関数
    template<typename System>
    state_type solve(const state_type& initial_state, System system, SolverType solver_type) {
        state_type current_state = initial_state;
        double t = t_start;

        while (t < t_end && !system.done(current_state)) {
            switch (solver_type) {
                case SolverType::EULER:
                    euler_step(current_state, system, t);
                    break;
                case SolverType::RK4:
                    rk4_step(current_state, system, t);
                    break;
            }
            t += dt;
        }

        return current_state;
    }

private:
    double t_start;
    double t_end;
    double dt;

    // オイラー法による1ステップの計算
    template<typename System>
    void euler_step(state_type& x, System& system, double t) {
        state_type dxdt(x.size());
        system(x, dxdt, t);
        
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += dt * dxdt[i];
        }
    }

    // RK4による1ステップの計算
    template<typename System>
    void rk4_step(state_type& x, System& system, double t) {
        state_type k1(x.size()), k2(x.size()), k3(x.size()), k4(x.size());
        state_type temp(x.size());

        // k1の計算
        system(x, k1, t);
        
        // k2の計算
        for (size_t i = 0; i < x.size(); ++i) {
            temp[i] = x[i] + dt * k1[i] / 2.0;
        }
        system(temp, k2, t + dt/2.0);

        // k3の計算
        for (size_t i = 0; i < x.size(); ++i) {
            temp[i] = x[i] + dt * k2[i] / 2.0;
        }
        system(temp, k3, t + dt/2.0);

        // k4の計算
        for (size_t i = 0; i < x.size(); ++i) {
            temp[i] = x[i] + dt * k3[i];
        }
        system(temp, k4, t + dt);

        // 状態の更新
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
        }
    }
};

// 使用例
struct HarmonicOscillator {
    void operator()(const state_type& x, state_type& dxdt, const double/*  t */) {
        dxdt[0] = x[1];          // dx/dt = v
        dxdt[1] = -x[0];         // dv/dt = -x
    }
    
    static bool done(const state_type& x) {
        return false;  // 終了条件なし
    }
};

int main() {
    // 初期条件の設定
    state_type initial_state = {1.0, 0.0};  // 位置=1.0, 速度=0.0
    
    // ソルバーの初期化
    ODESolver solver(0.0, 10.0, 0.01);  // 開始時間=0.0, 終了時間=10.0, dt=0.01
    
    // システムの定義
    HarmonicOscillator system;
    
    // オイラー法で解く
    auto result_euler = solver.solve(initial_state, system, ODESolver::SolverType::EULER);
    
    // RK4で解く
    auto result_rk4 = solver.solve(initial_state, system, ODESolver::SolverType::RK4);
    
    return 0;
}
