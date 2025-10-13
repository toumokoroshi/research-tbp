#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

using state_type = std::vector<double>;

// ソルバーのパラメータを集約する構造体
struct SolverParameters {
    double t_start;
    double t_end;
    double dt;
    bool store_history = true;
};

// 進捗通知用のインターフェース
class ProgressCallback {
   public:
    virtual void onProgress(double progress,
                            const state_type& current_state) = 0;
    virtual ~ProgressCallback() = default;
};

class ODESolver {
   public:
    enum class SolverType { RK4 };

    // コンストラクタでパラメータのバリデーションを行う
    explicit ODESolver(const SolverParameters& params) : params_(params) {
        validateParameters();
    }

    template <typename System>
    void rk4_step(state_type& x,
                  const System& system,
                  double t,
                  double dt) const {
        const size_t size = x.size();
        state_type k1(size), k2(size), k3(size), k4(size), temp(size);

        // k1の計算
        system(x, k1, t);

        // k2の計算
        for (size_t i = 0; i < size; ++i) {
            temp[i] = x[i] + dt * k1[i] / 2.0;
        }
        system(temp, k2, t + dt / 2.0);

        // k3の計算
        for (size_t i = 0; i < size; ++i) {
            temp[i] = x[i] + dt * k2[i] / 2.0;
        }
        system(temp, k3, t + dt / 2.0);

        // k4の計算
        for (size_t i = 0; i < size; ++i) {
            temp[i] = x[i] + dt * k3[i];
        }
        system(temp, k4, t + dt);

        // 状態の更新
        for (size_t i = 0; i < size; ++i) {
            x[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }

    // メインの解法関数
    template <typename System>
    std::vector<state_type> solve(const state_type& initial_state,
                                  const System& system,
                                  SolverType solver_type,
                                  ProgressCallback* callback = nullptr) const {
        validateInitialState(initial_state);

        state_type current_state = initial_state;
        std::vector<state_type> history;
        if (params_.store_history) {
            history.push_back(initial_state);
        }

        double t = params_.t_start;
        double current_dt = params_.dt;

        try {
            while (t < params_.t_end && !system.done(current_state)) {
                switch (solver_type) {
                    case SolverType::RK4:
                        rk4_step(current_state, system, t, current_dt);
                        break;
                    default: throw std::runtime_error("Unknown solver type");
                }

                checkNumericalStability(current_state);

                if (params_.store_history) {
                    history.push_back(current_state);
                }

                t += current_dt;

                if (callback) {
                    double progress = (t - params_.t_start) /
                                      (params_.t_end - params_.t_start);
                    callback->onProgress(progress, current_state);
                }
            }
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Solver failed: ") + e.what());
        }

        return history;
    }

   private:
    SolverParameters params_;

    void validateParameters() const {
        if (params_.dt <= 0 || params_.t_end <= params_.t_start) {
            throw std::invalid_argument("Invalid time parameters");
        }
        if (params_.tolerance <= 0) {
            throw std::invalid_argument("Invalid tolerance value");
        }
    }

    void validateInitialState(const state_type& state) const {
        if (state.empty()) {
            throw std::invalid_argument("Initial state is empty");
        }
    }

    void checkNumericalStability(const state_type& state) const {
        for (const auto& value : state) {
            if (std::isnan(value) || std::isinf(value)) {
                throw std::runtime_error("Numerical instability detected");
            }
        }
    }
};

// 使用例
class SimpleSystem {
   public:
    void operator()(const state_type& x, state_type& dxdt, double t) const {
        dxdt[0] = x[1];
        dxdt[1] = -x[0];  // 単振動の方程式
    }

    bool done(const state_type& x) const {
        return false;  // 終了条件
    }
};

class SimpleCallback : public ProgressCallback {
   public:
    void onProgress(double progress, const state_type& current_state) override {
        std::cout << "Progress: " << progress * 100 << "% Complete\n";
    }
};
