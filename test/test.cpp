#include <array>
#include <boost/numeric/odeint.hpp>
#include <boost/range/algorithm.hpp>
#include <iostream>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::array<double, 2> state_type;

typedef runge_kutta_cash_karp54<state_type> error_stepper_type;

class Myode {
   public:
    void operator()(const state_type& x,
                    state_type& dxdt,
                    const double t) {  // ばねの運動方程式
        dxdt[0] = x[1];
        dxdt[1] = -x[0];
    }
    static bool done(const state_type& x) { return x[0] <= 0.0; }
};

int main(int argc, char* argv[]) {
    Myode ode;
    state_type x = {1.0, 0.0};
    controlled_runge_kutta<runge_kutta_fehlberg78<state_type>> stepper;
    auto iter = boost::find_if(
        make_adaptive_range(stepper, ode, x, 0.0, 20.0, 0.001), Myode::done);

    cout << "integration stopped at" << iter->first << " with x = " << iter->second[0] << endl;
    return 1;
}