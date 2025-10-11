#include <array>
#include <boost/numeric/odeint.hpp>
#include <chrono>
#include <cmath>
#include <crtbp.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <utils.hpp>
#include <vector3d.hpp>
#include <vector>


double calc_r1(const Point3D& point, const double mu) {
    return std::sqrt(std::pow(point.x + mu, 2.) + std::pow(point.y, 2.) + std::pow(point.z, 2.));
}

double calc_r2(const Point3D& point, const double mu) {
    return std::sqrt(std::pow(point.x - 1. + mu, 2.) + std::pow(point.y, 2.) +
                     std::pow(point.z, 2.));
}

double calc_v_abs(const Point3D& point, const double mu, const double JACOBI_INTEGRAL) {
    double r1 = calc_r1(point, mu);
    double r2 = calc_r2(point, mu);
    return std::sqrt(point.x * point.x + point.y * point.y + 2. * (1. - mu) / r1 + 2. * mu / r2 +
                     mu * (1. - mu) - JACOBI_INTEGRAL);
}

double calc_jacobi_integral(const std::array<double, 6>& state, const double mu) {
    const double x = state[0];
    const double y = state[1];
    const double z = state[2];
    const double vx = state[3];
    const double vy = state[4];
    const double vz = state[5];

    // 第一質点からの距離
    const double r1 = std::sqrt(std::pow(x + mu, 2) + y * y + z * z);

    // 第二質点からの距離
    const double r2 = std::sqrt(std::pow(x - 1. + mu, 2) + y * y + z * z);

    // ヤコビ積分の計算
    return x * x + y * y - (vx * vx + vy * vy + vz * vz) + 2. * (1. - mu) / r1 + 2. * mu / r2 +
           mu * (1. - mu);
}

Vector3d calc_velocity(const Point3D& point,
                       const double v_abs,
                       const double mu,
                       const double inclination,
                       const double OMEGA,
                       const double theta) {
    // 回転軸と回転角からクオータニオン経由で回転行列を生成
    auto create_rot_matrix = [](const Vector3d& unit_n,
                                double theta_) -> std::array<std::array<double, 3>, 3> {
        double half_theta = theta_ / 2.0;
        double q0 = std::cos(half_theta);
        double sin_half_theta = std::sin(half_theta);
        double q1 = unit_n.x() * sin_half_theta;
        double q2 = unit_n.y() * sin_half_theta;
        double q3 = unit_n.z() * sin_half_theta;

        double q0q0 = q0 * q0;
        double q1q1 = q1 * q1;
        double q2q2 = q2 * q2;
        double q3q3 = q3 * q3;
        double q0q1 = q0 * q1;
        double q0q2 = q0 * q2;
        double q0q3 = q0 * q3;
        double q1q2 = q1 * q2;
        double q1q3 = q1 * q3;
        double q2q3 = q2 * q3;

        std::array<std::array<double, 3>, 3> rot_matrix = {
            {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 - q0q3), 2.0 * (q1q3 + q0q2)},
             {2.0 * (q1q2 + q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 - q0q1)},
             {2.0 * (q1q3 - q0q2), 2.0 * (q2q3 + q0q1), q0q0 - q1q1 - q2q2 + q3q3}}
        };
        return rot_matrix;
    };

    // 入力されたベクトルを入力された行列で変換
    auto convert = [](std::array<std::array<double, 3>, 3> convert_matrix,
                      const Vector3d& v) -> Vector3d {
        return {convert_matrix[0][0] * v.x() + convert_matrix[0][1] * v.y() +
                    convert_matrix[0][2] * v.z(),
                convert_matrix[1][0] * v.x() + convert_matrix[1][1] * v.y() +
                    convert_matrix[1][2] * v.z(),
                convert_matrix[2][0] * v.x() + convert_matrix[2][1] * v.y() +
                    convert_matrix[2][2] * v.z()};
    };
    double vx_, vy_, vz_ = 0;

    // inclinationとOMEGAを用いて軌道面の法線ベクトルを計算
    Vector3d normal_vector{std::sin(inclination) * std::cos(OMEGA),
                           std::sin(inclination) * std::sin(OMEGA),
                           std::cos(inclination)};
    // 法線ベクトルと位置ベクトルの外積を計算
    Vector3d r2_vector{point.x - 1. + mu, point.y, point.z};
    Vector3d h_vector = r2_vector.gaiseki(normal_vector);
    Vector3d normalized_h_vector = h_vector.normalise();

    // 速度ベクトルを計算
    vx_ = v_abs * normalized_h_vector.x();
    vy_ = v_abs * normalized_h_vector.y();
    vz_ = v_abs * normalized_h_vector.z();

    if (theta == 0.0)
        return Vector3d(vx_, vy_, vz_);
    else {
        // クオータニオンを用いて速度ベクトルをnormal_vector周りにthetaだけ回転させる
        std::array<std::array<double, 3>, 3> rot_matrix = create_rot_matrix(normal_vector, theta);
        Vector3d velocity{vx_, vy_, vz_};
        Vector3d rotated_velocity = convert(rot_matrix, velocity);

        return rotated_velocity;
    }
}

struct EoM_crtbp {
   public:
    using state = std::array<double, 6>;
    EoM_crtbp() = default;

    void operator()(const state& x, state& dxdt, const double /* t */) const {
        const double x_ = x[0];
        const double y_ = x[1];
        const double z_ = x[2];
        const double vx_ = x[3];
        const double vy_ = x[4];
        const double vz_ = x[5];

        const double r1 = std::sqrt(std::pow(x_ + 3.003e-6, 2) + y_ * y_ + z_ * z_);
        const double r2 = std::sqrt(std::pow(x_ - 1. + 3.003e-6, 2) + y_ * y_ + z_ * z_);
        dxdt[0] = vx_;
        dxdt[1] = vy_;
        dxdt[2] = vz_;
        dxdt[3] = 2. * vy_ + x_ - (1. - 3.003e-6) * (x_ + 3.003e-6) / std::pow(r1, 3) -
                  3.003e-6 * (x_ - 1. + 3.003e-6) / std::pow(r2, 3);
        dxdt[4] = -2. * vx_ + y_ - (1. - 3.003e-6) * y_ / std::pow(r1, 3) -
                  3.003e-6 * y_ / std::pow(r2, 3);
        dxdt[5] = -(1. - 3.003e-6) * z_ / std::pow(r1, 3) - 3.003e-6 * z_ / std::pow(r2, 3);
    }
};

struct observer {
   public:
    using state = EoM_crtbp::state;
    std::vector<std::array<double, 8>>& m_out;
    double init_jacobi_;
    double mu = 3.003e-6;
    observer(std::vector<std::array<double, 8>>& out, double init_jacobi)
        : m_out(out), init_jacobi_(init_jacobi) {}

    void operator()(const state& x, double t) const {
        double jacobi_integral = calc_jacobi_integral(x, mu);
        m_out.push_back({t, x[0], x[1], x[2], x[3], x[4], x[5], init_jacobi_ - jacobi_integral});
    }
};

int main() {
    std::ifstream ifs("config/point.txt");
    if (!ifs) {
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }
    std::string str;
    double x, y, z, vx, vy, vz = 0;
    double CALC_TIMESTEP = 0;
    double SALI_CALCTIME_THRESHOLD = 0;
    double SOI_RADIUS = 0;
    double FOREBIDDEN_AREA_RADIUS = 0;
    double JACOBI_INTEGRAL = 0;
    double inclination = 0;
    double OMEGA = 0;
    double THETA = 0;
    std::cout << "<> [config] select the mode" << std::endl;
    std::cout << "<>    1: calculate the initial velocity with jacobi constant and inital position"
              << std::endl;
    std::cout << "<>    2: load the initial velocity from config file" << std::endl;
    int mode;
    std::cin >> mode;

    std::cout << std::setprecision(15);
    while (std::getline(ifs, str)) {
        if (str.find("x") != std::string::npos) {
            x = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        x : " << x << std::endl;
        } else if (str.find("y") != std::string::npos) {
            y = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        y : " << y << std::endl;
        } else if (str.find("z") != std::string::npos) {
            z = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        z : " << z << std::endl;
        } else if (str.find("vx") != std::string::npos) {
            vx = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        vx : " << vx << std::endl;
        } else if (str.find("vy") != std::string::npos) {
            vy = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        vy : " << vy << std::endl;
        } else if (str.find("vz") != std::string::npos) {
            vz = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        vz : " << vz << std::endl;
        } else if (str.find("CALC TIMESTEP") != std::string::npos) {
            CALC_TIMESTEP = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        CALC TIMESTEP : " << CALC_TIMESTEP << std::endl;
        } else if (str.find("SALI CALCTIME THRESHOLD") != std::string::npos) {
            SALI_CALCTIME_THRESHOLD = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        SALI CALCTIME THRESHOLD : " << SALI_CALCTIME_THRESHOLD
                      << std::endl;
        } else if (str.find("RADIUS OF SOI") != std::string::npos) {
            SOI_RADIUS = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        SOI RADIUS : " << SOI_RADIUS << std::endl;
        } else if (str.find("RADIUS OF FOREBIDDEN AREA") != std::string::npos) {
            FOREBIDDEN_AREA_RADIUS = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        RADIUS OF FOREBIDDEN AREA : " << FOREBIDDEN_AREA_RADIUS
                      << std::endl;
        } else if (str.find("JACOBI INTEGRAL") != std::string::npos) {
            JACOBI_INTEGRAL = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        JACOBI INTEGRAL : " << JACOBI_INTEGRAL << std::endl;
        } else if (str.find("INCLINATION AGAINST XY PLANE(deg)") != std::string::npos) {
            inclination = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        INCLINATION(deg) : " << inclination << std::endl;
            inclination = inclination * std::acos(-1) / 180.;  // deg to rad
        } else if (str.find("LONGTITUDE AGAINST X AXIS+(deg)") != std::string::npos) {
            OMEGA = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        LONGTITUDE(deg) : " << OMEGA << std::endl;
            OMEGA = OMEGA * std::acos(-1) / 180.;  // deg to rad
        } else if (str.find("DEGREE FROM TANGENT") != std::string::npos) {
            THETA = std::stod(str.substr(str.find("=") + 1));
            std::cout << "<>        DEGREE FROM TANGENT(deg) : " << THETA << std::endl;
            THETA = THETA * std::acos(-1) / 180.;  // deg to rad
        }
    }
    ifs.close();
    Point3D point(x, y, z);
    double mu = 3.003e-6;
    double v_abs = calc_v_abs(point, mu, JACOBI_INTEGRAL);
    double perturbation = 1e-10;
    Vector3d velocity = calc_velocity(point, v_abs, mu, inclination, OMEGA, THETA);
    if (mode == 1) {
        vx = velocity.x();
        vy = velocity.y();
        vz = velocity.z();
    }

    // 速度と位置の内積を計算
    double inner_product = velocity.naiseki(Vector3d(point.x - 1. + mu, point.y, point.z));
    // 表示
    std::cout << "x: " << point.x << std::endl;
    std::cout << "y: " << point.y << std::endl;
    std::cout << "z: " << point.z << std::endl;
    std::cout << "vx: " << vx << std::endl;
    std::cout << "vy: " << vy << std::endl;
    std::cout << "vz: " << vz << std::endl;

    CRTBP crtbp(point.x, point.y, point.z, velocity.x(), velocity.y(), velocity.z(), mu);

    __KEYWAIT__
    std::vector<std::array<double, 8>> trajectory =
        crtbp.trajectory_propagation(SALI_CALCTIME_THRESHOLD);

    EoM_crtbp system;
    EoM_crtbp::state _x = {point.x, point.y, point.z, vx, vy, vz};
    std::vector<std::array<double, 8>> out;

    boost::numeric::odeint::controlled_runge_kutta<
        boost::numeric::odeint::runge_kutta_fehlberg78<EoM_crtbp::state>> stepper;
    boost::numeric::odeint::


    boost::numeric::odeint::integrate_const(
        stepper, system, _x, 0.0, SALI_CALCTIME_THRESHOLD, CALC_TIMESTEP, observer(out, JACOBI_INTEGRAL));

    std::string filename = "results/trajectory_" + getcurrent_date() + ".dat";
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }
    ofs << std::fixed<<std::setprecision(15);
    for (const auto& state : trajectory) {
        for (const auto& elem : state) {
            ofs << elem << " ";
        }
        ofs << std::endl;
    }
    ofs << "\n\n";
    for (const auto& state : out) {
        for (const auto& elem : state) {
            ofs << elem << " ";
        }
        ofs << std::endl;
    }

    ofs.close();

    return 0;
}