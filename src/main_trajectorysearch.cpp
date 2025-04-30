#include <omp.h>

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

double potential(double x, double y, double z, double mu) {
    auto calc_r1 = [&x, &y, &z, &mu]() -> double {
        return std::sqrt(std::pow(x + mu, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
    };
    auto calc_r2 = [&x, &y, &z, &mu]() -> double {
        return std::sqrt(std::pow(x - 1. + mu, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
    };

    double r1 = calc_r1();
    double r2 = calc_r2();

    return 0.5 * (x * x + y * y) + (1. - mu) / r1 + mu / r2;
}

std::vector<Point3D> create_debugmesh(const double ROI_radius,
                                      const int divisions,
                                      const Point3D& center) {
    std::vector<Point3D> meshPoints;
    if (divisions <= 0) return meshPoints;

    double step = (2.0 * ROI_radius) / (divisions - 1);
    double radiusSquared = ROI_radius * ROI_radius;

    for (int i = 0; i < divisions; ++i) {
        double x = center.x - ROI_radius + i * step;
        double y = center.y - ROI_radius + i * step;
        double z = center.z - ROI_radius + i * step;

        double x_ = std::sqrt(radiusSquared - y * y);
        double y_ = std::sqrt(radiusSquared - z * z);
        double z_ = std::sqrt(radiusSquared - x * x);

        meshPoints.emplace_back(x + center.x, center.y, center.z);
        meshPoints.emplace_back(center.x, y + center.y, center.z);
        meshPoints.emplace_back(center.x, center.y, z + center.z);

        meshPoints.emplace_back(x_ + center.x, y + center.y, center.z);
        meshPoints.emplace_back(-x_ + center.x, y + center.y, center.z);
        meshPoints.emplace_back(center.x, y_ + center.y, z + center.z);
        meshPoints.emplace_back(center.x, -y_ + center.y, z + center.z);
        meshPoints.emplace_back(x + center.x, center.y, z_ + center.z);
        meshPoints.emplace_back(x + center.x, center.y, -z_ + center.z);
    }

    std::sort(meshPoints.begin(), meshPoints.end(), [](const Point3D& a, const Point3D& b) {
        if (a.z != b.z) return a.z < b.z;
        if (a.y != b.y) return a.y < b.y;
        return a.x < b.x;
    });

    return meshPoints;
}

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
double calc_SALI(const std::array<double, 6>& ref_state,
                 const std::array<double, 6>& perturbed_state1,
                 const std::array<double, 6>& perturbed_state2) {
    std::array<double, 3> q_ref{ref_state[0], ref_state[1], ref_state[2]};
    std::array<double, 3> q_pertubed1{
        perturbed_state1[0], perturbed_state1[1], perturbed_state1[2]};
    std::array<double, 3> q_pertubed2{
        perturbed_state2[0], perturbed_state2[1], perturbed_state2[2]};
    std::array<double, 3> p_ref{
        ref_state[3] - ref_state[0], ref_state[4] + ref_state[1], ref_state[5]};
    std::array<double, 3> p_pertubed1{perturbed_state1[3] - perturbed_state1[0],
                                      perturbed_state1[4] + perturbed_state1[1],
                                      perturbed_state1[5]};
    std::array<double, 3> p_pertubed2{perturbed_state2[3] - perturbed_state2[0],
                                      perturbed_state2[4] + perturbed_state2[1],
                                      perturbed_state2[5]};
    Vector3d deviation_vec1{perturbed_state1[0] - ref_state[0],
                            perturbed_state1[1] - ref_state[1],
                            perturbed_state1[2] - ref_state[2]};
    Vector3d deviation_vec2{perturbed_state2[0] - ref_state[0],
                            perturbed_state2[1] - ref_state[1],
                            perturbed_state2[2] - ref_state[2]};

    Vector3d normalized_dev_vec1 = deviation_vec1.normalise();
    Vector3d normalized_dev_vec2 = deviation_vec2.normalise();

    Vector3d sa = normalized_dev_vec1 - normalized_dev_vec2;
    Vector3d wa = normalized_dev_vec1 + normalized_dev_vec2;

    double sa_norm = sa.magnitude();
    double wa_norm = wa.magnitude();

    // SALIの計算
    return std::min(sa_norm, wa_norm);
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
    double CALC_TIMESTEP = 0;
    double SALI_CALCTIME_THRESHOLD = 0;
    double SOI_RADIUS = 0;
    double FOREBIDDEN_AREA_RADIUS = 0;
    double JACOBI_INTEGRAL = 0;
    double inclination = 0;
    double OMEGA = 0;
    double THETA = 0;

    std::cout << std::setprecision(15);
    while (std::getline(ifs, str)) {
        if (str.find("CALC TIMESTEP") != std::string::npos) {
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

    __KEYWAIT__
    constexpr double mu = 3.003e-6;
    constexpr double eps = 1e-6;

    std::vector<Point3D> zvcPoints;
    std::vector<Point3D> meshPoints = createSphereMesh(0.01, 100, Point3D(1 - mu, 0, 0));


    for (const auto& point : meshPoints) {
        double zvp = 2 * potential(point.x, point.y, point.z, mu) - JACOBI_INTEGRAL;
        if (zvp < eps&& zvp > 0) {
            zvcPoints.push_back({point.x, point.y, point.z});
        }
    }

    double progress = 0;
    int totalIterations = zvcPoints.size();
    int count = 0;
    for (const auto& point : zvcPoints) {
        if (calc_r2(point, mu) < FOREBIDDEN_AREA_RADIUS) continue;
        displayProgressBar(progress);
        double v_abs = calc_v_abs(point, mu, JACOBI_INTEGRAL);
        double perturbation = 1e-10;
        Vector3d velocity = calc_velocity(point, v_abs, mu, inclination, OMEGA, THETA);
        if (v_abs <= 0) continue;
        double vx, vy, vz;
        vx = velocity.x();
        vy = velocity.y();
        vz = velocity.z();

        // CRTBP crtbp(point.x, point.y, point.z, velocity.x(), velocity.y(), velocity.z(), mu);

        // std::vector<std::array<double, 8>> trajectory =
        //     crtbp.trajectory_propagation(SALI_CALCTIME_THRESHOLD);

        EoM_crtbp system;
        EoM_crtbp::state ref = {point.x, point.y, point.z, vx, vy, vz};
        EoM_crtbp::state x_perturbed{point.x + perturbation, point.y, point.z, vx, vy, vz};
        EoM_crtbp::state y_perturbed{point.x, point.y + perturbation, point.z, vx, vy, vz};

        std::vector<std::array<double, 8>> out_ref;
        std::vector<std::array<double, 8>> out_xperturbed;
        std::vector<std::array<double, 8>> out_yperturbed;

        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_fehlberg78<EoM_crtbp::state>>
            stepper;

        boost::numeric::odeint::integrate_const(stepper,
                                                system,
                                                ref,
                                                0.0,
                                                SALI_CALCTIME_THRESHOLD,
                                                CALC_TIMESTEP,
                                                observer(out_ref, JACOBI_INTEGRAL));
        boost::numeric::odeint::integrate_const(stepper,
                                                system,
                                                x_perturbed,
                                                0.0,
                                                SALI_CALCTIME_THRESHOLD,
                                                CALC_TIMESTEP,
                                                observer(out_xperturbed, JACOBI_INTEGRAL));
        boost::numeric::odeint::integrate_const(stepper,
                                                system,
                                                y_perturbed,
                                                0.0,
                                                SALI_CALCTIME_THRESHOLD,
                                                CALC_TIMESTEP,
                                                observer(out_yperturbed, JACOBI_INTEGRAL));

        std::vector<double> SALI_data;
        SALI_data.reserve(out_ref.size());
        for (size_t i = 0; i < out_ref.size(); i++) {
            std::array<double, 6> ref_ = {out_ref[i][1],
                                          out_ref[i][2],
                                          out_ref[i][3],
                                          out_ref[i][4],
                                          out_ref[i][5],
                                          out_ref[i][6]};
            std::array<double, 6> x_perturbed_ = {out_xperturbed[i][1],
                                                  out_xperturbed[i][2],
                                                  out_xperturbed[i][3],
                                                  out_xperturbed[i][4],
                                                  out_xperturbed[i][5],
                                                  out_xperturbed[i][6]};
            std::array<double, 6> y_perturbed_ = {out_yperturbed[i][1],
                                                  out_yperturbed[i][2],
                                                  out_yperturbed[i][3],
                                                  out_yperturbed[i][4],
                                                  out_yperturbed[i][5],
                                                  out_yperturbed[i][6]};
            double SALIxy = calc_SALI(ref_, x_perturbed_, y_perturbed_);
            SALI_data.push_back(SALIxy);
        }
        // file output
        std::string filename = "results/debug/" + std::to_string(count) + ".dat";
        std::ofstream ofs1(filename);
        if (!ofs1) {
            std::cerr << "Failed to open file." << std::endl;
            return -1;
        }
        ofs1 << "t, x, y, z, vx, vy, vz, jacobi integral difference, SALI" << std::endl;
        ofs1 << std::fixed << std::setprecision(15);
        int countt = 0;
        for (const auto& data : out_ref) {
            ofs1 << data[0] << ", " << data[1] << ", " << data[2] << ", " << data[3] << ", "
                 << data[4] << ", " << data[5] << ", " << data[6] << ", " << data[7] << ", "
                 << SALI_data[countt] << std::endl;
            countt++;
        }

        // ofs1 << std::endl << std::endl;

        // for (const auto& data : out_xperturbed) {
        //     ofs1 << data[0] << ", " << data[1] << ", " << data[2] << ", " << data[3] << ", "
        //          << data[4] << ", " << data[5] << ", " << data[6] << ", " << data[7] <<
        //          std::endl;
        // }

        // ofs1 << std::endl << std::endl;

        // for (const auto& data : out_yperturbed) {
        //     ofs1 << data[0] << ", " << data[1] << ", " << data[2] << ", " << data[3] << ", "
        //          << data[4] << ", " << data[5] << ", " << data[6] << ", " << data[7] <<
        //          std::endl;
        // }

        ofs1.close();

        progress = (count + 1.0) / totalIterations;
        displayProgressBar(progress);
        count++;
    }

    return 0;
}