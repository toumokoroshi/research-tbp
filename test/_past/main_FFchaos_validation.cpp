#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <functional> // for std::function
#include <limits>    // for std::numeric_limits

#include <Eigen/Dense>

// Google C++ Style Guideに準拠した命名規則等を使用

// --- 定数 ---
namespace constants {
// 地球重力定数 (m^3/s^2) WGS84
constexpr double kMu = 3.986004418e14;
// 地球半径 (m) WGS84
constexpr double kRe = 6378137.0;
// J2 帯状調和係数 (無次元) WGS84
constexpr double kJ2 = 1.08262668e-3;
// 円周率
constexpr double kPi = M_PI;
// 度からラジアンへの変換係数
constexpr double kDegToRad = kPi / 180.0;
} // namespace constants

// --- 構造体 ---

// 軌道要素 (ケプラー要素)
struct OrbitalElements {
    double semi_major_axis;      // a (m) 軌道長半径
    double eccentricity;         // e (無次元) 離心率
    double inclination;          // i (rad) 軌道傾斜角
    double raan;                 // Omega (rad) 昇交点赤経
    double argument_of_perigee; // omega (rad) 近地点引数
    double true_anomaly;         // nu (rad) 真近点角
};

// 衛星の状態ベクトル (ECI: Earth-Centered Inertial frame)
struct SatelliteState {
    Eigen::Vector3d position; // r (m) 位置ベクトル
    Eigen::Vector3d velocity; // v (m/s) 速度ベクトル
};

// --- 座標変換関数 ---
namespace coordinate_transforms {

// ケプラー要素からECI状態ベクトルへ変換
// Input: oe (OrbitalElements)
// Output: SatelliteState in ECI frame
SatelliteState keplerian_to_eci(const OrbitalElements& oe) {
    using namespace constants;
    const double a = oe.semi_major_axis;
    const double e = oe.eccentricity;
    const double i = oe.inclination;
    const double Omega = oe.raan;
    const double omega = oe.argument_of_perigee;
    const double nu = oe.true_anomaly;

    // 軌道面内の位置と速度 (PQW frame)
    // p = a * (1 - e^2) : Semi-latus rectum
    const double p = a * (1.0 - e * e);
    // r_norm = p / (1 + e*cos(nu))
    const double r_norm = p / (1.0 + e * std::cos(nu));
    // r_pqw = [r_norm*cos(nu), r_norm*sin(nu), 0]^T
    const Eigen::Vector3d r_pqw(r_norm * std::cos(nu), r_norm * std::sin(nu), 0.0);
    // v_pqw = [-sqrt(mu/p)*sin(nu), sqrt(mu/p)*(e + cos(nu)), 0]^T
    const Eigen::Vector3d v_pqw(
        -std::sqrt(kMu / p) * std::sin(nu),
        std::sqrt(kMu / p) * (e + std::cos(nu)),
        0.0
    );

    // 回転行列 (PQW -> ECI)
    // R = Rz(Omega) * Rx(i) * Rz(omega)
    Eigen::Matrix3d rot_matrix;
    rot_matrix = Eigen::AngleAxisd(Omega, Eigen::Vector3d::UnitZ())
               * Eigen::AngleAxisd(i, Eigen::Vector3d::UnitX())
               * Eigen::AngleAxisd(omega, Eigen::Vector3d::UnitZ());

    SatelliteState state;
    // r_eci = R * r_pqw
    state.position = rot_matrix * r_pqw;
    // v_eci = R * v_pqw
    state.velocity = rot_matrix * v_pqw;

    return state;
}

// ECI状態ベクトルからケプラー要素へ変換 (特異点処理は簡略化)
// Input: state (SatelliteState in ECI frame)
// Output: oe (OrbitalElements)
OrbitalElements eci_to_keplerian(const SatelliteState& state) {
    using namespace constants;
    const Eigen::Vector3d r_vec = state.position;
    const Eigen::Vector3d v_vec = state.velocity;
    const double r = r_vec.norm();
    const double v = v_vec.norm();

    // h = r × v : Angular momentum vector
    const Eigen::Vector3d h_vec = r_vec.cross(v_vec);
    const double h = h_vec.norm();
    // n = k × h : Node vector (k = [0, 0, 1])
    const Eigen::Vector3d n_vec = Eigen::Vector3d::UnitZ().cross(h_vec);
    const double n = n_vec.norm();

    // e_vec = ( (v^2 - mu/r)*r - (r·v)*v ) / mu : Eccentricity vector
    const Eigen::Vector3d e_vec = ((v * v - kMu / r) * r_vec - r_vec.dot(v_vec) * v_vec) / kMu;
    const double e = e_vec.norm();

    // xi = v^2 / 2 - mu / r : Specific mechanical energy
    const double xi = v * v / 2.0 - kMu / r;

    OrbitalElements oe;
    // a = -mu / (2*xi) for elliptical orbits (e < 1)
    if (std::abs(e - 1.0) > 1e-9) {
        oe.semi_major_axis = -kMu / (2.0 * xi);
    } else { // Parabolic (e ≈ 1)
        oe.semi_major_axis = std::numeric_limits<double>::infinity();
        // Note: For hyperbolic (e > 1), a is negative by this formula.
    }
    oe.eccentricity = e;

    // i = acos(h_z / h) : Inclination
    oe.inclination = std::acos(h_vec.z() / h);

    // Omega = acos(n_x / n) : RAAN (Longitude of Ascending Node)
    // Adjust quadrant based on n_y
    if (n > 1e-9) { // Non-equatorial
        oe.raan = std::acos(n_vec.x() / n);
        if (n_vec.y() < 0.0) {
            oe.raan = 2.0 * kPi - oe.raan;
        }
    } else { // Equatorial (node vector is undefined or zero)
        oe.raan = 0.0; // Convention: Set RAAN to 0
    }

    // omega = acos( (n · e_vec) / (n * e) ) : Argument of Perigee
    // Adjust quadrant based on e_z
    if (e > 1e-9) { // Non-circular
        if (n > 1e-9) { // Non-equatorial, non-circular
            oe.argument_of_perigee = std::acos(n_vec.dot(e_vec) / (n * e));
             if (e_vec.z() < 0.0) {
                oe.argument_of_perigee = 2.0 * kPi - oe.argument_of_perigee;
            }
        } else { // Equatorial, non-circular (use longitude of perigee)
             // Longitude of perigee = atan2(e_y, e_x) adjusted for retrograde
             oe.argument_of_perigee = std::atan2(e_vec.y(), e_vec.x());
             if (h_vec.z() < 0.0) { // Retrograde equatorial
                oe.argument_of_perigee = 2.0 * kPi - oe.argument_of_perigee;
             }
        }
    } else { // Circular (perigee is undefined)
        oe.argument_of_perigee = 0.0; // Convention: Set argument of perigee to 0
    }

    // nu = acos( (e_vec · r) / (e * r) ) : True Anomaly
    // Adjust quadrant based on radial velocity (r · v)
    if (e > 1e-9) { // Non-circular
        oe.true_anomaly = std::acos(e_vec.dot(r_vec) / (e * r));
        // If r_vec.dot(v_vec) < 0, satellite is moving towards perigee
        if (r_vec.dot(v_vec) < 0.0) {
            oe.true_anomaly = 2.0 * kPi - oe.true_anomaly;
        }
    } else { // Circular (true anomaly is undefined, use argument of latitude or true longitude)
        if (n > 1e-9) { // Circular, non-equatorial: Use Argument of Latitude u = omega + nu
             // u = acos( (n · r) / (n * r) ), adjusted by r_z
             double arg_lat = std::acos(n_vec.dot(r_vec) / (n * r));
             if (r_vec.z() < 0.0) {
                 arg_lat = 2.0 * kPi - arg_lat;
             }
             oe.true_anomaly = arg_lat; // Store argument of latitude here
        } else { // Circular, equatorial: Use True Longitude l = Omega + omega + nu
             // l = atan2(r_y, r_x), adjusted by h_z for retrograde
             double true_lon = std::atan2(r_vec.y(), r_vec.x());
             if(h_vec.z() < 0.0) { // Retrograde
                 true_lon = 2.0 * kPi - true_lon;
             }
             oe.true_anomaly = true_lon; // Store true longitude here
        }
    }

    // Normalize angles to [0, 2*pi) if desired (optional)
    auto normalize_angle = [](double angle) {
        angle = std::fmod(angle, 2.0 * kPi);
        return (angle < 0.0) ? angle + 2.0 * kPi : angle;
    };
    oe.raan = normalize_angle(oe.raan);
    oe.argument_of_perigee = normalize_angle(oe.argument_of_perigee);
    oe.true_anomaly = normalize_angle(oe.true_anomaly);

    return oe;
}


// 主衛星に対する従衛星の相対状態ベクトル (ECI) を RTN 座標系に変換
// Input: relative_pos_eci (relative position in ECI),
//        chief_pos_eci (chief position in ECI),
//        chief_vel_eci (chief velocity in ECI)
// Output: relative_pos_rtn (relative position in RTN)
Eigen::Vector3d relative_eci_to_rtn(const Eigen::Vector3d& relative_pos_eci,
                                     const Eigen::Vector3d& chief_pos_eci,
                                     const Eigen::Vector3d& chief_vel_eci) {
    // RTN 基底ベクトル (unit vectors) を計算
    // R_unit = r_chief / |r_chief|
    Eigen::Vector3d r_unit = chief_pos_eci.normalized();
    // N_unit = (r_chief × v_chief) / |r_chief × v_chief|
    Eigen::Vector3d n_unit = chief_pos_eci.cross(chief_vel_eci).normalized();
    // T_unit = N_unit × R_unit (Completes the right-handed triad)
    Eigen::Vector3d t_unit = n_unit.cross(r_unit).normalized();

    // 座標変換行列 (ECI -> RTN): Rows are the RTN basis vectors in ECI
    // R_eci_to_rtn = [ R_unit^T ]
    //                [ T_unit^T ]
    //                [ N_unit^T ]
    Eigen::Matrix3d rot_eci_to_rtn;
    rot_eci_to_rtn.row(0) = r_unit; // R direction
    rot_eci_to_rtn.row(1) = t_unit; // T direction
    rot_eci_to_rtn.row(2) = n_unit; // N direction

    // 相対位置を RTN 座標系に変換
    // r_rel_rtn = R_eci_to_rtn * r_rel_eci
    return rot_eci_to_rtn * relative_pos_eci;
}

// --- この関数は不要になったため削除 ---
// // RTN 相対速度から ECI 相対速度へ変換 (簡単のため RTN 座標系の回転は無視)
// // ... (旧 relative_rtn_vel_to_eci_approx 関数のコード) ...
// --- ここまで削除 ---

} // namespace coordinate_transforms

// --- 摂動計算関数 ---
namespace perturbations {

// J2 摂動による加速度を計算
// Input: position (position vector in ECI)
// Output: accel_j2 (J2 perturbation acceleration in ECI)
Eigen::Vector3d calculate_j2_acceleration(const Eigen::Vector3d& position) {
    using namespace constants;
    const double r = position.norm();
    const double r2 = r * r;
    const double r5 = r2 * r2 * r;
    const double x = position.x();
    const double y = position.y();
    const double z = position.z();
    const double z2 = z * z;

    // Factor = -1.5 * mu * J2 * Re^2 / r^5
    const double factor = -1.5 * kMu * kJ2 * kRe * kRe / r5;

    // J2 acceleration components
    // a_j2_x = factor * x * (1 - 5*z^2/r^2)
    // a_j2_y = factor * y * (1 - 5*z^2/r^2)
    // a_j2_z = factor * z * (3 - 5*z^2/r^2)
    Eigen::Vector3d accel_j2;
    accel_j2.x() = factor * x * (1.0 - 5.0 * z2 / r2);
    accel_j2.y() = factor * y * (1.0 - 5.0 * z2 / r2);
    accel_j2.z() = factor * z * (3.0 - 5.0 * z2 / r2);

    return accel_j2;
}

// 他の摂動（大気抵抗など）が必要な場合はここに追加

} // namespace perturbations

// --- 運動方程式 ---
namespace dynamics {

// 衛星に働く加速度を計算（主重力 + 摂動）
// Input: state (SatelliteState in ECI)
// Output: acceleration (total acceleration in ECI)
Eigen::Vector3d calculate_acceleration(const SatelliteState& state) {
    using namespace constants;
    using namespace perturbations;

    const Eigen::Vector3d r_vec = state.position;
    const double r = r_vec.norm();
    const double r3 = r * r * r;

    // 主重力項 (Two-body gravity)
    // a_gravity = -mu * r / r^3
    Eigen::Vector3d accel_gravity = -kMu * r_vec / r3;

    // J2 摂動項
    Eigen::Vector3d accel_j2 = calculate_j2_acceleration(r_vec);

    // 他の摂動項があれば追加
    // Eigen::Vector3d accel_drag = calculate_drag_acceleration(state, atmospheric_model);
    // Eigen::Vector3d accel_srp = calculate_srp_acceleration(state, sun_position);
    // Eigen::Vector3d accel_third_body = calculate_third_body_acceleration(state, moon_position, sun_position);

    // 合算加速度
    // a_total = a_gravity + a_j2 + a_drag + ...
    return accel_gravity + accel_j2; // + accel_drag + ...
}

// 運動方程式の右辺 dy/dt = f(t, y)
// y = [rx, ry, rz, vx, vy, vz]^T (state vector in ECI)
// Output: dydt = [vx, vy, vz, ax, ay, az]^T
Eigen::VectorXd ode_function(double /*t*/, const Eigen::VectorXd& y) {
    // 状態ベクトル y から位置と速度を抽出
    SatelliteState current_state;
    // r = y[0:2]
    current_state.position = y.segment<3>(0);
    // v = y[3:5]
    current_state.velocity = y.segment<3>(3);

    // 加速度を計算
    // a = calculate_acceleration(r, v)
    Eigen::Vector3d acceleration = calculate_acceleration(current_state);

    // 結果ベクトル dy/dt = [v, a]^T を作成
    Eigen::VectorXd dydt(6);
    // dydt[0:2] = v
    dydt.segment<3>(0) = current_state.velocity;
    // dydt[3:5] = a
    dydt.segment<3>(3) = acceleration;

    return dydt;
}

} // namespace dynamics

// --- 数値積分器 ---
class RungeKutta4 {
public:
    // 1ステップ積分を実行
    // y_next = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    // k1 = f(t, y)
    // k2 = f(t + dt/2, y + dt/2 * k1)
    // k3 = f(t + dt/2, y + dt/2 * k2)
    // k4 = f(t + dt, y + dt * k3)
    Eigen::VectorXd step(double t, const Eigen::VectorXd& y, double dt,
                         const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>& func) {
        Eigen::VectorXd k1 = func(t, y);
        Eigen::VectorXd k2 = func(t + 0.5 * dt, y + 0.5 * dt * k1);
        Eigen::VectorXd k3 = func(t + 0.5 * dt, y + 0.5 * dt * k2);
        Eigen::VectorXd k4 = func(t + dt, y + dt * k3);

        return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
};

// --- シミュレーションクラス ---
class OrbitSimulation {
public:
    OrbitSimulation(const SatelliteState& chief_initial_state,
                    const SatelliteState& deputy_initial_state,
                    double simulation_time, double time_step)
        : chief_state_(chief_initial_state),
          deputy_state_(deputy_initial_state),
          simulation_time_(simulation_time),
          dt_(time_step),
          current_time_(0.0) {}

    // シミュレーションを実行し、結果をファイルに出力
    void run_and_output(const std::string& output_filename) {
        std::ofstream output_file(output_filename);
        if (!output_file.is_open()) {
            std::cerr << "Error: Cannot open output file " << output_filename << std::endl;
            return;
        }

        // ヘッダーを出力
        output_file << "Time(s),RelPosX_R(m),RelPosY_T(m),RelPosZ_N(m)\n";
        output_file << std::fixed << std::setprecision(6); // 精度設定

        // 初期状態の相対位置を出力
        output_relative_state(output_file);

        // 時間積分ループ
        while (current_time_ < simulation_time_) {
            // 時間ステップが終了時間を超えないように調整
            double current_dt = dt_;
            if (current_time_ + current_dt > simulation_time_) {
                current_dt = simulation_time_ - current_time_;
            }
             // 非常に小さいステップは避ける（無限ループ防止）
             if (current_dt < 1e-9) break;

            // 数値積分器で次の状態を計算
            Eigen::VectorXd chief_y(6);
            chief_y << chief_state_.position, chief_state_.velocity;
            Eigen::VectorXd next_chief_y = integrator_.step(current_time_, chief_y, current_dt, dynamics::ode_function);
            chief_state_.position = next_chief_y.segment<3>(0);
            chief_state_.velocity = next_chief_y.segment<3>(3);

            Eigen::VectorXd deputy_y(6);
            deputy_y << deputy_state_.position, deputy_state_.velocity;
            Eigen::VectorXd next_deputy_y = integrator_.step(current_time_, deputy_y, current_dt, dynamics::ode_function);
            deputy_state_.position = next_deputy_y.segment<3>(0);
            deputy_state_.velocity = next_deputy_y.segment<3>(3);

            // 時間を更新
            current_time_ += current_dt;

            // 現在の相対状態を出力
            output_relative_state(output_file);

        }

        output_file.close();
        std::cout << "Simulation finished. Results saved to " << output_filename << std::endl;
    }

private:
    // 現在の相対状態をファイルに出力
    void output_relative_state(std::ofstream& file) {
        using namespace coordinate_transforms;
        // ECIでの相対位置
        // r_rel_eci = r_deputy_eci - r_chief_eci
        Eigen::Vector3d relative_pos_eci = deputy_state_.position - chief_state_.position;
        // RTN座標系に変換
        // r_rel_rtn = R_eci_to_rtn * r_rel_eci
        Eigen::Vector3d relative_pos_rtn = relative_eci_to_rtn(relative_pos_eci,
                                                               chief_state_.position,
                                                               chief_state_.velocity);
        // ファイルに出力 (Time, R, T, N)
        file << current_time_ << ","
             << relative_pos_rtn.x() << "," // R component
             << relative_pos_rtn.y() << "," // T component
             << relative_pos_rtn.z() << "\n"; // N component
    }

    SatelliteState chief_state_;
    SatelliteState deputy_state_;
    double simulation_time_;
    double dt_;
    double current_time_;
    RungeKutta4 integrator_;
};


// --- 太陽同期軌道の例 ---
// Calculates approximate Keplerian elements for a sun-synchronous orbit
// Input: altitude_km (altitude above Earth's surface in km)
//        local_time_asc_node (desired local time at ascending node in hours - currently unused for RAAN calculation)
// Output: sso_oe (OrbitalElements for the SSO)
OrbitalElements get_sun_synchronous_orbit(double altitude_km, double /* local_time_asc_node */ = 10.5) {
    using namespace constants;
    double altitude_m = altitude_km * 1000.0;
    OrbitalElements sso_oe;

    // a = Re + altitude
    sso_oe.semi_major_axis = kRe + altitude_m;
    // Assume near-circular orbit
    sso_oe.eccentricity = 0.001;
    // Assume argument of perigee is 0 for near-circular orbit
    sso_oe.argument_of_perigee = 0.0;
    // Start at true anomaly 0 (e.g., perigee, though meaning is weak for near-circular)
    sso_oe.true_anomaly = 0.0;

    // Calculate inclination for sun-synchronous condition
    // n = sqrt(mu / a^3) : Mean motion
    double n = std::sqrt(kMu / std::pow(sso_oe.semi_major_axis, 3));
    // p = a * (1 - e^2) : Semi-latus rectum
    double p = sso_oe.semi_major_axis * (1.0 - sso_oe.eccentricity * sso_oe.eccentricity);
    // Omega_dot_target = 2*pi / (sidereal_year_in_seconds) ≈ 2*pi / (365.256 * 24 * 3600)
    // More commonly approximated using mean solar days: 2*pi / (365.24219 * 24 * 3600) rad/s
    double raan_rate_target = 2.0 * kPi / (365.24219 * 24.0 * 3600.0);

    // J2-induced RAAN precession rate: Omega_dot = -1.5 * n * J2 * (Re/p)^2 * cos(i)
    // Solve for cos(i): cos(i) = - Omega_dot_target / (1.5 * n * J2 * (Re/p)^2)
    double cos_i = -raan_rate_target / (1.5 * n * kJ2 * std::pow(kRe / p, 2));

    // Check if a valid inclination exists
    if (std::abs(cos_i) > 1.0) {
        std::cerr << "Warning: Cannot achieve sun-synchronous orbit for the given altitude. "
                  << "Calculated cos(i) = " << cos_i << ". Using a default inclination." << std::endl;
        // Use a typical SSO inclination as fallback
        sso_oe.inclination = 97.6 * kDegToRad; // Example for ~550km
    } else {
        // i = acos(cos_i)
        sso_oe.inclination = std::acos(cos_i);
    }

    // RAAN depends on the desired local time of ascending node (LTAN) and the date/time.
    // Calculating the exact RAAN requires knowing the sun's position relative to the vernal equinox.
    // Here, we simplify and set RAAN = 0. For a specific LTAN, this needs adjustment.
    sso_oe.raan = 0.0 * kDegToRad;

    std::cout << "--- Sun-Synchronous Orbit Parameters (approx) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Altitude: " << altitude_km << " km" << std::endl;
    std::cout << std::setprecision(0);
    std::cout << "Semi-major axis (a): " << sso_oe.semi_major_axis << " m" << std::endl;
    std::cout << std::setprecision(5);
    std::cout << "Eccentricity (e): " << sso_oe.eccentricity << std::endl;
    std::cout << std::setprecision(3);
    std::cout << "Inclination (i): " << sso_oe.inclination / kDegToRad << " deg" << std::endl;
    std::cout << "RAAN (Omega): " << sso_oe.raan / kDegToRad << " deg (Set to 0, LTAN not used)" << std::endl;
    std::cout << "Arg Perigee (omega): " << sso_oe.argument_of_perigee / kDegToRad << " deg" << std::endl;
    std::cout << "True Anomaly (nu): " << sso_oe.true_anomaly / kDegToRad << " deg" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;


    return sso_oe;
}

// --- メイン関数 ---
int main() {
    using namespace constants;
    using namespace coordinate_transforms;

    // --- シミュレーション設定 ---
    const double kSimulationAltitudeKm = 550.0; // 地球低軌道高度 [km]
    // 約2周回 (Period ≈ 2*pi*sqrt(a^3/mu))
    const double semi_major_axis_m = kRe + kSimulationAltitudeKm * 1000.0;
    const double orbital_period_s = 2.0 * kPi * std::sqrt(std::pow(semi_major_axis_m, 3) / kMu);
    const double kSimulationDuration = 2.0 * orbital_period_s; // シミュレーション時間 [s]
    const double kTimeStep = 1.0;             // 時間刻み [s]
    const std::string kOutputFilename = "relative_orbit_rtn_corrected.csv";

    std::cout << "Simulation Altitude: " << kSimulationAltitudeKm << " km" << std::endl;
    std::cout << "Approx Orbital Period: " << orbital_period_s / 60.0 << " min" << std::endl;
    std::cout << "Simulation Duration: " << kSimulationDuration / 60.0 << " min (" << kSimulationDuration << " s)" << std::endl;
    std::cout << "Time Step: " << kTimeStep << " s" << std::endl;


    // --- 主衛星の初期値設定 ---
    // 例: 高度550kmの太陽同期軌道
    OrbitalElements chief_oe = get_sun_synchronous_orbit(kSimulationAltitudeKm);

    // 主衛星の初期状態 (ECI) を計算
    SatelliteState chief_initial_state = keplerian_to_eci(chief_oe);

    std::cout << "\n--- Chief Satellite Initial State (ECI) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Position (m): [" << chief_initial_state.position.x() << ", "
              << chief_initial_state.position.y() << ", " << chief_initial_state.position.z() << "]" << std::endl;
    std::cout << "Velocity (m/s): [" << chief_initial_state.velocity.x() << ", "
              << chief_initial_state.velocity.y() << ", " << chief_initial_state.velocity.z() << "]" << std::endl;
    std::cout << "------------------------------------------" << std::endl;


    // --- 従衛星の初期値設定 ---
    // 主衛星に対する相対位置 (RTN) [m]
    Eigen::Vector3d relative_position_rtn(0, 10, 0); // 例: 少し前、少し上、少し右
    // 主衛星に対する相対速度 (RTN) で指定 [m/s]
    // R: Radial, T: Tangential (along-track), N: Normal (cross-track)
    Eigen::Vector3d relative_velocity_rtn(0.0, 0, 0); // 例: わずかに外向き、少し速く、わずかに北向きへ

    // ----- 相対位置・速度をECIに変換 (厳密な方法) -----

    // 1. RTN基底ベクトルとRTN->ECI回転行列を計算
    // R_unit = r_chief / |r_chief|
    Eigen::Vector3d r_unit = chief_initial_state.position.normalized();
    // N_unit = (r_chief × v_chief) / |r_chief × v_chief|
    Eigen::Vector3d n_unit = chief_initial_state.position.cross(chief_initial_state.velocity).normalized();
    // T_unit = N_unit × R_unit
    Eigen::Vector3d t_unit = n_unit.cross(r_unit);
    // R_rtn_to_eci = [ R_unit, T_unit, N_unit ] (Columns are RTN basis vectors in ECI)
    Eigen::Matrix3d rot_rtn_to_eci;
    rot_rtn_to_eci.col(0) = r_unit;
    rot_rtn_to_eci.col(1) = t_unit;
    rot_rtn_to_eci.col(2) = n_unit;

    // 2. ECI座標系での初期相対位置を計算
    // r_rel_eci = R_rtn_eci * r_rel_rtn
    Eigen::Vector3d relative_pos_eci = rot_rtn_to_eci * relative_position_rtn;

    // 3. RTN座標系のECIに対する角速度ベクトルを計算
    // ω_rtn_eci = (r_chief_eci × v_chief_eci) / |r_chief_eci|^2
    // (Note: This is the instantaneous angular velocity of the RTN frame w.r.t. ECI frame)
    Eigen::Vector3d omega_rtn_eci = chief_initial_state.position.cross(chief_initial_state.velocity)
                                    / chief_initial_state.position.squaredNorm();

    // 4. 厳密なECI座標系での初期相対速度を計算 (輸送定理)
    // v_rel_eci = R_rtn_eci * v_rel_rtn + ω_rtn_eci × r_rel_eci
    Eigen::Vector3d relative_vel_eci = rot_rtn_to_eci * relative_velocity_rtn
                                      + omega_rtn_eci.cross(relative_pos_eci);

    // ----- ここまでが厳密な変換 -----

    // 従衛星の初期状態 (ECI) を計算
    SatelliteState deputy_initial_state;
    // r_deputy_eci = r_chief_eci + r_rel_eci
    deputy_initial_state.position = chief_initial_state.position + relative_pos_eci;
    // v_deputy_eci = v_chief_eci + v_rel_eci
    deputy_initial_state.velocity = chief_initial_state.velocity + relative_vel_eci;

    std::cout << "\n--- Deputy Satellite Initial State (ECI) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Initial Relative Position (RTN, m): [" << relative_position_rtn.x() << ", "
              << relative_position_rtn.y() << ", " << relative_position_rtn.z() << "]" << std::endl;
    std::cout << "Initial Relative Velocity (RTN, m/s): [" << relative_velocity_rtn.x() << ", "
              << relative_velocity_rtn.y() << ", " << relative_velocity_rtn.z() << "]" << std::endl;
    std::cout << "Initial Relative Position (ECI, m): [" << relative_pos_eci.x() << ", "
              << relative_pos_eci.y() << ", " << relative_pos_eci.z() << "]" << std::endl;
    std::cout << "Initial Relative Velocity (ECI, m/s): [" << relative_vel_eci.x() << ", "
              << relative_vel_eci.y() << ", " << relative_vel_eci.z() << "]" << std::endl;
    std::cout << "Position (m): [" << deputy_initial_state.position.x() << ", "
              << deputy_initial_state.position.y() << ", " << deputy_initial_state.position.z() << "]" << std::endl;
    std::cout << "Velocity (m/s): [" << deputy_initial_state.velocity.x() << ", "
              << deputy_initial_state.velocity.y() << ", " << deputy_initial_state.velocity.z() << "]" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;


    // --- シミュレーション実行 ---
    OrbitSimulation simulation(chief_initial_state, deputy_initial_state,
                             kSimulationDuration, kTimeStep);
    simulation.run_and_output(kOutputFilename);

    return 0;
}