#ifndef BARYCENTRIC_CONVERSION_HPP_
#define BARYCENTRIC_CONVERSION_HPP_

#include <cmath>
#include <stdexcept>

#include "rtbp.hpp"

namespace crtbp {

/**
 * @brief Type of orbital anomaly stored in KeplerianElements::anomaly.
 */
enum class AnomalyType { kTrue, kMean, kEccentric };

/**
 * @brief Simple container for Keplerian orbital elements referenced to the SSB inertial frame.
 *
 * Angles are expected in radians. Length units must be consistent with the supplied GM
 * (e.g., meters when GM is given in m^3/s^2). Velocities derived from these elements therefore
 * come out in the corresponding length/time units (e.g., m/s).
 */
template <typename ScalarType>
struct KeplerianElements {
  ScalarType semi_major_axis;         ///< a (length)
  ScalarType eccentricity;            ///< e [0, 1)
  ScalarType inclination;             ///< i (rad)
  ScalarType raan;                    ///< RAAN, Omega (rad)
  ScalarType argument_of_periapsis;   ///< omega (rad)
  ScalarType anomaly;                 ///< true/mean/eccentric anomaly according to AnomalyType (rad)
};

namespace detail {

template <typename ScalarType>
constexpr ScalarType Pi() {
  return static_cast<ScalarType>(3.141592653589793238462643383279502884L);
}

template <typename ScalarType>
constexpr ScalarType TwoPi() {
  return static_cast<ScalarType>(6.283185307179586476925286766559005768L);
}

template <typename ScalarType>
ScalarType NormalizeAngle(ScalarType angle) {
  const ScalarType two_pi = TwoPi<ScalarType>();
  angle = std::fmod(angle, two_pi);
  const ScalarType pi = Pi<ScalarType>();
  if (angle > pi) angle -= two_pi;
  if (angle < -pi) angle += two_pi;
  return angle;
}

// Solve Kepler's equation for elliptic orbits using Newton-Raphson on the eccentric anomaly.
template <typename ScalarType>
ScalarType SolveKeplerEquation(ScalarType mean_anomaly, ScalarType e, int max_iter = 50,
                               ScalarType tol = 1e-12) {
  if (e < 0 || e >= 1) {
    throw std::invalid_argument("Eccentricity must satisfy 0 <= e < 1 for Kepler solver.");
  }
  ScalarType M = NormalizeAngle(mean_anomaly);
  ScalarType E = std::abs(M) < Pi<ScalarType>() ? M : (M > 0 ? Pi<ScalarType>() : -Pi<ScalarType>());

  for (int i = 0; i < max_iter; ++i) {
    ScalarType f = E - e * std::sin(E) - M;
    ScalarType fp = static_cast<ScalarType>(1) - e * std::cos(E);
    ScalarType delta = f / fp;
    E -= delta;
    if (std::abs(delta) < tol) break;
  }
  return E;
}

template <typename ScalarType>
ScalarType TrueAnomalyFrom(ScalarType anomaly, ScalarType e, AnomalyType type) {
  if (type == AnomalyType::kTrue) {
    return anomaly;
  }
  const ScalarType pi = Pi<ScalarType>();
  if (e < 0 || e >= 1) {
    throw std::invalid_argument("Eccentricity must satisfy 0 <= e < 1.");
  }
  if (type == AnomalyType::kEccentric) {
    const ScalarType cos_E = std::cos(anomaly);
    const ScalarType sin_E = std::sin(anomaly);
    return std::atan2(std::sqrt(static_cast<ScalarType>(1) - e * e) * sin_E,
                      cos_E - e);
  }
  // type == kMean
  ScalarType E = SolveKeplerEquation(anomaly, e);
  return TrueAnomalyFrom(E, e, AnomalyType::kEccentric);
}

// Convert metric (m, m/s) or otherwise consistent units into [AU, AU/day].
template <typename ScalarType>
my_type::State<ScalarType> ToAUPerDay(const my_type::State<ScalarType>& state_metric,
                                      const param::AstroConstants<ScalarType>& astro_params) {
  if (astro_params.au == 0) {
    throw std::runtime_error("AstroConstants::au cannot be zero.");
  }
  constexpr ScalarType kSecondsPerDay = static_cast<ScalarType>(86400.0);
  const ScalarType pos_scale = static_cast<ScalarType>(1.0) / astro_params.au;
  const ScalarType vel_scale = kSecondsPerDay / astro_params.au;

  return {state_metric.x * pos_scale,
          state_metric.y * pos_scale,
          state_metric.z * pos_scale,
          state_metric.vx * vel_scale,
          state_metric.vy * vel_scale,
          state_metric.vz * vel_scale};
}

}  // namespace detail

/**
 * @brief Convert Keplerian orbital elements to a Cartesian state vector in the inertial frame.
 *
 * @param elements Barycentric elements (angles in rad, semi-major axis with same length unit as GM).
 * @param gm_central Gravitational parameter of the central body, consistent units with "a".
 * @param type Type of anomaly stored in elements.anomaly (true/mean/eccentric).
 * @return my_type::State<ScalarType> Position/velocity in the inertial frame.
 */
template <typename ScalarType>
my_type::State<ScalarType> ElementsToStateInertial(const KeplerianElements<ScalarType>& elements,
                                                   ScalarType gm_central,
                                                   AnomalyType type = AnomalyType::kTrue) {
  if (gm_central <= 0) {
    throw std::invalid_argument("gm_central must be positive.");
  }
  if (elements.semi_major_axis <= 0) {
    throw std::invalid_argument("Semi-major axis must be positive.");
  }
  if (elements.eccentricity < 0 || elements.eccentricity >= 1) {
    throw std::invalid_argument("Eccentricity must satisfy 0 <= e < 1 for elliptic orbits.");
  }

  const ScalarType nu = detail::TrueAnomalyFrom(elements.anomaly, elements.eccentricity, type);
  const ScalarType cos_nu = std::cos(nu);
  const ScalarType sin_nu = std::sin(nu);
  const ScalarType p = elements.semi_major_axis *
                       (static_cast<ScalarType>(1) - elements.eccentricity * elements.eccentricity);
  const ScalarType r_pf = p / (static_cast<ScalarType>(1) + elements.eccentricity * cos_nu);
  const ScalarType sqrt_mu_over_p = std::sqrt(gm_central / p);

  // Perifocal position/velocity
  const ScalarType x_pf = r_pf * cos_nu;
  const ScalarType y_pf = r_pf * sin_nu;
  const ScalarType vx_pf = -sqrt_mu_over_p * sin_nu;
  const ScalarType vy_pf = sqrt_mu_over_p * (elements.eccentricity + cos_nu);

  // Rotation matrix R = R3(raan) * R1(i) * R3(arg_periapsis)
  const ScalarType cO = std::cos(elements.raan);
  const ScalarType sO = std::sin(elements.raan);
  const ScalarType ci = std::cos(elements.inclination);
  const ScalarType si = std::sin(elements.inclination);
  const ScalarType cw = std::cos(elements.argument_of_periapsis);
  const ScalarType sw = std::sin(elements.argument_of_periapsis);

  const ScalarType r11 = cO * cw - sO * sw * ci;
  const ScalarType r12 = -cO * sw - sO * cw * ci;
  const ScalarType r13 = sO * si;
  const ScalarType r21 = sO * cw + cO * sw * ci;
  const ScalarType r22 = -sO * sw + cO * cw * ci;
  const ScalarType r23 = -cO * si;
  const ScalarType r31 = sw * si;
  const ScalarType r32 = cw * si;
  const ScalarType r33 = ci;

  const Vector3d<ScalarType> r_I{r11 * x_pf + r12 * y_pf, r21 * x_pf + r22 * y_pf,
                                 r31 * x_pf + r32 * y_pf};
  const Vector3d<ScalarType> v_I{r11 * vx_pf + r12 * vy_pf, r21 * vx_pf + r22 * vy_pf,
                                 r31 * vx_pf + r32 * vy_pf};

  return {r_I.x(), r_I.y(), r_I.z(), v_I.x(), v_I.y(), v_I.z()};
}

/**
 * @brief Convert SSB inertial states (meters, m/s) to CRTBP rotating ND coordinates.
 *
 * This is a thin wrapper around ConvertInertial2RotatingV2 that applies the necessary
 * unit conversion (m -> AU, m/s -> AU/day) before calling the existing converter.
 */
template <typename ScalarType>
my_type::State<ScalarType> ConvertSSBStateToRotating(
    const my_type::State<ScalarType>& target_state_metric,
    const my_type::State<ScalarType>& p2_state_metric,
    const param::AstroConstants<ScalarType>& astro_params) {
  const auto target_state_AU = detail::ToAUPerDay(target_state_metric, astro_params);
  const auto p2_state_AU = detail::ToAUPerDay(p2_state_metric, astro_params);
  return ConvertInertial2RotatingV2(target_state_AU, p2_state_AU, astro_params);
}

/**
 * @brief Convert SSB Keplerian elements directly to CRTBP rotating ND coordinates.
 *
 * A Keplerian two-body assumption is used to create Cartesian SSB states for both the
 * small body and the secondary (p2) around the supplied central mass. The output is in the
 * rotating frame non-dimensionalized by the primary-secondary distance and mean motion.
 *
 * @param target_elements Keplerian elements of the object of interest (SSB inertial frame).
 * @param p2_elements Keplerian elements of the secondary (e.g., Earth) in the same frame.
 * @param astro_params Needed for AU/GM values used inside ConvertInertial2RotatingV2.
 * @param gm_central Gravitational parameter of the central mass for the Keplerian model (m^3/s^2).
 * @param type The anomaly type stored in the element sets (true/mean/eccentric).
 */
template <typename ScalarType>
my_type::State<ScalarType> ConvertSSBElementsToRotating(
    const KeplerianElements<ScalarType>& target_elements,
    const KeplerianElements<ScalarType>& p2_elements,
    const param::AstroConstants<ScalarType>& astro_params, ScalarType gm_central,
    AnomalyType type = AnomalyType::kTrue) {
  const auto target_state_metric = ElementsToStateInertial(target_elements, gm_central, type);
  const auto p2_state_metric = ElementsToStateInertial(p2_elements, gm_central, type);
  return ConvertSSBStateToRotating(target_state_metric, p2_state_metric, astro_params);
}

}  // namespace crtbp

#endif  // BARYCENTRIC_CONVERSION_HPP_
