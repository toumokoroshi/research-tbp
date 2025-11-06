#ifndef VECTOR3D_HPP_
#define VECTOR3D_HPP_

#include <cmath>
#include <stdexcept>

template <typename T>
class Vector3d {
 private:
  T x_;
  T y_;
  T z_;

 public:
  Vector3d() : x_(0), y_(0), z_(0) {}

  Vector3d(T x, T y, T z) : x_(x), y_(y), z_(z) {}

  Vector3d(const Vector3d& other) : x_(other.x_), y_(other.y_), z_(other.z_) {}

  // getter
  T x() const { return x_; }
  T y() const { return y_; }
  T z() const { return z_; }

  Vector3d& operator=(const Vector3d& other) {
    x_ = other.x_;
    y_ = other.y_;
    z_ = other.z_;
    return *this;
  }

  Vector3d(Vector3d&& other) noexcept : x_(other.x_), y_(other.y_), z_(other.z_) {
    other.x_ = 0;
    other.y_ = 0;
    other.z_ = 0;
  }

  Vector3d& operator=(Vector3d&& other) noexcept {
    x_ = other.x_;
    y_ = other.y_;
    z_ = other.z_;
    other.x_ = 0;
    other.y_ = 0;
    other.z_ = 0;
    return *this;
  }

  double magnitude() const {
    T x = x_;
    T y = y_;
    T z = z_;
    T norm = std::sqrt(x * x + y * y + z * z);
    return norm;
  }

  Vector3d normalise() const {
    Vector3d result_buf(x_, y_, z_);
    T norm = result_buf.magnitude();
    if (norm == 0) {
      throw std::runtime_error("normalized: zero magnitude");
    }
    Vector3d result(x_ / norm, y_ / norm, z_ / norm);
    return result;
  }

  Vector3d gaiseki(const Vector3d& other) const {
    return Vector3d(y_ * other.z_ - z_ * other.y_, z_ * other.x_ - x_ * other.z_,
                    x_ * other.y_ - y_ * other.x_);
  }

  double naiseki(const Vector3d& other) const {
    return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
  }

  Vector3d operator+(const Vector3d& other) const {
    return Vector3d(x_ + other.x_, y_ + other.y_, z_ + other.z_);
  }

  Vector3d operator-(const Vector3d& other) const {
    return Vector3d(x_ - other.x_, y_ - other.y_, z_ - other.z_);
  }

  Vector3d operator*(T scalar) const { return Vector3d(x_ * scalar, y_ * scalar, z_ * scalar); }

  Vector3d operator/(T scalar) const {
    if (scalar == 0) {
      throw std::runtime_error("operator/: division by zero");
    }
    return Vector3d(x_ / scalar, y_ / scalar, z_ / scalar);
  }
};

#endif  // VECTOR3D_HPP_
