#ifndef VECTOR3D_HPP_
#define VECTOR3D_HPP_

#include <cmath>
#include <iostream>
#include <stdexcept>

class Vector3d {
   private:
    double x_;
    double y_;
    double z_;

   public:
    Vector3d(double x, double y, double z) : x_(x), y_(y), z_(z) {}

    Vector3d(const Vector3d& other) : x_(other.x_), y_(other.y_), z_(other.z_) {}

    // getter
    double x() const { return x_; }
    double y() const { return y_; }
    double z() const { return z_; }

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
        double x = x_;
        double y = y_;
        double z = z_;
        double norm = std::sqrt(x * x + y * y + z * z);
        return norm;
    }

    Vector3d normalise() const {
        Vector3d result_buf(x_, y_, z_);
        double norm = result_buf.magnitude();
        if (norm == 0) {
            throw std::runtime_error("normalized: zero magnitude");
        }
        Vector3d result(x_ / norm, y_ / norm, z_ / norm);
        return result;
    }

    Vector3d gaiseki(const Vector3d& other) const {
        return Vector3d(y_ * other.z_ - z_ * other.y_,
                        z_ * other.x_ - x_ * other.z_,
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

    Vector3d operator*(double scalar) const {
        return Vector3d(x_ * scalar, y_ * scalar, z_ * scalar);
    }

    Vector3d operator/(double scalar) const {
        if (scalar == 0) {
            throw std::runtime_error("operator/: division by zero");
        }
        return Vector3d(x_ / scalar, y_ / scalar, z_ / scalar);
    }
};

#endif  // VECTOR3D_HPP_
