/**
 * @file quaternion.h
 * @pre matrix.h
 * @brief 四元数类
 * @author tiny_fish
 * @version 0.0.1
 */

#pragma once
#include "math/matrix.h"

namespace math {
    class quaternion : public matrix <4, 1> {
        using quat_t = matrix <4, 1>;
        static constexpr float eps = 1e-12f;

        [[nodiscard]] float norm_square() const {
            return data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3];
        }
    public:
        explicit quaternion(const quat_t &mat) : quat_t(mat) { }
        explicit quaternion(float real = 1, float imag_x = 0, float imag_y = 0, float imag_z = 0) : quat_t({real, imag_x, imag_y, imag_z}) {}
        ~quaternion() = default;

        float &operator[] (const int &idx) {
            BSP_ASSERT(0 <= idx and idx < 4);
            return data[idx];
        }
        const float &operator[] (const int &idx) const {
            BSP_ASSERT(0 <= idx and idx < 4);
            return data[idx];
        }

        using quat_t::operator=;

        quaternion operator+ (const quaternion &oth) const {
            return quaternion(data[0] + oth.data[0], data[1] + oth.data[1], data[2] + oth.data[2], data[3] + oth.data[3]);
        }
        quaternion operator- () const {
            return quaternion(-data[0], -data[1], -data[2], -data[3]);
        }
        quaternion operator- (const quaternion &oth) const {
            return quaternion(data[0] - oth.data[0], data[1] - oth.data[1], data[2] - oth.data[2], data[3] - oth.data[3]);
        }

        quaternion operator* (const float &val) const {
            return quaternion(data[0] * val, data[1] * val, data[2] * val, data[3] * val);
        }
        friend quaternion operator* (const float &val, const quaternion &quat) {
            return quat * val;
        }

        quaternion operator/ (const float val) const {
            BSP_ASSERT(std::abs(val) > eps);
            return *this * (1.0f / val);
        }

        quaternion operator* (const quaternion &oth) const {
            const float q0q0 = data[0] * oth.data[0];
            const float q0q1 = data[0] * oth.data[1];
            const float q0q2 = data[0] * oth.data[2];
            const float q0q3 = data[0] * oth.data[3];
            const float q1q0 = data[1] * oth.data[0];
            const float q1q1 = data[1] * oth.data[1];
            const float q1q2 = data[1] * oth.data[2];
            const float q1q3 = data[1] * oth.data[3];
            const float q2q0 = data[2] * oth.data[0];
            const float q2q1 = data[2] * oth.data[1];
            const float q2q2 = data[2] * oth.data[2];
            const float q2q3 = data[2] * oth.data[3];
            const float q3q0 = data[3] * oth.data[0];
            const float q3q1 = data[3] * oth.data[1];
            const float q3q2 = data[3] * oth.data[2];
            const float q3q3 = data[3] * oth.data[3];
            return quaternion(
                q0q0 - q1q1 - q2q2 - q3q3,
                q0q1 + q1q0 + q2q3 - q3q2,
                q0q2 - q1q3 + q2q0 + q3q1,
                q0q3 + q1q2 - q2q1 + q3q0
            );
        }

        quaternion operator/ (const quaternion &oth) const {
            return *this * oth.inv();
        }

        quaternion &operator+= (const quaternion &oth) {
            data[0] += oth.data[0];
            data[1] += oth.data[1];
            data[2] += oth.data[2];
            data[3] += oth.data[3];
            return *this;
        }
        quaternion &operator-= (const quaternion &oth) {
            data[0] -= oth.data[0];
            data[1] -= oth.data[1];
            data[2] -= oth.data[2];
            data[3] -= oth.data[3];
            return *this;
        }
        quaternion &operator*= (const float &val) {
            data[0] *= val;
            data[1] *= val;
            data[2] *= val;
            data[3] *= val;
            return *this;
        }
        quaternion &operator*= (const quaternion &oth) {
            *this = *this * oth;
            return *this;
        }
        quaternion &operator/= (const quaternion &oth) {
            *this = *this / oth;
            return *this;
        }

        using quat_t::operator==, quat_t::operator!=;

        [[nodiscard]] float real() const {
            return data[0];
        }
        [[nodiscard]] matrix <3, 1> imag() const {
            return matrix <3, 1> ({data[1], data[2], data[3]});
        }
        [[nodiscard]] quaternion conjugate() const {
            return quaternion(data[0], -data[1], -data[2], -data[3]);
        };
        [[nodiscard]] quaternion inv() const {
            const float square = norm_square();
            BSP_ASSERT(square > eps);
            if (square <= eps) return quaternion();
            const float reciprocal = 1.0f / square;
            return quaternion(
                data[0] * reciprocal,
                -data[1] * reciprocal,
                -data[2] * reciprocal,
                -data[3] * reciprocal
            );
        }

        // yaw pitch roll
        [[nodiscard]] matrix<3, 1> get_euler_angle() const {
            matrix<3, 1> ret;
            const float square = norm_square();
            if (square <= eps) return ret;

            const float q0 = data[0], q1 = data[1], q2 = data[2], q3 = data[3];
            ret[0] = std::atan2(2.0f * (q0 * q3 + q1 * q2), square - 2.0f * (q2 * q2 + q3 * q3));
            const float sin_pitch = 2.0f * (q0 * q2 - q3 * q1) / square;
            ret[1] = std::asin(std::clamp(sin_pitch, -1.0f, 1.0f));
            ret[2] = std::atan2(2.0f * (q0 * q1 + q2 * q3), square - 2.0f * (q1 * q1 + q2 * q2));
            return ret;
        }

        [[nodiscard]] matrix <4, 4> get_self_matrix() const {
            return matrix <4, 4> {
                data[0], -data[1], -data[2], -data[3],
                data[1], data[0], -data[3], data[2],
                data[2], data[3], data[0], -data[1],
                data[3], -data[2], data[1], data[0]
            };
        }

        [[nodiscard]] matrix <3, 3> get_rotation_matrix() const {
            const float square = norm_square();
            if (square <= eps) return matrix <3, 3>::eye();

            const float scale = 2.0f / square;
            const float q0 = data[0], q1 = data[1], q2 = data[2], q3 = data[3];
            const float q1s = q1 * scale, q2s = q2 * scale, q3s = q3 * scale;

            return math::matrix <3, 3> {
                1.f - q2 * q2s - q3 * q3s, q1 * q2s - q0 * q3s, q1 * q3s + q0 * q2s,
                q1 * q2s + q0 * q3s, 1.f - q1 * q1s - q3 * q3s, q2 * q3s - q0 * q1s,
                q1 * q3s - q0 * q2s, q2 * q3s + q0 * q1s, 1.f - q1 * q1s - q2 * q2s
            };
        }
    };
}
