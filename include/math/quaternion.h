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
        static constexpr float eps = 1e-6;
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

        using quat_t::operator=, quat_t::operator+, quat_t::operator-, quat_t::operator/;

        quaternion operator* (const float &val) const {
            return quaternion(data[0] * val, data[1] * val, data[2] * val, data[3] * val);
        }
        friend quaternion operator* (const float &val, const quaternion &quat) {
            return quat * val;
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

        quaternion operator/ (const quaternion &oth) {
            float modulus_square = oth.data[0] * oth.data[0] + oth.data[1] * oth.data[1] + oth.data[2] * oth.data[2] + oth.data[3] * oth.data[3];
            BSP_ASSERT(modulus_square > eps);
            modulus_square = 1.0f / modulus_square;
            float q0q0 = data[0] *  oth.data[0];
            float q0q1 = data[0] * -oth.data[1];
            float q0q2 = data[0] * -oth.data[2];
            float q0q3 = data[0] * -oth.data[3];
            float q1q0 = data[1] *  oth.data[0];
            float q1q1 = data[1] * -oth.data[1];
            float q1q2 = data[1] * -oth.data[2];
            float q1q3 = data[1] * -oth.data[3];
            float q2q0 = data[2] *  oth.data[0];
            float q2q1 = data[2] * -oth.data[1];
            float q2q2 = data[2] * -oth.data[2];
            float q2q3 = data[2] * -oth.data[3];
            float q3q0 = data[3] *  oth.data[0];
            float q3q1 = data[3] * -oth.data[1];
            float q3q2 = data[3] * -oth.data[2];
            float q3q3 = data[3] * -oth.data[3];
            return quaternion(
                (q0q0 - q1q1 - q2q2 - q3q3) * modulus_square,
                (q0q1 + q1q0 + q2q3 - q3q2) * modulus_square,
                (q0q2 - q1q3 + q2q0 + q3q1) * modulus_square,
                (q0q3 + q1q2 - q2q1 + q3q0) * modulus_square
            );
        }

        using quat_t::operator+=, quat_t::operator-=;
        quaternion &operator*= (const float &val) {
            data[0] *= val;
            data[1] *= val;
            data[2] *= val;
            data[3] *= val;
            return *this;
        }
        quaternion &operator*= (const quaternion &oth) {
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
            data[0] = q0q0 - q1q1 - q2q2 - q3q3;
            data[1] = q0q1 + q1q0 + q2q3 - q3q2;
            data[2] = q0q2 - q1q3 + q2q0 + q3q1;
            data[3] = q0q3 + q1q2 - q2q1 + q3q0;
            return *this;
        }
        quaternion &operator/= (const quaternion &oth) {
            float modulus_square = oth.data[0] * oth.data[0] + oth.data[1] * oth.data[1] + oth.data[2] * oth.data[2] + oth.data[3] * oth.data[3];
            BSP_ASSERT(modulus_square > eps);
            modulus_square = 1.0f / modulus_square;
            float q0q0 = data[0] *  oth.data[0];
            float q0q1 = data[0] * -oth.data[1];
            float q0q2 = data[0] * -oth.data[2];
            float q0q3 = data[0] * -oth.data[3];
            float q1q0 = data[1] *  oth.data[0];
            float q1q1 = data[1] * -oth.data[1];
            float q1q2 = data[1] * -oth.data[2];
            float q1q3 = data[1] * -oth.data[3];
            float q2q0 = data[2] *  oth.data[0];
            float q2q1 = data[2] * -oth.data[1];
            float q2q2 = data[2] * -oth.data[2];
            float q2q3 = data[2] * -oth.data[3];
            float q3q0 = data[3] *  oth.data[0];
            float q3q1 = data[3] * -oth.data[1];
            float q3q2 = data[3] * -oth.data[2];
            float q3q3 = data[3] * -oth.data[3];
            data[0] = (q0q0 - q1q1 - q2q2 - q3q3) * modulus_square;
            data[1] = (q0q1 + q1q0 + q2q3 - q3q2) * modulus_square;
            data[2] = (q0q2 - q1q3 + q2q0 + q3q1) * modulus_square;
            data[3] = (q0q3 + q1q2 - q2q1 + q3q0) * modulus_square;
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
            float modulus_square = data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3];
            if (modulus_square < eps) {
                return quaternion(1, 0, 0, 0);
            }
            modulus_square = 1.0f / modulus_square;
            return quaternion(
                data[0] * modulus_square,
                -data[1] * modulus_square,
                -data[2] * modulus_square,
                -data[3] * modulus_square
            );
        }

        // yaw pitch roll
        [[nodiscard]] matrix<3, 1> get_euler_angle() const {
            matrix<3, 1> ret;
            float modulus_square = data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3];
            if (modulus_square < eps) {
                return ret;
            }
            float modulus = 1.0f / std::sqrt(modulus_square);
            float q0 = data[0] * modulus, q1 = data[1] * modulus, q2 = data[2] * modulus, q3 = data[3] * modulus;

            ret[0] = std::atan2(2.0f * (q0 * q3 + q1 * q2), 1.0f - 2.0f * (q2 * q2 + q3 * q3));
            float sin_pitch = 2.0f * (q0 * q2 - q3 * q1);
            if (std::abs(sin_pitch) > 1.0f) {
                ret[1] = sin_pitch > 0 ? M_PI / 2.0f : -M_PI / 2.0f;
            } else {
                ret[1] = std::asin(sin_pitch);
            }
            ret[2] = std::atan2(2.0f * (q0 * q1 + q2 * q3), 1.0f - 2.0f * (q1 * q1 + q2 * q2));
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
            float modulus_square = data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3];
            if (modulus_square < eps) {
                return matrix <3, 3>::eye();
            }
            float modulus = 1.0f / std::sqrt(modulus_square);
            float q0 = data[0] * modulus, q1 = data[1] * modulus, q2 = data[2] * modulus, q3 = data[3] * modulus;

            return math::matrix <3, 3> {
                1.f - 2.f * (q2 * q2 + q3 * q3), 2.f * (q1 * q2 - q0 * q3), 2.f * (q1 * q3 + q0 * q2),
                2.f * (q1 * q2 + q0 * q3), 1.f - 2.f * (q1 * q1 + q3 * q3), 2.f * (q2 * q3 - q0 * q1),
                2.f * (q1 * q3 - q0 * q2), 2.f * (q2 * q3 + q0 * q1), 1.f - 2.f * (q1 * q1 + q2 * q2)
            };
        }
    };
}