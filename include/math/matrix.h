/**
 * @file matrix.h
 * @brief 矩阵类
 * @author tiny_fish
 * @version 0.0.1
 */

#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <type_traits>
#include <utility>

#include "bsp/def.h"

namespace math {
    template <int _row, int _col>
    class matrix {
        template <int, int>
        friend class matrix;

        static constexpr float eps = 1e-10f;

        [[nodiscard]] constexpr bool is_identity() const noexcept {
            if constexpr (_row != _col) return false;
            for (int r = 0; r < _row; r++) {
                for (int c = 0; c < _col; c++) {
                    if (data[r * _col + c] != (r == c ? 1.f : 0.f)) return false;
                }
            }
            return true;
        }
    public:
        static constexpr int rows() noexcept { return _row; }
        static constexpr int cols() noexcept { return _col; }

        explicit matrix(const float arr[_row * _col]) { load(arr); }
        explicit matrix(const std::array <float, _row * _col> &arr) { load(arr); }
        template <typename T = float>
        explicit matrix(const T &val = 0) {
            data.fill(0);
            if (val == 0) return;
            for (int i = 0; i < std::min(_row, _col); i++) {
                data[i * _col + i] = val;
            }
        }
        matrix(const matrix&) = default;
        matrix(const std::initializer_list <float> &val) {
            BSP_ASSERT(val.size() == _row * _col);
            auto p = val.begin();
            for (int i = 0; i < _row * _col; i++) data[i] = *(p++);
        }
        ~matrix() = default;

        static std::pair <int, int> size() noexcept {
            return std::make_pair(_row, _col);
        }

        // 取值
        float &operator() (const int &r, const int &c) {
            BSP_ASSERT(0 <= r and r < _row and 0 <= c and c < _col);
            return data[r * _col + c];
        }
        const float &operator() (const int &r, const int &c) const {
            BSP_ASSERT(0 <= r and r < _row and 0 <= c and c < _col);
            return data[r * _col + c];
        }

        // NOTE: 20260614 修改 [i] 为取 data 中的第 i 个元素, 上面括号取值不变
        float &operator[] (int idx) {
            BSP_ASSERT(0 <= idx and idx < _row * _col);
            return data[idx];
        }
        const float &operator[] (int idx) const {
            BSP_ASSERT(0 <= idx and idx < _row * _col);
            return data[idx];
        }

        // 基本运算
        matrix& operator=(const matrix&) = default;
        template <typename T>
        constexpr matrix operator+ (T val) const {
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = data[i] + val;
            }
            return ret;
        }
        constexpr matrix operator+ (const matrix &oth) const {
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = data[i] + oth.data[i];
            }
            return ret;
        }
        constexpr matrix operator- () const {
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = -data[i];
            }
            return ret;
        }
        constexpr matrix operator- (const float &val) const {
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = data[i] - val;
            }
            return ret;
        }
        constexpr matrix operator- (const matrix &oth) const {
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = data[i] - oth.data[i];
            }
            return ret;
        }
        constexpr matrix operator* (const float &val) const {
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = data[i] * val;
            }
            return ret;
        }
        constexpr friend matrix operator* (const float &val, const matrix &oth) {
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = val * oth.data[i];
            }
            return ret;
        }
        template <int el>
        constexpr matrix <_row, el> operator* (const matrix <_col, el> &oth) const {
            if constexpr (_row == _col) {
                if (is_identity()) return oth;
            }

            matrix <_row, el> ret(0);
            for (int i = 0; i < _row; i++) {
                const int dst = i * el;
                const int lhs = i * _col;
                for (int k = 0; k < _col; k++) {
                    const float val = data[lhs + k];
                    const int rhs = k * el;
                    for (int j = 0; j < el; j++) {
                        ret.data[dst + j] += val * oth.data[rhs + j];
                    }
                }
            }
            return ret;
        }
        constexpr matrix operator/ (const float &val) const {
            BSP_ASSERT(std::abs(val) > eps);
            const float reciprocal = 1.0f / val;
            matrix ret;
            for (int i = 0; i < _row * _col; i++) {
                ret.data[i] = data[i] * reciprocal;
            }
            return ret;
        }
        constexpr matrix &operator+= (float val) {
            for (int i = 0; i < _row * _col; i++) data[i] += val;
            return *this;
        }
        constexpr matrix &operator+= (const matrix &oth) {
            for (int i = 0; i < _row * _col; i++) {
                data[i] += oth.data[i];
            }
            return *this;
        }
        constexpr matrix &operator-= (const float &val) {
            for (int i = 0; i < _row * _col; i++) {
                data[i] -= val;
            }
            return *this;
        }
        constexpr matrix &operator-= (const matrix &oth) {
            for (int i = 0; i < _row * _col; i++) {
                data[i] -= oth.data[i];
            }
            return *this;
        }
        constexpr matrix &operator*= (const float &val) {
            for (int i = 0; i < _row * _col; i++) {
                data[i] *= val;
            }
            return *this;
        }
        constexpr matrix &operator/= (const float &val) {
            BSP_ASSERT(std::abs(val) > eps);
            const float reciprocal = 1.0f / val;
            for (int i = 0; i < _row * _col; i++) {
                data[i] *= reciprocal;
            }
            return *this;
        }
        constexpr bool operator== (const matrix &oth) const {
            for (int i = 0; i < _row * _col; i++) {
                if (std::abs(data[i] - oth.data[i]) > eps) return false;
            }
            return true;
        }
        constexpr bool operator!= (const matrix &oth) const {
            return !(*this == oth);
        }
        template <int _r, int _c>
        constexpr std::enable_if_t <_r == _row, matrix <_row, _col + _c>> operator| (const matrix <_r, _c> &oth) const {
            matrix <_row, _col + _c> ret;
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _col + _c; j++) {
                    if (j < _col)
                        ret(i, j) = (*this)(i, j);
                    else
                        ret(i, j) = oth(i, j - _col);
                }
            }
            return ret;
        }
        constexpr matrix <_col, _row> T() const {
            matrix <_col, _row> ret;
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _col; j++) {
                    ret(j, i) = (*this)(i, j);
                }
            }
            return ret;
        }
        [[nodiscard]] float trace() const {
            float ret = 0;
            for (int i = 0; i < std::min(_row, _col); i++) {
                ret += (*this)(i, i);
            }
            return ret;
        }
        template <int _r = _row, int _c = _col>
        std::enable_if_t <_r == _c and _r <= 3, matrix> inv() const {
            if constexpr (_r == 3) {
                const float a = data[0], b = data[1], c = data[2];
                const float d = data[3], e = data[4], f = data[5];
                const float g = data[6], h = data[7], i = data[8];

                // 余子式（cofactors，未转置）
                const float C00 =  (e * i - f * h);
                const float C01 = -(d * i - f * g);
                const float C02 =  (d * h - e * g);

                const float C10 = -(b * i - c * h);
                const float C11 =  (a * i - c * g);
                const float C12 = -(a * h - b * g);

                const float C20 =  (b * f - c * e);
                const float C21 = -(a * f - c * d);
                const float C22 =  (a * e - b * d);

                // 行列式（按第一行展开）
                float det = a * C00 + b * C01 + c * C02;

                if (std::abs(det) < eps) {
                    // 此时认为矩阵不可逆
                    BSP_ASSERT(false);
                    return matrix(1);
                }
                const float id = 1.0f / det;

                // 逆矩阵 = adj(A) / det = Cofactor^T / det
                matrix ret; // 3x3
                ret.data[0] = C00 * id;  ret.data[1] = C10 * id;  ret.data[2] = C20 * id;
                ret.data[3] = C01 * id;  ret.data[4] = C11 * id;  ret.data[5] = C21 * id;
                ret.data[6] = C02 * id;  ret.data[7] = C12 * id;  ret.data[8] = C22 * id;
                return ret;
            }
            if constexpr (_r == 2) {
                const float a = data[0], b = data[1], c = data[2], d = data[3];
                const float det = a * d - b * c;
                if (std::abs(det) < eps) {
                    BSP_ASSERT(false);
                    return matrix(1);
                }
                return matrix({d, -b, -c, a}) / det;
            }
            if (std::abs(data[0]) < eps) {
                BSP_ASSERT(false);
                return matrix({1});
            }
            return matrix({1 / data[0]});
        }

        template <int _r = _row, int _c = _col>
        std::enable_if_t <_r == _c and (_r > 3), matrix> inv() const {
            constexpr int n = _r;

            float A[n*n], I[n*n];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    A[i*n + j] = data[i * n + j];
                    I[i*n + j] = (i == j) ? 1.0f : 0.0f;
                }
            }

            auto swap_rows = [&](int r1, int r2) {
                if (r1 == r2) return;
                for (int j = 0; j < n; ++j) {
                    std::swap(A[r1*n + j], A[r2*n + j]);
                    std::swap(I[r1*n + j], I[r2*n + j]);
                }
            };

            for (int k = 0; k < n; ++k) {
                // 选主元
                int p = k;
                float maxv = std::abs(A[k*n + k]);
                for (int i = k + 1; i < n; ++i) {
                    float v = std::abs(A[i*n + k]);
                    if (v > maxv) { maxv = v; p = i; }
                }
                swap_rows(k, p);

                float diag = A[k*n + k];

                // 如果主元太小，当作矩阵奇异
                if (std::abs(diag) < eps) {
                    BSP_ASSERT(false);
                    return matrix(1);
                }

                float inv_diag = 1.0f / diag;
                for (int j = 0; j < n; ++j) {
                    A[k*n + j] *= inv_diag;
                    I[k*n + j] *= inv_diag;
                }

                for (int i = 0; i < n; ++i) {
                    if (i == k) continue;
                    float f = A[i*n + k];
                    for (int j = 0; j < n; ++j) {
                        A[i*n + j] -= f * A[k*n + j];
                        I[i*n + j] -= f * I[k*n + j];
                    }
                }
            }

            matrix ret;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    ret.data[i * n + j] = I[i*n + j];
            return ret;
        }

        void load(const float *ptr) {
            memcpy(data.data(), ptr, _row * _col * sizeof(float));
        }

        void load(const std::array <float, _row * _col> &arr) {
            data = arr;
        }

        void save(float *ptr) const {
            memcpy(ptr, data.data(), _row * _col * sizeof(float));
        }

        void save(std::array <float, _row * _col> &arr) const {
            arr = data;
        }

        static constexpr matrix eye() {
            return matrix(1);
        }
        static constexpr matrix zeros() {
            return matrix(0);
        }
        static constexpr matrix ones() {
            matrix ret;
            ret.data.fill(1);
            return ret;
        }

        // Rot matrix
        template <int _r = _row, int _c = _col, int _dim = _r>
        static constexpr std::enable_if_t <_r == _c, matrix> rot_z(float rad) {
            static_assert(_dim >= 3 && _dim <= 4, "dimension can not bigger than '4' or smaller than '3'");
            auto temp = matrix<_dim, _dim>::eye();
            temp(0, 0) = cosf(rad); temp(0, 1) = -sinf(rad);
            temp(1, 0) = sinf(rad); temp(1, 1) = cosf(rad);
            return temp;
        }

        template <int _r = _row, int _c = _col, int _dim = _r>
        static constexpr std::enable_if_t <_r == _c, matrix> rot_y(float rad) {
            static_assert(_dim >= 3 && _dim <= 4, "dimension can not bigger than '4' or smaller than '3'");
            auto temp = matrix<_dim, _dim>::eye();
            temp(0, 0) = cosf(rad); temp(0, 2) = sinf(rad);
            temp(2, 0) = -sinf(rad); temp(2, 2) = cosf(rad);
            return temp;
        }

        template <int _r = _row, int _c = _col, int _dim = _r>
        static constexpr std::enable_if_t <_r == _c, matrix> rot_x(float rad) {
            static_assert(_dim >= 3 && _dim <= 4, "dimension can not bigger than '4' or smaller than '3'");
            auto temp = matrix<_dim, _dim>::eye();
            temp(1, 1) = cosf(rad); temp(1, 2) = -sinf(rad);
            temp(2, 1) = sinf(rad); temp(2, 2) = cosf(rad);
            return temp;
        }

        template <int _r = _row, int _c = _col, int _dim = _r>
        static constexpr std::enable_if_t <_r == _c, matrix> rot_euler(float rol, float pit, float yaw) {
            static_assert(_dim >= 3 && _dim <= 4, "dimension can not bigger than '4' or smaller than '3'");
            auto matrix_rol = matrix::rot_x <_dim> (rol);
            auto matrix_pit = matrix::rot_y <_dim> (pit);
            auto matrix_yaw = matrix::rot_z <_dim> (yaw);
            auto matrix_euler = matrix_yaw * matrix_pit * matrix_rol;
            return matrix_euler;
        }

    protected:
        std::array <float, _row * _col> data;
    };
}
