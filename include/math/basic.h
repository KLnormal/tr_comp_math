//
// Created by fish on 2025/10/7.
//

#pragma once

namespace math {
    inline void endian_reverse_16(void *addr) {
        const auto p = static_cast<uint16_t *>(addr);
        *p = __builtin_bswap16(*p);
    }
    inline void endian_reverse_32(void *addr) {
        const auto p = static_cast<uint32_t *>(addr);
        *p = __builtin_bswap32(*p);
    }
    template <typename T>
    T pow(T x, unsigned int y) {
        T ret(1);
        while (y) {
            if (y & 1) ret = ret * x;
            x = x * x;
            y >>= 1;
        }
        return ret;
    }
}