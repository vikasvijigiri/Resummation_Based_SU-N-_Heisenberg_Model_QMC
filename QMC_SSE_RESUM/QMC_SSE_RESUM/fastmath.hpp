#include <cstdint>
#include <string.h>

bool myIsnan(float v) {
    std::uint32_t i;
    memcpy(&i, &v, 4);
    return ((i&0x7f800000)==0x7f800000)&&(i&0x7fffff);
}

bool myIsnan(double v) {
    std::uint64_t i;
    memcpy(&i, &v, 8);
    return ((i&0x7ff0000000000000)==0x7ff0000000000000)&&(i&0xfffffffffffff);
}
