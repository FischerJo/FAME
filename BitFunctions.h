#ifndef BITFUNCTIONS_H
#define BITFUNCTIONS_H

#include "CONST.h"

namespace BitFun {

// METHOD FOR REVERSING bits OF 64bit
// following TAOCP by Knuth (Vol.4;1A Section 7.1.3 p.13 code from (70)

// constants for bitswap
// constexpr uint64_t mu0 = 0x5555555555555555ULL;
// constexpr uint64_t c4 =  0x0300c0303030c303ULL;
// constexpr uint64_t c8 =  0x00c0300c03f0003fULL;
// constexpr uint64_t c20 = 0x00000ffc00003fffULL;

// reverses the order of bits in x
// inline uint64_t rev64(uint64_t x)
// {
//
//     x = ((x >> 1) & mu0) | ((x & mu0) << 1);
//     uint64_t tmp = (x ^ (x >> 4)) & c4;
//     x = x ^ tmp ^ (tmp << 4);
//     tmp = (x ^ (x >> 8)) & c8;
//     x = x ^ tmp ^ (tmp << 8);
//     tmp = (x ^ (x >> 20)) & c20;
//     x = x ^ tmp ^ (tmp << 20);
//     return ((x >> 34) | (x << 30));
//
// }

// METHOD FOR REVERSING bits OF 64bit KEEPING PAIRS OF BITS INTACT (i.e. our encoding)
// so x_63 x_62 x_61 ... x_3 x_2 x_1 x_0
// to x_1 x_0 x_3 x_2 ... x_63 x_62
// following TAOCP by Knuth (Vol.4;1A Section 7.1.3 p.13 code from (50)
// adapted code to keep pairs intact (i.e. skip first swap)
//

// constants for bitswap
constexpr uint64_t mu1 = 0x3333333333333333ULL;
constexpr uint64_t mu2 = 0x0f0f0f0f0f0f0f0fULL;
constexpr uint64_t mu3 = 0x00ff00ff00ff00ffULL;
constexpr uint64_t mu4 = 0x0000ffff0000ffffULL;

inline uint64_t rev64(uint64_t x)
{
    x = ((x >> 2) & mu1) | ((x & mu1) << 2);
    x = ((x >> 4) & mu2) | ((x & mu2) << 4);
    x = ((x >> 8) & mu3) | ((x & mu3) << 8);
    x = ((x >> 16) & mu4) | ((x & mu4) << 16);
    return (x >> 32) | (x << 32);

}

// compute the reverse complement of DNA bitstring with encoding
// A -> 00
// C -> 01
// G -> 10
// T -> 11
//
// Note that A XOR 11 == 11 == T
//           T XOR 11 == 00 == A
//           C XOR 11 == 10 == G
//           G XOR 11 == 01 == C
inline uint64_t revKmer(const uint64_t kmer)
{
    // empty bits (i.e. bits that are untouched because kmer is not as big as 64 bits)
    constexpr unsigned int emptyBits = (64 - (2 * MyConst::KMERLEN));
    return (BitFun::rev64(kmer ^ MyConst::KMERMASK) >> emptyBits);
}


} // end namespace BitFun

#endif /* BITFUNCTIONS_H */
