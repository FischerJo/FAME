#ifndef BITFUNCTIONS_H
#define BITFUNCTIONS_H

#include "CONST.h"

namespace BitFun {

// METHOD FOR REVERSING bits OF 64bit
// following TAOCP by Knuth (Vol.4;1A Section 7.1.3 p.13 code from (70)

// constants for bitswap
constexpr uint64_t mu0 = 0x5555555555555555ULL;
constexpr uint64_t c4 =  0x0300c0303030c303ULL;
constexpr uint64_t c8 =  0x00c0300c03f0003fULL;
constexpr uint64_t c20 = 0x00000ffc00003fffULL;

// reverses the order of bits in x
inline uint64_t rev64(uint64_t x)
{

    x = ((x >> 1) & mu0) | ((x & mu0) << 1);
    uint64_t tmp = (x ^ (x >> 4)) & c4;
    x = x ^ tmp ^ (tmp << 4);
    tmp = (x ^ (x >> 8)) & c8;
    x = x ^ tmp ^ (tmp << 8);
    tmp = (x ^ (x >> 20)) & c20;
    x = x ^ tmp ^ (tmp << 20);
    return ((x >> 34) | (x << 30));

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
inline uint64_t revKmer(uint64_t kmer)
{
    return kmer ^ MyConst::KMERMASK;
}


} // end namespace BitFun

#endif /* BITFUNCTIONS_H */
