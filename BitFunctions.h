//	Metal - A fast methylation alignment and calling tool for WGBS data.
//	Copyright (C) 2017  Jonas Fischer
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	Jonas Fischer	jonaspost@web.de

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
namespace detail {
    constexpr uint64_t mu1 = 0x3333333333333333ULL;
    constexpr uint64_t mu2 = 0x0f0f0f0f0f0f0f0fULL;
    constexpr uint64_t mu3 = 0x00ff00ff00ff00ffULL;
    constexpr uint64_t mu4 = 0x0000ffff0000ffffULL;
}
inline uint64_t rev64(uint64_t x)
{
    x = ((x >> 2) & detail::mu1) | ((x & detail::mu1) << 2);
    x = ((x >> 4) & detail::mu2) | ((x & detail::mu2) << 4);
    x = ((x >> 8) & detail::mu3) | ((x & detail::mu3) << 8);
    x = ((x >> 16) & detail::mu4) | ((x & detail::mu4) << 16);
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
inline uint64_t revKmer(uint64_t kmer)
{
    // empty bits (i.e. bits that are untouched because kmer is not as big as 64 bits)
    constexpr unsigned int emptyBits = (64 - (2 * MyConst::KMERLEN));
    return (BitFun::rev64(kmer ^ MyConst::KMERMASK) >> emptyBits);
}


// compute a bitmask from a given sequence encoded with alphabet above
//
// NOTE:    if sequence is smaller then 64 bit, all remaining bits will
//          be 1 after code execution
inline uint64_t getMask(uint64_t seq)
{
    const uint64_t x = seq & 0x5555555555555555ULL;
    uint64_t y = seq & 0xaaaaaaaaaaaaaaaaULL;
    y = (y >> 1) & x;
    return ~((x^y) << 1);
}

// return the encoding of the DNA letter c according to encoding scheme
// reverse complement version of c with suffix 'Rev'
// A -> 00
// C -> 01
// G -> 10
// T -> 11
//
inline uint64_t getBitRepr(const char& c)
{

    switch (c)
    {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
    }
    std::cerr << "Cannot compute bit representation of letter \"" << c << "\"\n\n";
    return 0;
}
inline uint64_t getBitReprRev(const char& c)
{

    switch (c)
    {
        case 'A':
            return 3;
        case 'C':
            return 2;
        case 'G':
            return 1;
        case 'T':
            return 0;
    }
    std::cerr << "Cannot compute bit representation of letter \"" << c << "\"\n\n";
    return 0;
}

} // end namespace BitFun

#endif /* BITFUNCTIONS_H */
