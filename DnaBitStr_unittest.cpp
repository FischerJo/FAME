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

#include "gtest/gtest.h"
#include "DnaBitStr.h"


// Test setting and reading a simple 32bp sequence of As
TEST(DnaBitStr_test, setSimple1)
{
    std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const uint64_t bitEnc = 0x0000000000000000ULL;
    const uint64_t bitMask = 0xffffffffffffffffULL;
    const uint64_t revBitEnc = 0xffffffffffffffffULL;
    const uint64_t revBitMask = bitMask;

    DnaBitStr bitstr(32);
    bitstr.setBitStrN(std::move(seq), 0);

    uint64_t bits, bitsRev, mask, maskRev;
    // align kmer to lower bits
    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    for (unsigned int i = 0; i < (32 - MyConst::KMERLEN + 1); ++i)
    {

        bits = bitstr.getSeqKmer(i);
        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (revBitMask << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, maskRev);
    }
}

// Test setting and reading a simple 32bp sequence of Cs
TEST(DnaBitStr_test, setSimple2)
{
    std::string seq = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    const uint64_t bitEnc = 0x5555555555555555ULL;
    const uint64_t bitMask = bitEnc;
    const uint64_t revBitEnc = 0xaaaaaaaaaaaaaaaaULL;
    const uint64_t revBitMask = 0xffffffffffffffffULL;

    DnaBitStr bitstr(32);
    bitstr.setBitStrN(std::move(seq), 0);

    uint64_t bits, bitsRev, mask, maskRev;
    // align kmer to lower bits
    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    for (unsigned int i = 0; i < (32 - MyConst::KMERLEN + 1); ++i)
    {
        bits = bitstr.getSeqKmer(i);
        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (revBitMask << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, maskRev);
    }
}

// Test setting and reading a simple 32bp sequence of Gs
TEST(DnaBitStr_test, setSimple3)
{
    std::string seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    const uint64_t bitEnc = 0xaaaaaaaaaaaaaaaaULL;
    const uint64_t bitMask = 0xffffffffffffffffULL;
    const uint64_t revBitEnc = 0x5555555555555555ULL;
    const uint64_t revBitMask = revBitEnc;

    DnaBitStr bitstr(32);
    bitstr.setBitStrN(std::move(seq), 0);

    uint64_t bits, bitsRev, mask, maskRev;
    // align kmer to lower bits
    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    for (unsigned int i = 0; i < (32 - MyConst::KMERLEN + 1); ++i)
    {

        bits = bitstr.getSeqKmer(i);
        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (revBitMask << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, maskRev);
    }
}

// Test setting and reading a simple 32bp sequence of Ts
TEST(DnaBitStr_test, setSimple4)
{
    std::string seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const uint64_t bitMask = 0xffffffffffffffffULL;
    const uint64_t bitEnc = bitMask;
    const uint64_t revBitEnc = 0x0000000000000000ULL;
    const uint64_t revBitMask = bitMask;

    DnaBitStr bitstr(32);
    bitstr.setBitStrN(std::move(seq), 0);

    uint64_t bits, bitsRev, mask, maskRev;
    // align kmer to lower bits
    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    for (unsigned int i = 0; i < (32 - MyConst::KMERLEN + 1); ++i)
    {

        bits = bitstr.getSeqKmer(i);
        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (revBitMask << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, maskRev);
    }
}

// Test setting and reading a 32bp sequence
TEST(DnaBitStr_test, setComplex1)
{
    // encoding: 1001 1111 0010 0111 1000 0001 1011 0010 0111 1000 1101 1011 1000 1110 0110 1110
    std::string seq = "GCTTAGCTGAACGTAGCTGATCGTGATGCGTG";
    const uint64_t bitEnc = 0x9f2781b278db8e6eULL;
    const uint64_t bitMask = 0xdff7fdff7fdfff7fULL;
    // reverse complement: CACGCATCACGATCAGCTACGTTCAGCTAAGC
    // encoding: 0100 0110 0100 1101 0001 1000 1101 0010 0111 0001 1011 1101 0010 0111 0000 1001
    const uint64_t revBitEnc = 0x464d18d271bd2709ULL;
    const uint64_t revBitMask = 0x777ddfdf7dfdf7fdULL;

    DnaBitStr bitstr(32);
    bitstr.setBitStrN(std::move(seq), 0);

    uint64_t bits, bitsRev, mask, maskRev;
    // align kmer to lower bits
    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    for (unsigned int i = 0; i < (32 - MyConst::KMERLEN + 1); ++i)
    {

        bits = bitstr.getSeqKmer(i);
        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (revBitMask << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, maskRev);
    }


}

// Test setting and reading a simple 64bp sequence
TEST(DnaBitStr_test, setLong1)
{

    std::string seq0 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    std::string seq1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const uint64_t bitEnc = 0x0000000000000000ULL;
    const uint64_t bitMask = 0x000000ffffffffffULL;
    const uint64_t revBitEnc = bitMask;

    DnaBitStr bitstr(64);
    bitstr.setBitStrN(std::move(seq0), 0);
    bitstr.setBitStrN(std::move(seq1), 1);

    for (unsigned int i = 15; i <= 33; ++i)
    {

        uint64_t bits = bitstr.getSeqKmer(i);
        ASSERT_EQ(bitEnc, bits);

        uint64_t bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ(revBitEnc, bitsRev);

        uint64_t mask = bitstr.getMaskKmer(i);
        ASSERT_EQ(bitMask, mask);

        uint64_t maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ(bitMask, maskRev);
    }

}

// Test setting and reading a simple 64bp sequence
TEST(DnaBitStr_test, setLong2)
{

    std::string seq0 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    std::string seq1 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    const uint64_t bitEnc = 0x0000005555555555ULL;
    const uint64_t bitMask = bitEnc;
    const uint64_t revBitEnc = 0x000000aaaaaaaaaaULL;
    const uint64_t revBitMask = 0x000000ffffffffffULL;

    DnaBitStr bitstr(64);
    bitstr.setBitStrN(std::move(seq0), 0);
    bitstr.setBitStrN(std::move(seq1), 1);

    for (unsigned int i = 15; i <= 33; ++i)
    {

        uint64_t bits = bitstr.getSeqKmer(i);
        ASSERT_EQ(bitEnc, bits);

        uint64_t bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ(revBitEnc, bitsRev);

        uint64_t mask = bitstr.getMaskKmer(i);
        ASSERT_EQ(bitMask, mask);

        uint64_t maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ(revBitMask, maskRev);
    }

}

//Test setting and reading last bit segment
TEST(DnaBitStr_test, setLastSimple)
{
    // sequence of 21 bp
    // encoding: 0111 0111 0011 0110 1100 1011 0110 1100 1000 1101 1000
    std::string seq = "CTCTATCGTAGTCGTAGATCG";
    const uint64_t bitEnc = 0x7736cb6c8d800000ULL;
    // mask: 0111 0111 1111 0111 1111 1111 0111 1111 1111 1101 1100
    const uint64_t bitMask = 0x77f7ff7ffdc00000ULL;
    // sequence: CGATCTACGACTACGATAGAG
    // encoding: 0001 1000 1101 1100 0110 0001 1100 0110 0011 0010 0010
    const uint64_t revBitEnc = 0x0000018dc61c6322ULL;
    // mask: 0001 1111 1101 1111 0111 1101 1111 0111 1111 1111 1111
    const uint64_t revBitMask = 0x000001fdf7df7fffULL;

    DnaBitStr bitstr(21);
    bitstr.setBitStrLast(std::move(seq));

    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    uint64_t bits, bitsRev, mask, maskRev;
    for (unsigned int i = 0; i < 2; ++i)
    {

        bits = bitstr.getSeqKmer(i);
        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (revBitMask << 2*(32 - MyConst::KMERLEN - i)) >> shiftRight, maskRev);
    }

}

// Test setting and reading last bit segment with overlap
TEST(DnaBitStr_test, setLastOverlap)
{

    // pos 19                              |
    // pos 24                                   |
    std::string seq0 = "TGACCGTTCACCAATTATAGCGCTAAATGCTA";
    std::string seq1 = "AGTATGCAGCCC";

    DnaBitStr bitstr(44);
    bitstr.setBitStrN(std::move(seq0), 0);
    bitstr.setBitStrLast(std::move(seq1));

    const uint64_t bitEnc19 = 0x00000099c0e70b39ULL;
    const uint64_t bitMask19 = 0x000000ddfff7fffdULL;
    // GCATACTTAGCATTTAGCGC
    const uint64_t revBitEnc19 = 0x000000931f24fc99ULL;
    const uint64_t revBitMask19 = 0x000000dfdff7ffddULL;

    const uint64_t bitEnc24 = 0x000000039c2ce495ULL;
    const uint64_t bitMask24 = 0x000000ffdffff7d5ULL;
    // GGGCTGCATACTTAGCATTT
    const uint64_t revBitEnc24 = 0x000000a9e4c7c93fULL;
    const uint64_t revBitMask24 = 0x000000fdf7f7fdffULL;

    ASSERT_EQ(bitEnc19, bitstr.getSeqKmer(19));
    ASSERT_EQ(bitMask19, bitstr.getMaskKmer(19));
    ASSERT_EQ(revBitEnc19, bitstr.getSeqKmerRev(19));
    ASSERT_EQ(revBitMask19, bitstr.getMaskKmerRev(19));

    ASSERT_EQ(bitEnc24, bitstr.getSeqKmer(24));
    ASSERT_EQ(bitMask24, bitstr.getMaskKmer(24));
    ASSERT_EQ(revBitEnc24, bitstr.getSeqKmerRev(24));
    ASSERT_EQ(revBitMask24, bitstr.getMaskKmerRev(24));
}
