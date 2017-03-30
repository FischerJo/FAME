
#include "gtest/gtest.h"
#include "BitFunctions.h"

// ---------------------------------
// ----- TESTS FOR BIT REVERSE -----
// ---------------------------------

TEST(BitFunctions, reverseSimple1)
{
    uint64_t bitstr = 0xffffffffffffffffULL;
    ASSERT_EQ(bitstr, BitFun::rev64(bitstr));
}

TEST(BitFunctions, reverseSimple2)
{
    // sequence of 010101.....
    uint64_t bitstr = 0x5555555555555555ULL;
    // result should be the same since we keep pairs of bits intact
    ASSERT_EQ(bitstr, BitFun::rev64(bitstr));
}

TEST(BitFunctions, reverseSimple3)
{
    // sequence of 101010.....
    uint64_t bitstr = 0xaaaaaaaaaaaaaaaaULL;
    // result should be the same since we keep pairs of bits intact
    ASSERT_EQ(bitstr, BitFun::rev64(bitstr));
}

TEST(BitFunctions, reverseSimple4)
{
    uint64_t bitstr = 0x00000000ffffffffULL;
    uint64_t bitstrRes = 0xffffffff00000000ULL;
    ASSERT_EQ(bitstrRes, BitFun::rev64(bitstr));
}

TEST(BitFunctions, reverseHard1)
{
    uint64_t bitstr = 0xeaca0f9800fac0bbULL;
    uint64_t bitstrRes = 0xee03af0026f0a3abULL;
    ASSERT_EQ(bitstrRes, BitFun::rev64(bitstr));
}

TEST(BitFunctions, reverseHard2)
{
    uint64_t bitstr = 0xab582f34b1eaa7c8ULL;
    uint64_t bitstrRes = 0x23daab4e1cf825eaULL;
    ASSERT_EQ(bitstrRes, BitFun::rev64(bitstr));
}

// test idempotence of function
TEST(BitFunctions, reverseIdemp)
{
    uint64_t bitstr = 0xeaca0f9800fac0bbULL;
    ASSERT_EQ(bitstr, BitFun::rev64(BitFun::rev64(bitstr)));
    bitstr = 0x00000000ffffffffULL;
    ASSERT_EQ(bitstr, BitFun::rev64(BitFun::rev64(bitstr)));
}


// ---------------------------------
// ---------------------------------
// ---------------------------------


// ---------------------------------
// ------ TESTS FOR REV KMER -------
// ---------------------------------

TEST(BitFunctions, kmer1)
{
    // kmer: CTTAACTGATCCTGTACTAT
    // encoding 0111 1100 0001 1110 0011 0101 1110 1100 0111 0011
    uint64_t kmer = 0x0000007c1e35ec73ULL;


    // reverse complement kmer: ATAGTACAGGATCAGTTAAG
    // encoding: 0011 0010 1100 0100 1010 0011 0100 1011 1100 0010
    uint64_t revComp = 0x00000032c4a34bc2ULL;
    ASSERT_EQ(revComp, BitFun::revKmer(kmer));
}

TEST(BitFunctions, kmer2)
{
    // kmer: TAGTGCTAGTCTGTAGCGAT
    // encoding: 1100 1011 1001 1100 1011 0111 1011 0010 0110 0011
    uint64_t kmer = 0x000000cb9cb7b263ULL;

    // reverse complement kmer: ATCGCTACAGACTAGCACTA
    // encoding: 0011 0110 0111 0001 0010 0001 1100 1001 0001 1100
    uint64_t revComp = 0x000000367121c91cULL;
    ASSERT_EQ(revComp, BitFun::revKmer(kmer));
}

// ---------------------------------
// ---------------------------------
// ---------------------------------
