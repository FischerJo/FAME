
#include "gtest/gtest.h"
#include "DnaBitStr.h"


// Test setting and reading a simple 32bp sequence of As
TEST(DnaBitStr_test, setSimple1)
{
    std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    uint64_t bitEnc = 0x0000000000000000ULL;
    uint64_t bitMask = 0xffffffffffffffffULL;
    uint64_t revBitEnc = 0xffffffffffffffffULL;

    DnaBitStr bitstr(32);
    bitstr.setBitStrN(seq, 0);

    uint64_t bits, bitsRev, mask, maskRev;
    // align kmer to lower bits
    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    for (unsigned int i = 0; i < (32 - MyConst::KMERLEN + 1); ++i)
    {

        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*i) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, maskRev);
    }


}

// Test setting and reading a 32bp sequence
TEST(DnaBitStr_test, setComplex1)
{
    // encoding: 1001 1111 0010 0111 1000 0001 1011 0010 0111 1000 1101 1011 1000 1110 0110 1110
    std::string seq = "GCTTAGCTGAACGTAGCTGATCGTGATGCGTG";
    uint64_t bitEnc = 0x9f2781b278db8d6dULL;
    uint64_t bitMask = 0xdff7fdff7fdfff7fULL;
    // reverse complement: CACGCATCACGATCAGCTACGTTCAGCTAAGC
    // encoding: 0100 0110 0100 1101 0001 1000 1101 0010 0111 0001 1011 1101 0010 0111 0000 1001
    uint64_t revBitEnc = 0x464d18d271bd2709ULL;
    uint64_t revBitMask = 0x777ddfdf7dfdf7fdULL;

    DnaBitStr bitstr(32);
    bitstr.setBitStrN(seq, 0);

    uint64_t bits, bitsRev, mask, maskRev;
    // align kmer to lower bits
    constexpr unsigned int shiftRight = (64 - (2 * MyConst::KMERLEN) );
    for (unsigned int i = 0; i < (32 - MyConst::KMERLEN + 1); ++i)
    {

        bits = bitstr.getSeqKmer(i);
        ASSERT_EQ( (bitEnc << 2*i) >> shiftRight, bits);

        bitsRev = bitstr.getSeqKmerRev(i);
        ASSERT_EQ( (revBitEnc << 2*i) >> shiftRight, bitsRev);

        mask = bitstr.getMaskKmer(i);
        ASSERT_EQ( (bitMask << 2*i) >> shiftRight, mask);

        maskRev = bitstr.getMaskKmerRev(i);
        ASSERT_EQ( (revBitMask << 2*i) >> shiftRight, maskRev);
    }


}
