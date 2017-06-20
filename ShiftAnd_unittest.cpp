
#include "gtest/gtest.h"


#include "ShiftAnd.h"


// test fixture
class ShiftAnd_test : public testing::Test {

    protected:


        virtual void SetUp()
        {

            lmap['A'] = 0;
            lmap['C'] = 1;
            lmap['G'] = 2;
            lmap['T'] = 3;
        }



        std::array<uint8_t, 256> lmap;

};


// tests the reset function for resetting the set of active states
// which are represented by bitvectors internally
TEST_F(ShiftAnd_test, reset)
{

    std::string seq = "ACCATGTGACTGCATG";

    ShiftAnd<0> sa0(seq, lmap);

    sa0.active[0].B_0 = 15;
    sa0.active[0].B_1 = 1;

    sa0.reset();

    ASSERT_EQ(1, sa0.active[0].B_0);
    ASSERT_EQ(0, sa0.active[0].B_1);

    ShiftAnd<2> sa2(seq, lmap);

    sa2.active[0].B_0 = 15;
    sa2.active[0].B_1 = 1;
    sa2.active[1].B_0 = 1;
    sa2.active[1].B_1 = 1;
    sa2.active[2].B_0 = 0;
    sa2.active[2].B_1 = 1;

    sa2.reset();

    ASSERT_EQ(1, sa2.active[0].B_0);
    ASSERT_EQ(0, sa2.active[0].B_1);
    ASSERT_EQ(3, sa2.active[1].B_0);
    ASSERT_EQ(0, sa2.active[1].B_1);
    ASSERT_EQ(7, sa2.active[2].B_0);
    ASSERT_EQ(0, sa2.active[2].B_1);


}


// tests if the bitmasks are set correctly for the simple sequence
// ACACACCCC
TEST_F(ShiftAnd_test, simpleBitmasks)
{
    std::string seq = "ACACACCCC";

    ShiftAnd<1> sa1(seq, lmap);

    const uint64_t maskA_0 = sa1.masks[lmap['A']].B_0;
    const uint64_t maskA_1 = sa1.masks[lmap['A']].B_1;
    const uint64_t maskC_0 = sa1.masks[lmap['C']].B_0;
    const uint64_t maskC_1 = sa1.masks[lmap['C']].B_1;
    const uint64_t maskG_0 = sa1.masks[lmap['G']].B_0;
    const uint64_t maskG_1 = sa1.masks[lmap['G']].B_1;
    const uint64_t maskT_0 = sa1.masks[lmap['T']].B_0;
    const uint64_t maskT_1 = sa1.masks[lmap['T']].B_1;

    const uint64_t full = 0xffffffffffffffffULL;
    ASSERT_EQ(0xfffffffffffffc2bULL, maskA_0);
    ASSERT_EQ(0xffffffffffffffd5ULL, maskC_0);
    ASSERT_EQ(0xfffffffffffffc01ULL, maskG_0);
    ASSERT_EQ(0xfffffffffffffc01ULL, maskT_0);
    ASSERT_EQ(full, maskA_1);
    ASSERT_EQ(full, maskC_1);
    ASSERT_EQ(full, maskG_1);
    ASSERT_EQ(full, maskT_1);

    ASSERT_EQ(0x0000000000000200ULL, sa1.accepted.B_0);
    ASSERT_EQ(0, sa1.accepted.B_1);
}

// tests if the bitmasks are set correctly for the simple sequence
// AAAAA.......ACG
// where the last 3 letters are at position 63 to 65 in the sequence
// to test for correct splitting between the two uint64_t
TEST_F(ShiftAnd_test, simpleBitmasksOverflow)
{

    std::string seq = "AAAAAAAAAA";
    seq = seq + seq + seq + seq + seq + seq + "AAAACG";

    ShiftAnd<1> sa1(seq, lmap);

    const uint64_t maskA_0 = sa1.masks[lmap['A']].B_0;
    const uint64_t maskA_1 = sa1.masks[lmap['A']].B_1;
    const uint64_t maskC_0 = sa1.masks[lmap['C']].B_0;
    const uint64_t maskC_1 = sa1.masks[lmap['C']].B_1;
    const uint64_t maskG_0 = sa1.masks[lmap['G']].B_0;
    const uint64_t maskG_1 = sa1.masks[lmap['G']].B_1;
    const uint64_t maskT_0 = sa1.masks[lmap['T']].B_0;
    const uint64_t maskT_1 = sa1.masks[lmap['T']].B_1;

    const uint64_t full = 0xffffffffffffffffULL;
    ASSERT_EQ(full, maskA_0);
    ASSERT_EQ(1, maskC_0);
    ASSERT_EQ(1, maskG_0);
    ASSERT_EQ(1, maskT_0);
    ASSERT_EQ(0xfffffffffffffff9ULL, maskA_1);
    ASSERT_EQ(0xfffffffffffffffaULL, maskC_1);
    ASSERT_EQ(0xfffffffffffffffcULL, maskG_1);
    ASSERT_EQ(0xfffffffffffffff8ULL, maskT_1);

    ASSERT_EQ(0, sa1.accepted.B_0);
    ASSERT_EQ(0x0000000000000004ULL, sa1.accepted.B_1);
}

// tests if too long patterns are still correctly processed
TEST_F(ShiftAnd_test, bitmaskLongPattern)
{
    std::string seq = "GGGGGGGGGGGGGGGG";
    seq += seq;
    seq += seq;
    std::string seq2 = "CCCCCCCCCCCCCCCC";
    seq2 += seq2;
    seq2 += seq2;

    // after the following line, we reached 128 letters
    seq += seq2;

    seq += "AAA";

    ShiftAnd<2> sa2(seq, lmap);

    const uint64_t maskA_0 = sa2.masks[lmap['A']].B_0;
    const uint64_t maskA_1 = sa2.masks[lmap['A']].B_1;
    const uint64_t maskC_0 = sa2.masks[lmap['C']].B_0;
    const uint64_t maskC_1 = sa2.masks[lmap['C']].B_1;
    const uint64_t maskG_0 = sa2.masks[lmap['G']].B_0;
    const uint64_t maskG_1 = sa2.masks[lmap['G']].B_1;
    const uint64_t maskT_0 = sa2.masks[lmap['T']].B_0;
    const uint64_t maskT_1 = sa2.masks[lmap['T']].B_1;

    const uint64_t full = 0xffffffffffffffffULL;
    ASSERT_EQ(1, maskA_0);
    ASSERT_EQ(1, maskC_0);
    ASSERT_EQ(full, maskG_0);
    ASSERT_EQ(1, maskT_0);
    ASSERT_EQ(0, maskA_1);
    ASSERT_EQ(0xfffffffffffffffeULL, maskC_1);
    ASSERT_EQ(1, maskG_1);
    ASSERT_EQ(0, maskT_1);

    ASSERT_EQ(0, sa2.accepted.B_0);
    ASSERT_EQ(0x8000000000000000ULL, sa2.accepted.B_1);
}

// tests if the bitmasks are set correctly for the sequence
// ATTATTTCCC
// to look if WGBS specific behaviour is achieved
TEST_F(ShiftAnd_test, bitmaskBisulfite)
{

    std::string seq = "ATTATTTCCC";

    ShiftAnd<1> sa1(seq, lmap);

    const uint64_t maskA_0 = sa1.masks[lmap['A']].B_0;
    const uint64_t maskC_0 = sa1.masks[lmap['C']].B_0;
    const uint64_t maskG_0 = sa1.masks[lmap['G']].B_0;
    const uint64_t maskT_0 = sa1.masks[lmap['T']].B_0;

    ASSERT_EQ(0xfffffffffffff813ULL, maskA_0);
    ASSERT_EQ(0xffffffffffffffedULL, maskC_0);
    ASSERT_EQ(0xfffffffffffff801ULL, maskG_0);
    ASSERT_EQ(0xfffffffffffff8edULL, maskT_0);
}


// simple matching test with
// same sequence for text as for pattern
// and with one substitution
TEST_F(ShiftAnd_test, matching_same)
{

    // string of size 20
    std::string seq = "AAAAAAAAAAAAAAAAAAAA";
    seq = seq + seq + seq;
    seq += seq;
    seq += "TTTTTTT";
    std::vector<char> t (seq.begin(), seq.end());
    // don't allow errors
    ShiftAnd<0> sa0(seq, lmap);
    // allow single error
    ShiftAnd<1> sa1(seq, lmap);

    // query sequence
    std::vector<size_t> matchings0;
    std::vector<uint16_t> errors0;
    sa0.querySeq(t.begin(), t.end(), matchings0, errors0);
    std::vector<size_t> matchings1;
    std::vector<uint16_t> errors1;
    sa1.querySeq(t.begin(), t.end(), matchings1, errors1);

    // check the accepting masks
    ASSERT_EQ(0, sa0.accepted.B_0);
    ASSERT_EQ(0, sa1.accepted.B_0);
    ASSERT_EQ(0x8000000000000000ULL, sa0.accepted.B_1);
    ASSERT_EQ(0x8000000000000000ULL, sa1.accepted.B_1);

    // check size of matchings
    ASSERT_EQ(1, matchings0.size());
    ASSERT_EQ(2, matchings1.size());
    ASSERT_EQ(126, matchings0[0]);
    ASSERT_EQ(125, matchings1[0]);
    ASSERT_EQ(126, matchings1[1]);

    // check the mismatches
    ASSERT_EQ(1, errors0.size());
    ASSERT_EQ(2, errors1.size());
    ASSERT_EQ(0, errors0[0]);
    ASSERT_EQ(1, errors1[0]);
    ASSERT_EQ(0, errors1[1]);

    // Exchange a letter to produce mismatch
    t[5] = 'C';

    matchings0.clear();
    matchings1.clear();
    errors0.clear();
    errors1.clear();
    sa0.querySeq(t.begin(), t.end(), matchings0, errors0);
    sa1.querySeq(t.begin(), t.end(), matchings1, errors1);

    ASSERT_EQ(0, matchings0.size());
    ASSERT_EQ(1, matchings1.size());
    ASSERT_EQ(126, matchings1[0]);

    ASSERT_EQ(1, errors1[0]);

}


// simple matching test of small pattern in larger sequence
// p = AGGCGAGGC
// t = AGGCGAGGCGAAGCGAGGC
TEST_F(ShiftAnd_test, matching_smaller)
{

    // init pattern and text
    std::string p = "AGGCGAGGC";
    std::string tSeq = "AGGCGAGGCGAAGCGAGGC";
    std::vector<char> t (tSeq.begin(), tSeq.end());

    // allow single error
    ShiftAnd<1> sa1(p, lmap);

    // query the text to automata
    std::vector<size_t> matchings;
    std::vector<uint16_t> errors;
    sa1.querySeq(t.begin(), t.end(), matchings, errors);


    // we should have one full match at pos 8, one with deletion at pos 7,
    // one with insertion at pos 9,
    // one with substitution 13,
    // another one with deletion at 18
    ASSERT_EQ(5, matchings.size());
    ASSERT_EQ(5, errors.size());

    ASSERT_EQ(7, matchings[0]);
    ASSERT_EQ(1, errors[0]);
    ASSERT_EQ(8, matchings[1]);
    ASSERT_EQ(0, errors[1]);
    ASSERT_EQ(9, matchings[2]);
    ASSERT_EQ(1, errors[2]);
    ASSERT_EQ(13, matchings[3]);
    ASSERT_EQ(1, errors[3]);
    ASSERT_EQ(18, matchings[4]);
    ASSERT_EQ(1, errors[4]);

    // query only last part
    matchings.clear();
    errors.clear();
    sa1.querySeq(t.begin() + 11, t.end(), matchings, errors);

    ASSERT_EQ(7, matchings[0]);
    ASSERT_EQ(1, errors[0]);
}


// simple matching test for pattern with Ts against text with Cs
TEST_F(ShiftAnd_test, matching_bisulfite)
{

    std::string p = "AGGTTATTC";
    std::string tSeq = "AGGTCACCCAAT";
    std::vector<char> t (tSeq.begin(), tSeq.end());

    ShiftAnd<1> sa1(p, lmap);

    std::vector<size_t> matchings;
    std::vector<uint16_t> errors;
    sa1.querySeq(t.begin(), t.end(), matchings, errors);

    ASSERT_EQ(3, matchings.size());
    ASSERT_EQ(3, errors.size());

    ASSERT_EQ(7, matchings[0]);
    ASSERT_EQ(1, errors[0]);
    ASSERT_EQ(8, matchings[1]);
    ASSERT_EQ(0, errors[1]);
    ASSERT_EQ(9, matchings[2]);
    ASSERT_EQ(1, errors[2]);

    // introduce a substitution error
    t[4]= 'G';

    matchings.clear();
    errors.clear();
    sa1.querySeq(t.begin(), t.end(), matchings, errors);

    ASSERT_EQ(1, matchings.size());
    ASSERT_EQ(1, errors.size());

    ASSERT_EQ(8, matchings[0]);
    ASSERT_EQ(1, errors[0]);

}
