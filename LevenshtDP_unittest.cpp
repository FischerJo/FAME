
#include <limits>   // numeric limits (max)

#include "gtest/gtest.h"
#include "LevenshtDP.h"


// test if everything is initialized correctly with 2 letter string
TEST(init_test, twoLetterMismatches)
{

    struct Compi
    {
        uint8_t operator()(const char c1, const char c2)
        {
            return c1==c2 ? 0 : 1;
        }
    } comp;

    std::string s1 = "AG";

    std::string s2_c = "CTAG";
    const char* s2_correct = s2_c.data() + 3;
    std::string s2_ic = "CCTCG";
    const char* s2_incorrect = s2_ic.data() + 4;

    // test with no error
    LevenshtDP<uint16_t, 0> l0_corr(s1,s2_correct);
    LevenshtDP<uint16_t, 0> l0_incorr(s1,s2_incorrect);

    l0_corr.runDPFill<Compi>(comp);
    l0_incorr.runDPFill<Compi>(comp);

    ASSERT_EQ(0, l0_corr.getEditDist());
    ASSERT_EQ(1, l0_incorr.getEditDist());

    // test with one error
    LevenshtDP<uint16_t, 1> l1_corr(s1,s2_correct);
    LevenshtDP<uint16_t, 1> l1_incorr(s1,s2_incorrect);

    l1_corr.runDPFill<Compi>(comp);
    l1_incorr.runDPFill<Compi>(comp);

    ASSERT_EQ(0, l1_corr.getEditDist());
    ASSERT_EQ(1, l1_incorr.getEditDist());

    // test with two errors
    LevenshtDP<uint16_t, 2> l2_corr(s1,s2_correct);
    LevenshtDP<uint16_t, 2> l2_incorr(s1,s2_incorrect);

    l2_corr.runDPFill<Compi>(comp);
    l2_incorr.runDPFill<Compi>(comp);

    ASSERT_EQ(0, l2_corr.getEditDist());
    ASSERT_EQ(1, l2_incorr.getEditDist());

    std::vector<ERROR_T> errTypes_corr;
    l2_corr.backtrackDP<Compi>(comp, errTypes_corr);
    std::vector<ERROR_T> errTypes_incorr;
    l2_incorr.backtrackDP<Compi>(comp, errTypes_incorr);

    ASSERT_EQ(std::vector<ERROR_T>({MATCHING, MATCHING}), errTypes_corr);
    ASSERT_EQ(std::vector<ERROR_T>({INSERTION, MATCHING}), errTypes_incorr);

    // introduce another error
    s2_ic[4] = 'T';
    LevenshtDP<uint16_t, 2> l2_incorr2(s1,s2_incorrect);

    l2_incorr2.runDPFill<Compi>(comp);
    ASSERT_EQ(2, l2_incorr2.getEditDist());
}


// test simple sequence with 2 mismatches
TEST(matching_test, simpleSeq)
{
    struct Compi
    {
        uint8_t operator()(const char c1, const char c2)
        {
            return c1==c2 ? 0 : 1;
        }
    } comp;

    std::string s1 =   "ACTTTGACCG";
    std::string s2 = "AAACTGTGACTG";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data() + 11);

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MATCHING, MISMATCH, MATCHING, MATCHING, MATCHING, MATCHING, MISMATCH, MATCHING});
    std::vector<ERROR_T> backtrace;
    lev.backtrackDP<Compi>(comp, backtrace);

    ASSERT_EQ(trace, backtrace);

}
// test simple sequence for reverse query
TEST(matching_test, simpleSeqRev)
{
    struct Compi
    {
        uint8_t operator()(const char c1, const char c2)
        {
            switch (c1)
            {
                case 'A':
                    return c2 == 'T' ? 0 : 1;
                    break;

                case 'C':
                    return c2 == 'G' ? 0 : 1;
                    break;

                case 'G':
                    return c2 == 'C' ? 0 : 1;
                    break;

                case 'T':
                    return c2 == 'A' ? 0 : 1;
                    break;
            }
            return 1;
        }
    } comp;

    std::string s1 = "TCCAG";
    std::string s2 = "GACTGAA";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data() + 6);

    lev.runDPFillRev<Compi>(comp);

    ASSERT_EQ(1, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MATCHING, MISMATCH, MATCHING});
    std::vector<ERROR_T> backtrace;
    lev.backtrackDPRev<Compi>(comp, backtrace);

    ASSERT_EQ(trace, backtrace);

}

// test a more complex sequence with insertions
TEST(matching_test, complexSeqIns)
{
    struct Compi
    {
        uint8_t operator()(const char c1, const char c2)
        {
            return c1==c2 ? 0 : 1;
        }
    } comp;

    std::string s1 = "ACTAAGTGACTG";
    std::string s2 = "TTTTACTGTGACTG";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data() + 13);

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MATCHING, INSERTION, INSERTION, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING});
    std::vector<ERROR_T> backtrace;
    lev.backtrackDP<Compi>(comp, backtrace);

    ASSERT_EQ(trace, backtrace);
}

// test a more complex sequence with deletions
TEST(matching_test, complexSeqDel)
{
    struct Compi
    {
        uint8_t operator()(const char c1, const char c2)
        {
            return c1==c2 ? 0 : 1;
        }
    } comp;

    std::string s1 = "TGTGACTT";
    std::string s2 = "TGTGACTAAT";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data() + 9);

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, DELETION, DELETION, MATCHING});
    std::vector<ERROR_T> backtrace;
    lev.backtrackDP<Compi>(comp, backtrace);

    ASSERT_EQ(trace, backtrace);
}

// test sequence with different types of errors
TEST(matching_test, complexSeqAll)
{
    struct Compi
    {
        uint8_t operator()(const char c1, const char c2)
        {
            return c1==c2 ? 0 : 1;
        }
    } comp;

    std::string s1 = "ACGTGAACCTG";
    std::string s2 = "AAACGGTGAACTG";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data() + 12);

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MATCHING, DELETION, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, INSERTION, MATCHING, MATCHING});
    std::vector<ERROR_T> backtrace;
    lev.backtrackDP<Compi>(comp, backtrace);

    ASSERT_EQ(trace, backtrace);
}

// TODO test indels for reverse
//
//
//
//
