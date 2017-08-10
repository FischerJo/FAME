
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

    std::string s2_c = "AGCT";
    const char* s2_correct = s2_c.data();
    std::string s2_ic = "CGCCT";
    const char* s2_incorrect = s2_ic.data();

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

    std::vector<ERROR_T> errTypes_corr = l2_corr.backtrackDP<Compi>(comp);
    std::vector<ERROR_T> errTypes_incorr = l2_incorr.backtrackDP<Compi>(comp);

    ASSERT_EQ(std::vector<ERROR_T>({MATCHING, MATCHING, MISMATCH, MISMATCH}), errTypes_corr);
    ASSERT_EQ(std::vector<ERROR_T>({MATCHING, MISMATCH, MISMATCH, MISMATCH}), errTypes_incorr);

    // introduce another error
    s2_ic[1] = 'T';
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

    std::string s1 = "ACTTTGACCG";
    std::string s2 = "ACTGTGACTGAA";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MISMATCH, MATCHING, MATCHING, MATCHING, MATCHING, MISMATCH, MATCHING, MATCHING, MATCHING, MISMATCH, MISMATCH});

    ASSERT_EQ(trace, lev.backtrackDP<Compi>(comp));

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
    std::string s2 = "ACTGTGACTGAATT";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, INSERTION, INSERTION, MATCHING, MATCHING, MATCHING, MISMATCH, MISMATCH});

    ASSERT_EQ(trace, lev.backtrackDP<Compi>(comp));
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

    std::string s1 = "ACTGACTT";
    std::string s2 = "ACTGTGACTTAA";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, MATCHING, DELETION, DELETION, MATCHING, MATCHING});

    ASSERT_EQ(trace, lev.backtrackDP<Compi>(comp));
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

    std::string s1 = "ACTAGTGATTG";
    std::string s2 = "ACTGTGACTGAAT";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill<Compi>(comp);

    ASSERT_EQ(2, lev.getEditDist());

    std::vector<ERROR_T> trace ({MATCHING, MATCHING, MISMATCH, MATCHING, MATCHING, MATCHING, MATCHING, INSERTION, MATCHING, MATCHING, MATCHING, MISMATCH, MISMATCH});

    ASSERT_EQ(trace, lev.backtrackDP<Compi>(comp));
}
