
#include <limits>   // numeric limits (max)

#include "gtest/gtest.h"
#include "LevenshtDP.h"


// test if everything is initialized correctly with 2 letter string
TEST(init_test, twoLetterMismatches)
{
    std::string s1 = "AG";

    std::string s2_c = "AGCT";
    const char* s2_correct = s2_c.data();
    std::string s2_ic = "CGCCT";
    const char* s2_incorrect = s2_ic.data();

    // test with no error
    LevenshtDP<uint16_t, 0> l0_corr(s1,s2_correct);
    LevenshtDP<uint16_t, 0> l0_incorr(s1,s2_incorrect);

    l0_corr.runDPFill();
    l0_incorr.runDPFill();

    ASSERT_EQ(0, l0_corr.getEditDist());
    ASSERT_EQ(1, l0_incorr.getEditDist());

    // test with one error
    LevenshtDP<uint16_t, 1> l1_corr(s1,s2_correct);
    LevenshtDP<uint16_t, 1> l1_incorr(s1,s2_incorrect);

    l1_corr.runDPFill();
    l1_incorr.runDPFill();

    ASSERT_EQ(0, l1_corr.getEditDist());
    ASSERT_EQ(1, l1_incorr.getEditDist());

    // test with two errors
    LevenshtDP<uint16_t, 2> l2_corr(s1,s2_correct);
    LevenshtDP<uint16_t, 2> l2_incorr(s1,s2_incorrect);

    l2_corr.runDPFill();
    l2_incorr.runDPFill();

    ASSERT_EQ(0, l2_corr.getEditDist());
    ASSERT_EQ(1, l2_incorr.getEditDist());

    // introduce another error
    s2_ic[1] = 'T';
    LevenshtDP<uint16_t, 2> l2_incorr2(s1,s2_incorrect);

    l2_incorr2.runDPFill();
    ASSERT_EQ(2, l2_incorr2.getEditDist());
}


// test simple sequence with 2 mismatches
TEST(matching_test, simpleSeq)
{
    std::string s1 = "ACTTTGACCG";
    std::string s2 = "ACTGTGACTGAA";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill();

    ASSERT_EQ(2, lev.getEditDist());

}

// test a more complex sequence with insertions
TEST(matching_test, complexSeqIns)
{
    std::string s1 = "ACTAAGTGACTG";
    std::string s2 = "ACTGTGACTGAA";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill();

    ASSERT_EQ(2, lev.getEditDist());
}

// test a more complex sequence with deletions
TEST(matching_test, complexSeqDel)
{
    std::string s1 = "ACTGACTG";
    std::string s2 = "ACTGTGACTGAA";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill();

    ASSERT_EQ(2, lev.getEditDist());
}

// test sequence with different types of errors
TEST(matching_test, complexSeqAll)
{
    std::string s1 = "ACTAGTGATTG";
    std::string s2 = "ACTGTGACTGAA";

    LevenshtDP<uint16_t, 2> lev(s1, s2.data());

    lev.runDPFill();

    ASSERT_EQ(2, lev.getEditDist());
}
