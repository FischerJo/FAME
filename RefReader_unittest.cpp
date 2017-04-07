
#include "gtest/gtest.h"
#include "RefReader_istr.h"

// FOR THESE TESTS
//
//
// SET  READLEN=30 in CONST.h
//      KMERLEN=20

// Test if simple line is read correctly
TEST(readLine_test, simple)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seqS = "AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAA";
    std::vector<char> seq (seqS.begin(), seqS.end());
    std::string readSeqS = "TTT";
    std::vector<char> readSeq (readSeqS.begin(), readSeqS.end());
    bool lastC = false;
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);
    ASSERT_EQ(1, cpgTab.size());
    ASSERT_EQ(0, cpgTab[0].chrom);
    // + 3 offset because readSeq is already initialized with 3 letters
    ASSERT_EQ(31 - MyConst::READLEN + 2 + 3, cpgTab[0].pos);
    ASSERT_EQ(0, cpgStartTab.size());

    ASSERT_EQ(false, lastC);
}

// Test if lastC flag is set correctly when calling with active lastC
TEST(readLine_test, lastCPos)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seqS = "CGC";
    std::string readSeqS = "";
    std::vector<char> seq (seqS.begin(), seqS.end());
    std::vector<char> readSeq (readSeqS.begin(), readSeqS.end());
    bool lastC = true;
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);

    seqS = "CG";
    seq = std::vector<char>(seqS.begin(), seqS.end());
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(false, lastC);

}

// Test if lastC flag is set correctly when calling with inactive lastC
TEST(readLine_test, lastCNeg)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seqS = "CG";
    std::string readSeqS = "";
    std::vector<char> seq (seqS.begin(), seqS.end());
    std::vector<char> readSeq (readSeqS.begin(), readSeqS.end());
    bool lastC = false;
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(false, lastC);

    seqS = "CGC";
    seq = std::vector<char>(seqS.begin(), seqS.end());
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);

}

// Test if CpG is constructed according to lastC flag
TEST(readLine_test, firstCharG)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seqS = "GC";
    std::string readSeqS = "";
    std::vector<char> seq (seqS.begin(), seqS.end());
    std::vector<char> readSeq (readSeqS.begin(), readSeqS.end());
    bool lastC = false;
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(0, cpgStartTab.size());

    seqS = "G";
    seq = std::vector<char>(seqS.begin(), seqS.end());
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(1, cpgStartTab[0].pos);
    ASSERT_EQ(3, readSeq.size());
}

// Test if small letter characters are read correctly
//
TEST(readLine_test, smallChars)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seqS = "GCCTctaAGg";
    std::string readSeqS = "";
    std::vector<char> seq (seqS.begin(), seqS.end());
    std::vector<char> readSeq (readSeqS.begin(), readSeqS.end());
    bool lastC = false;
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(0, cpgStartTab.size());

    ASSERT_EQ(10, readSeq.size());
    seqS = "GCCTCTAAGG";
    std::vector<char> testSeq (seqS.begin(), seqS.end());
    ASSERT_EQ(testSeq, readSeq);
}

// test if small/ big letter CpG over 2 lines is read correctly
TEST(readLine_test, smallChars2)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seqS = "GC";
    std::string readSeqS = "";
    std::vector<char> seq (seqS.begin(), seqS.end());
    std::vector<char> readSeq (readSeqS.begin(), readSeqS.end());
    bool lastC = false;
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(0, cpgStartTab.size());

    seqS = "GC";
    seq = std::vector<char>(seqS.begin(), seqS.end());
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(1, cpgStartTab[0].pos);
    ASSERT_EQ(4, readSeq.size());

    seqS = "G";
    seq = std::vector<char>(seqS.begin(), seqS.end());
    readLine(seqS, lastC, 1, cpgTab, cpgStartTab, readSeq);
    ASSERT_EQ(false, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(2, cpgStartTab.size());
    ASSERT_EQ(3, cpgStartTab[1].pos);
    ASSERT_EQ(5, readSeq.size());
    readSeqS = "GCGCG";
    std::vector<char> readSeqTest (readSeqS.begin(), readSeqS.end());
    ASSERT_EQ(readSeqTest, readSeq);
}


// Test if CpG is constructed and sequence is read
// SimpleFile1.fa contains
//
// >C
// AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAA
//
TEST(readReference_test, SimpleFile1)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("SimpleFile1.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpG is correctly constructed
    ASSERT_EQ(0, cpgStartTab.size());
    ASSERT_EQ(1, cpgTab.size());
    ASSERT_EQ(0, cpgTab[0].chrom);
    ASSERT_EQ(31 - MyConst::READLEN + 2, cpgTab[0].pos);

    // test if full sequence is read (and nothing else)
    std::string seqS = "AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAA";
    std::vector<char> seq (seqS.begin(), seqS.end());
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(seq, genSeq[0]);

}

// Test if CpGs at start and end are read correctly
// SimpleFile2.fa contains
//
// >C
// AACGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTCGTGCTG
//
TEST(readReference_test, SimpleFile2)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("SimpleFile2.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpGs are correctly constructed
    ASSERT_EQ(1, cpgTab.size());
    ASSERT_EQ(0, cpgTab[0].chrom);
    ASSERT_EQ(38 - MyConst::READLEN + 2, cpgTab[0].pos);
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(0, cpgStartTab[0].chrom);
    ASSERT_EQ(2, cpgStartTab[0].pos);

    // test if full sequence is read (and nothing else)
    std::string seqS = "AACGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTCGTGCTG";
    std::vector<char> seq (seqS.begin(), seqS.end());
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(seq, genSeq[0]);
}

// Test if multiline sequence is read correctly
// MultilineFile1.fa contains
//
// >C
// AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAA
// AACGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTCGTGCTG
//
TEST(readReference_test, MultilineFile1)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("MultilineFile1.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpG is correctly constructed
    ASSERT_EQ(0, cpgStartTab.size());
    ASSERT_EQ(3, cpgTab.size());
    ASSERT_EQ(0, cpgTab[0].chrom);
    ASSERT_EQ(0, cpgTab[1].chrom);
    ASSERT_EQ(0, cpgTab[2].chrom);
    ASSERT_EQ(31 - MyConst::READLEN + 2, cpgTab[0].pos);
    ASSERT_EQ(66 - MyConst::READLEN + 2, cpgTab[1].pos);
    ASSERT_EQ(102 - MyConst::READLEN + 2, cpgTab[2].pos);

    // test if full sequence is read (and nothing else)
    std::string seqS = "AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAAAACGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTCGTGCTG";
    std::vector<char> seq (seqS.begin(), seqS.end());
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(seq, genSeq[0]);
}


// Test if multiline sequence with CG over multiple lines is read correctly
// MultilineFile2.fa contains
//
// >C
// AACTCGTTTGGTAATTGTGAATTATGGGGAAGGTGAATTAAATC
// GTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTAGTGCTG
//
TEST(readReference_test, MultilineFile2)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("MultilineFile2.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpG is correctly constructed
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(0, cpgStartTab[0].chrom);
    ASSERT_EQ(4, cpgStartTab[0].pos);
    ASSERT_EQ(1, cpgTab.size());
    ASSERT_EQ(0, cpgTab[0].chrom);
    ASSERT_EQ(43 - MyConst::READLEN + 2, cpgTab[0].pos);

    std::string seqS = "AACTCGTTTGGTAATTGTGAATTATGGGGAAGGTGAATTAAATCGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTAGTGCTG";
    std::vector<char> seq (seqS.begin(), seqS.end());
    // test if full sequence is read (and nothing else)
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(seq, genSeq[0]);
}



// Test if non primary assembly lines are skipped
// MultiseqFile1.fa contains
//
// >B
// ACTCTTCTCTGGGACTTCG
// >C
// AACTCTCCGAAA
// >B
// AACTCTCCCCCCCG
// >D
//
TEST(readReference_test, MultiseqFile1)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("MultiseqFile1.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpG is correctly constructed
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(0, cpgStartTab[0].chrom);
    ASSERT_EQ(7, cpgStartTab[0].pos);
    ASSERT_EQ(0, cpgTab.size());

    std::string seqS = "AACTCTCCGAAA";
    std::vector<char> seq (seqS.begin(), seqS.end());

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(seq, genSeq[0]);

}


// Test if multiple simple chromosome sequences with small letters are read correctly
// MultiseqFile2.fa contains
//
// >C
// AACTCTCCGAaa
// >C
// AACtcTcCGAAA
//
TEST(readReference_test, MultiseqFile2)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("MultiseqFile2.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpG is correctly constructed
    ASSERT_EQ(2, cpgStartTab.size());
    ASSERT_EQ(0, cpgStartTab[0].chrom);
    ASSERT_EQ(7, cpgStartTab[0].pos);
    ASSERT_EQ(1, cpgStartTab[1].chrom);
    ASSERT_EQ(7, cpgStartTab[1].pos);
    ASSERT_EQ(0, cpgTab.size());

    std::string seqS = "AACTCTCCGAAA";
    std::vector<char> seq (seqS.begin(), seqS.end());

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(2, genSeq.size());
    ASSERT_EQ(seq, genSeq[0]);
    ASSERT_EQ(seq, genSeq[1]);

}


// Test if empty line is skipped correctly when having multiple sequences
// MultiseqFile3.fa contains
//
// >C
// AACTCTCCGAAA
//
// >C
// AACTCTCCGAAA
//
TEST(readReference_test, MultiseqFile3)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("MultiseqFile3.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpG is correctly constructed
    ASSERT_EQ(2, cpgStartTab.size());
    ASSERT_EQ(0, cpgStartTab[0].chrom);
    ASSERT_EQ(7, cpgStartTab[0].pos);
    ASSERT_EQ(1, cpgStartTab[1].chrom);
    ASSERT_EQ(7, cpgStartTab[1].pos);
    ASSERT_EQ(0, cpgTab.size());

    std::string seqS = "AACTCTCCGAAA";
    std::vector<char> seq (seqS.begin(), seqS.end());

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(2, genSeq.size());
    ASSERT_EQ(seq, genSeq[0]);
    ASSERT_EQ(seq, genSeq[1]);

}



// Test if one line of one chromosome ends with a C and the next line of next chromosome with a G
// that this is not a CpG
// MultiSeqFile4.fa contains
//
// >C
// AACTCTCCTATGTGTAAATCCTTAGTTATACGTTTggttgtatatctc
// >C
// GAACTCTCCGAAAATACTTAAGCTTTGCTATAGTCCTG
//
TEST(readReference_test, MultiseqFile4)
{

    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char> > genSeq;
    readReference("MultiseqFile4.fa", cpgTab, cpgStartTab, genSeq);

    // test if CpG is correctly constructed
    ASSERT_EQ(1, cpgTab.size());
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(0, cpgTab[0].chrom);
    ASSERT_EQ(30 - MyConst::READLEN + 2, cpgTab[0].pos);
    ASSERT_EQ(1, cpgStartTab[0].chrom);
    ASSERT_EQ(8, cpgStartTab[0].pos);

    std::string seq1S = "AACTCTCCTATGTGTAAATCCTTAGTTATACGTTTGGTTGTATATCTC";
    std::vector<char> seq1 (seq1S.begin(), seq1S.end());
    std::string seq2S = "GAACTCTCCGAAAATACTTAAGCTTTGCTATAGTCCTG";
    std::vector<char> seq2 (seq2S.begin(), seq2S.end());

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(2, genSeq.size());
    ASSERT_EQ(seq1, genSeq[0]);
    ASSERT_EQ(seq2, genSeq[1]);

}

