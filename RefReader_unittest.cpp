
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
    std::string seq = "AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAA";
    std::string readSeq = "TTT";
    bool lastC = false;
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);
    ASSERT_EQ(1, cpgTab.size());
    ASSERT_EQ(0, cpgTab[0].chrom);
    // + 3 offset because readSeq is already initialized with 3 letters
    ASSERT_EQ(31 - MyConst::READLEN + 2 + 3, cpgTab[0].pos);
    ASSERT_EQ(0, cpgStartTab.size());

    EXPECT_EQ("TTT" + seq, readSeq);

    ASSERT_EQ(false, lastC);
}

// Test if lastC flag is set correctly when calling with active lastC
TEST(readLine_test, lastCPos)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seq = "CGC";
    std::string readSeq = "";
    bool lastC = true;
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);

    seq = "CG";
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(false, lastC);

}

// Test if lastC flag is set correctly when calling with inactive lastC
TEST(readLine_test, lastCNeg)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seq = "CG";
    std::string readSeq = "";
    bool lastC = false;
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(false, lastC);

    seq = "CGC";
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);

}

// Test if CpG is constructed according to lastC flag
TEST(readLine_test, firstCharG)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seq = "GC";
    std::string readSeq = "";
    bool lastC = false;
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(0, cpgStartTab.size());

    seq = "G";
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(1, cpgStartTab[0].pos);
    ASSERT_EQ(3, readSeq.size());
}

// Test if small letter characters are read correctly
TEST(readLine_test, smallChars)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seq = "GCCTctaAGg";
    std::string readSeq = "";
    bool lastC = false;
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(0, cpgStartTab.size());

    ASSERT_EQ(10, readSeq.size());
    ASSERT_EQ(seq, readSeq);
}

// test if small/ big letter CpG over 2 lines is read correctly
TEST(readLine_test, smallChars2)
{
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::string seq = "Gc";
    std::string readSeq = "";
    bool lastC = false;
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(0, cpgStartTab.size());

    seq = "GC";
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);

    ASSERT_EQ(true, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(1, cpgStartTab.size());
    ASSERT_EQ(1, cpgStartTab[0].pos);
    ASSERT_EQ(4, readSeq.size());

    seq = "g";
    readLine(seq, lastC, 1, cpgTab, cpgStartTab, readSeq);
    ASSERT_EQ(false, lastC);
    ASSERT_EQ(0, cpgTab.size());
    ASSERT_EQ(2, cpgStartTab.size());
    ASSERT_EQ(3, cpgStartTab[1].pos);
    ASSERT_EQ(5, readSeq.size());
    ASSERT_EQ("GcGCg", readSeq);
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
    std::string seq = "AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAA";
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[0]);

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
    std::string seq = "AACGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTCGTGCTG";
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[0]);
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
    std::string seq = "AACTCTTTGGTAATTGTGAATTATGGGGGGGCGAATTAAAAAACTGGGAAATGTGAACAACAAAAACGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTCGTGCTG";
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[0]);
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

    std::string seq = "AACTCGTTTGGTAATTGTGAATTATGGGGAAGGTGAATTAAATCGTTTAGTTGTGTTCTTATTTAAATGTGGTGTCTTTAGTGCTG";
    // test if full sequence is read (and nothing else)
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[0]);
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

    std::string seq = "AACTCTCCGAAA";

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(1, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[0]);

}


// Test if multiple simple chromosome sequences are read
// MultiseqFile2.fa contains
//
// >C
// AACTCTCCGAAA
// >C
// AACTCTCCGAAA
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

    std::string seq = "AACTCTCCGAAA";

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(2, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[0]);
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[1]);

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

    std::string seq = "AACTCTCCGAAA";

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(2, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[0]);
    ASSERT_EQ(std::vector<char>(seq.begin(),seq.end()), genSeq[1]);

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

    std::string seq1 = "AACTCTCCTATGTGTAAATCCTTAGTTATACGTTTggttgtatatctc";
    std::string seq2 = "GAACTCTCCGAAAATACTTAAGCTTTGCTATAGTCCTG";

    // test if full sequence is read (and nothing else)
    ASSERT_EQ(2, genSeq.size());
    ASSERT_EQ(std::vector<char>(seq1.begin(),seq1.end()), genSeq[0]);
    ASSERT_EQ(std::vector<char>(seq2.begin(),seq2.end()), genSeq[1]);

}

