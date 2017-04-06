
#include "gtest/gtest.h"


#include "RefGenome.h"


// TO RUN THE FOLLOWING TESTS
// MAKE REFGENOME FUNCTIONS PUBLIC
//
// SET READLEN TO 30
//     KMERLEN TO 20




// Test hashing of a simple sequence with one normal CpGs
TEST(RefGenome_test, simple1)
{

    // set up sequence container
    std::string seq = "ATGTTGCCTAATTTCACTATTCAGGGTTATACGCCTGGAATATTCTAGGATTCCTAGTCAATTTAT";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTTAATTTTATTATTTAGGGTTATATGTTTGGAATATTTTAGGATTTTTAGTTAATTTAT";
    // reverse sequence
    std::string revSeq = "ATAAATTGACTAGGAATCCTAGAATATTCCAGGCGTATAACCCTGAATAGTGAAATTAGGCAACAT";
    std::string redRevSeq = "ATAAATTGATTAGGAATTTTAGAATATTTTAGGTGTATAATTTTGAATAGTGAAATTAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (cpgTab, cpgStart, genSeq);

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.kmerTable.size(), 0);

    for (unsigned int i = 5; i <= (63 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = (2*MyConst::READLEN - 2) - i + 5 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        // lookup if kmer is present in hash table
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer kRev = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, kRev.cpg);
        ASSERT_EQ(off, kRev.offset);
    }

    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer k = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, k.cpg);
        ASSERT_EQ(off | STRANDMASK, k.offset);
    }

    // test if nothing else was hashed
    for (unsigned int i = 0; i < ref.kmerTable.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.kmerTable[i].len);

    }
}



// Test with N before and after CpG
TEST(RefGenome_test, simpleWithN)
{

    // set up sequence container
    std::string seq = "ATGTTGCCTNATTTCACTATTCAGGGTTATACGCCTGGAATATTCTAGGANTCCTAGTCAATTTAT";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTTNATTTTATTATTTAGGGTTATATGTTTGGAATATTTTAGGANTTTTAGTTAATTTAT";
    // reverse sequence
    std::string revSeq = "ATAAATTGACTAGGANTCCTAGAATATTCCAGGCGTATAACCCTGAATAGTGAAATNAGGCAACAT";
    std::string redRevSeq = "ATAAATTGATTAGGANTTTTAGAATATTTTAGGTGTATAATTTTGAATAGTGAAATNAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (cpgTab, cpgStart, genSeq);

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.kmerTable.size(), 0);

    for (unsigned int i = 16; i <= (56 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = (2*MyConst::READLEN - 2) - i + 5 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        // lookup if kmer is present in hash table
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer kRev = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, kRev.cpg);
        ASSERT_EQ(off, kRev.offset);
    }

    for (unsigned int i = 10; i <= (50 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer k = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, k.cpg);
        ASSERT_EQ(off | STRANDMASK, k.offset);
    }

    // test if nothing else was hashed
    for (unsigned int i = 0; i < ref.kmerTable.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.kmerTable[i].len);

    }
}


// Test with N before and after CpG such that there is no kmer
TEST(RefGenome_test, simpleTooShort)
{

    // set up sequence container
    std::string seq = "ATGTTGCCTATTTCACTATTCAGGGTTNTACGCCTGGANTATTCTAGGATCCTAGTCAATTTAT";
    // reverse sequence
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (cpgTab, cpgStart, genSeq);


    // test if nothing else was hashed
    for (unsigned int i = 0; i < ref.kmerTable.size(); ++i)
    {

        ASSERT_EQ(0, ref.kmerTable[i].len);

    }
}

// Test hashing of CpG at start of sequence
TEST(RefGenome_test, simpleAtStart)
{

    // set up sequence container
    //                                                   |
    std::string seq = "ATGTTCGCCTATTTCACTATTCAGGGTTATAGCCTGGA";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTTGTTTATTTTATTATTTAGGGTTATAGTTTGGA";
    // reverse sequence
    std::string revSeq = "TCCAGGCTATAACCCTGAATAGTGAAATAGGCGAACAT";
    std::string redRevSeq = "TTTAGGTTATAATTTTGAATAGTGAAATAGGTGAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;

    std::vector<struct CpG> cpgStart;
    cpgStart.emplace_back(0, 5);

    RefGenome ref (cpgTab, cpgStart, genSeq);

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.kmerTable.size(), 0);

    for (unsigned int i = 3; i <= (38 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = (MyConst::READLEN + 5) - i + 3 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        // lookup if kmer is present in hash table
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer kRev = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, kRev.cpg ^ MyConst::INDMASK);
        ASSERT_EQ(MyConst::INDMASK, kRev.cpg & MyConst::INDMASK);
        ASSERT_EQ(off, kRev.offset);
    }

    for (unsigned int i = 0; i <= (35 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = i;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer k = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, k.cpg ^ MyConst::INDMASK);
        ASSERT_EQ(MyConst::INDMASK, k.cpg & MyConst::INDMASK);
        ASSERT_EQ(off | STRANDMASK, k.offset);
    }

    // test if nothing else was hashed
    for (unsigned int i = 0; i < ref.kmerTable.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.kmerTable[i].len);

    }
}


// Test hashing of CpG at end of sequence
TEST(RefGenome_test, simpleAtEnd)
{

    // set up sequence container
    std::string seq = "ATGTTGCCTAATTTCACTATTCAGGGTTATACGCCT";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTTAATTTTATTATTTAGGGTTATATGTTT";
    // reverse sequence
    std::string revSeq = "AGGCGTATAACCCTGAATAGTGAAATTAGGCAACAT";
    std::string redRevSeq = "AGGTGTATAATTTTGAATAGTGAAATTAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (cpgTab, cpgStart, genSeq);

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.kmerTable.size(), 0);

    for (unsigned int i = 0; i <= (33 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = (MyConst::READLEN + 3) - i - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        // lookup if kmer is present in hash table
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer kRev = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, kRev.cpg);
        ASSERT_EQ(off, kRev.offset);
    }

    for (unsigned int i = 3; i <= (36 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer k = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, k.cpg);
        ASSERT_EQ(off | STRANDMASK, k.offset);
    }

    // test if nothing else was hashed
    for (unsigned int i = 0; i < ref.kmerTable.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.kmerTable[i].len);

    }
}


// Test hashing of CpG at end of sequence with Ns
TEST(RefGenome_test, simpleAtEndN)
{

    // set up sequence container
    std::string seq = "ATGTTGCCNAATTTCACTATTCAGGGTTATACGCNT";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTNAATTTTATTATTTAGGGTTATATGTNT";
    // reverse sequence
    std::string revSeq = "ANGCGTATAACCCTGAATAGTGAAATTNGGCAACAT";
    std::string redRevSeq = "ANGTGTATAATTTTGAATAGTGAAATTNGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (cpgTab, cpgStart, genSeq);

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.kmerTable.size(), 0);

    for (unsigned int i = 2; i <= (27 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = (MyConst::READLEN + 3) - i - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        // lookup if kmer is present in hash table
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer kRev = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, kRev.cpg);
        ASSERT_EQ(off, kRev.offset);
    }

    for (unsigned int i = 9; i <= (34 - MyConst::KMERLEN); ++i)
    {

        uint16_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = kmerInd[hash % kmerInd.size()]++;
        ASSERT_LE(index + 1, ref.kmerTable[hash % ref.kmerTable.size()].len);
        const struct kmer k = ref.kmerTable[hash % ref.kmerTable.size()].collis[index];
        ASSERT_EQ(0, k.cpg);
        ASSERT_EQ(off | STRANDMASK, k.offset);
    }

    // test if nothing else was hashed
    for (unsigned int i = 0; i < ref.kmerTable.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.kmerTable[i].len);

    }
}
