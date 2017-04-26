
#include "gtest/gtest.h"


#include "RefGenome.h"


// TO RUN THE FOLLOWING TESTS
// MAKE REFGENOME FUNCTIONS PUBLIC
//
// SET READLEN TO 30
//     KMERLEN TO 20
//     HTABSIZE TO 2 << 17


TEST(RefGenome_test, kmerDef)
{

    KMER::kmer k = std::move(KMER::constructKmer(1, 3, 7));

    ASSERT_EQ(3, KMER::getMetaCpG(k));
    ASSERT_EQ(1, KMER::isStartCpG(k));
    ASSERT_EQ(7, KMER::getOffset(k));

}


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

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    // test if metaCpG is constructed correctly
    ASSERT_EQ(0, ref.metaStartCpGs.size());
    ASSERT_EQ(1, ref.metaCpGs.size());
    ASSERT_EQ(0, ref.metaCpGs[0].start);
    ASSERT_EQ(0, ref.metaCpGs[0].end);

    // test if hashTable is initialized correctly
    ASSERT_EQ(2*(2*MyConst::READLEN - MyConst::KMERLEN - 1), ref.kmerTable.size());
    ASSERT_EQ(2*(2*MyConst::READLEN - MyConst::KMERLEN - 1), ref.strandTable.size());

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.tabIndex.size(), 0);

    // count the hashes
    for (unsigned int i = 5; i <= (63 - MyConst::KMERLEN); ++i)
    {


        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }

    unsigned long sum = 0;
    // test if hash offsets are correct
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        sum += kmerInd[i];
        kmerInd[i] = sum;

    }

    // test if everything is placed correctly
    for (unsigned int i = 5; i <= (63 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = (2*MyConst::READLEN - 2) - i + 5 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }

    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.tabIndex[i]);

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

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);



    // test if metaCpG is constructed correctly
    ASSERT_EQ(0, ref.metaStartCpGs.size());
    ASSERT_EQ(1, ref.metaCpGs.size());
    ASSERT_EQ(0, ref.metaCpGs[0].start);
    ASSERT_EQ(0, ref.metaCpGs[0].end);

    // test if hashTable is initialized correctly
    ASSERT_EQ(2*21, ref.kmerTable.size());
    ASSERT_EQ(2*21, ref.strandTable.size());

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.tabIndex.size(), 0);

    // count the hashes
    for (unsigned int i = 16; i <= (56 - MyConst::KMERLEN); ++i)
    {


        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 10; i <= (50 - MyConst::KMERLEN); ++i)
    {

        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }

    unsigned long sum = 0;
    // test if hash offsets are correct
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        sum += kmerInd[i];
        kmerInd[i] = sum;

    }

    // test if everything is placed correctly
    for (unsigned int i = 16; i <= (56 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = (2*MyConst::READLEN - 2) - i + 5 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }

    for (unsigned int i = 10; i <= (50 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.tabIndex[i]);

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

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    ASSERT_EQ(0, ref.metaStartCpGs.size());
    ASSERT_EQ(1, ref.metaCpGs.size());

    // test if nothing else was hashed
    ASSERT_EQ(0, ref.kmerTable.size());
    ASSERT_EQ(0, ref.strandTable.size());
}


// Test with multiple non overlapping CpGs in one Meta CpG
TEST(RefGenome_test, multiCpGNoOverlap)
{
    // set up sequence container
    //                    |                           ||                           |  |                           ||                           |
    std::string seq = "ATGTTGCCTAATTTCACTATTCAGGGTTATACGCCTGGAATATTCTAGGATTCCTAGTCAATTTATGCCATTAGATACTAGCTGTGACCCTCGAGCTGCTGGGAGCTATTGCATGGGTAGT";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTTAATTTTATTATTTAGGGTTATATGTTTGGAATATTTTAGGATTTTTAGTTAATTTATGTTATTAGATATTAGTTGTGATTTTTGAGTTGTTGGGAGTTATTGTATGGGTAGT";
    // reverse sequence
    //                    |                           ||                           |  |                           ||                           |
    std::string revSeq = "ACTACCCATGCAATAGCTCCCAGCAGCTCGAGGGTCACAGCTAGTATCTAATGGCATAAATTGACTAGGAATCCTAGAATATTCCAGGCGTATAACCCTGAATAGTGAAATTAGGCAACAT";
    std::string redRevSeq = "ATTATTTATGTAATAGTTTTTAGTAGTTTGAGGGTTATAGTTAGTATTTAATGGTATAAATTGATTAGGAATTTTAGAATATTTTAGGTGTATAATTTTGAATAGTGAAATTAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);
    cpgTab.emplace_back(0, 63);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    // test if metaCpG is constructed correctly
    ASSERT_EQ(0, ref.metaStartCpGs.size());
    ASSERT_EQ(1, ref.metaCpGs.size());
    ASSERT_EQ(0, ref.metaCpGs[0].start);
    ASSERT_EQ(1, ref.metaCpGs[0].end);

    // test if hashTable is initialized correctly
    ASSERT_EQ(4*(2*MyConst::READLEN - MyConst::KMERLEN - 1), ref.kmerTable.size());
    ASSERT_EQ(4*(2*MyConst::READLEN - MyConst::KMERLEN - 1), ref.strandTable.size());

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.tabIndex.size(), 0);

    // count the hashes
    for (unsigned int i = 60; i <= (118 - MyConst::KMERLEN); ++i)
    {
        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 0; i <= (58 - MyConst::KMERLEN); ++i)
    {
        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 63; i <= (121 - MyConst::KMERLEN); ++i)
    {
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }

    unsigned long sum = 0;
    // test if hash offsets are correct
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        sum += kmerInd[i];
        kmerInd[i] = sum;

    }

    // test if everything is placed correctly
    for (unsigned int i = 60; i <= (118 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = (2*MyConst::READLEN - 2) - i + 60 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }
    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }

    for (unsigned int i = 0; i <= (58 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = 60 + (2*MyConst::READLEN - 2) - i - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }
    for (unsigned int i = 63; i <= (121 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.tabIndex[i]);

    }

}


// Test with multiple overlapping CpGs in one Meta CpG
TEST(RefGenome_test, multiCpGOverlap)
{
    // set up sequence container
    //                    |     |                     ||    ||                     |     |
    std::string seq = "ATGTTGCCTAATTTCACTATTCAGGGTTATACGCCTGCGAATATTCTAGGATTCCTAGTCAATTTAT";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTTAATTTTATTATTTAGGGTTATATGTTTGTGAATATTTTAGGATTTTTAGTTAATTTAT";
    // reverse sequence
    std::string revSeq = "ATAAATTGACTAGGAATCCTAGAATATTCGCAGGCGTATAACCCTGAATAGTGAAATTAGGCAACAT";
    std::string redRevSeq = "ATAAATTGATTAGGAATTTTAGAATATTTGTAGGTGTATAATTTTGAATAGTGAAATTAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);
    cpgTab.emplace_back(0, 9);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    // test if metaCpG is constructed correctly
    ASSERT_EQ(0, ref.metaStartCpGs.size());
    ASSERT_EQ(1, ref.metaCpGs.size());
    ASSERT_EQ(0, ref.metaCpGs[0].start);
    ASSERT_EQ(1, ref.metaCpGs[0].end);

    // test if hashTable is initialized correctly
    ASSERT_EQ(2*(2*MyConst::READLEN - MyConst::KMERLEN - 1) + 12, ref.kmerTable.size());

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.tabIndex.size(), 0);

    // count the hashes
    for (unsigned int i = 3; i <= (67 - MyConst::KMERLEN); ++i)
    {
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 0; i <= (64 - MyConst::KMERLEN); ++i)
    {
        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }

    unsigned long sum = 0;
    // test if hash offsets are correct
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        sum += kmerInd[i];
        kmerInd[i] = sum;

    }

    // test if everything is placed correctly
    for (unsigned int i = 6; i <= (64 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = (2*MyConst::READLEN - 2) - i + 6 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
    }
    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }


    for (unsigned int i = 0; i < 6; ++i)
    {

        uint32_t off = 6 + (2*MyConst::READLEN - 2) - i - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }
    for (unsigned int i = (62 - MyConst::KMERLEN); i <= (67 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.tabIndex[i]);

    }

}

// Test with multiple Meta CpGs (on two chromosomes)
TEST(RefGenome_test, multiMeta)
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
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);
    cpgTab.emplace_back(1, 3);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    // test if metaCpG is constructed correctly
    ASSERT_EQ(0, ref.metaStartCpGs.size());
    ASSERT_EQ(2, ref.metaCpGs.size());
    ASSERT_EQ(0, ref.metaCpGs[0].start);
    ASSERT_EQ(0, ref.metaCpGs[0].end);
    ASSERT_EQ(1, ref.metaCpGs[1].start);
    ASSERT_EQ(1, ref.metaCpGs[1].end);

    // test if hashTable is initialized correctly
    ASSERT_EQ(4*(2*MyConst::READLEN - MyConst::KMERLEN - 1), ref.kmerTable.size());
    ASSERT_EQ(4*(2*MyConst::READLEN - MyConst::KMERLEN - 1), ref.strandTable.size());

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.tabIndex.size(), 0);

    // count the hashes
    for (unsigned int i = 5; i <= (63 - MyConst::KMERLEN); ++i)
    {


        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        kmerInd[hash % kmerInd.size()] += 2;
    }
    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        kmerInd[hash % kmerInd.size()] += 2;
    }

    unsigned long sum = 0;
    // test if hash offsets are correct
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        sum += kmerInd[i];
        kmerInd[i] = sum;

    }

    // test if everything is placed correctly for first meta
    for (unsigned int i = 5; i <= (63 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = (2*MyConst::READLEN - 2) - i + 5 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }

    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    // test if everything is placed correctly for second meta
    for (unsigned int i = 5; i <= (63 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = (2*MyConst::READLEN - 2) - i + 5 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(1, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }

    for (unsigned int i = 3; i <= (61 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(1, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.tabIndex[i]);

    }

}


// Test hashing of CpG at end of sequence
TEST(RefGenome_test, simpleAtEnd)
{
    // set up sequence container
    std::string seq = "ATGTTGCCTAATTTCACTATTCAGGGTTATACGCCTGGAA";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTTAATTTTATTATTTAGGGTTATATGTTTGGAA";
    // reverse sequence
    std::string revSeq = "TTCCAGGCGTATAACCCTGAATAGTGAAATTAGGCAACAT";
    std::string redRevSeq = "TTTTAGGTGTATAATTTTGAATAGTGAAATTAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    // test if metaCpG is constructed correctly
    ASSERT_EQ(0, ref.metaStartCpGs.size());
    ASSERT_EQ(1, ref.metaCpGs.size());
    ASSERT_EQ(0, ref.metaCpGs[0].start);
    ASSERT_EQ(0, ref.metaCpGs[0].end);

    // test if hashTable is initialized correctly
    ASSERT_EQ(36, ref.kmerTable.size());
    ASSERT_EQ(36, ref.strandTable.size());

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.tabIndex.size(), 0);

    // count the hashes
    for (unsigned int i = 0; i <= (37 - MyConst::KMERLEN); ++i)
    {


        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 3; i <= (40 - MyConst::KMERLEN); ++i)
    {

        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }

    unsigned long sum = 0;
    // test if hash offsets are correct
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        sum += kmerInd[i];
        kmerInd[i] = sum;

    }

    // test if everything is placed correctly
    for (unsigned int i = 0; i <= (37 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = (MyConst::READLEN + 7) - i - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }

    for (unsigned int i = 3; i <= (40 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i - 3;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.tabIndex[i]);

    }

}


// Test hashing of CpG at start of CpG
TEST(RefGenome_test, simpleAtStart)
{

    // set up sequence container
    //                 |         ||                           |
    std::string seq = "CAGGGTTATACGCCTGGAATATTCTAGGATTCCTAGTCAATTTAT";
    // sequence with reduced alphabet
    std::string redSeq = "TAGGGTTATATGTTTGGAATATTTTAGGATTTTTAGTTAATTTAT";
    // reverse sequence
    std::string revSeq = "ATAAATTGACTAGGAATCCTAGAATATTCCAGGCGTATAACCCTG";
    std::string redRevSeq = "ATAAATTGATTAGGAATTTTAGAATATTTTAGGTGTATAATTTTG";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;

    std::vector<struct CpG> cpgStart;
    cpgStart.emplace_back(0, 10);

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    // test if metaCpG is constructed correctly
    ASSERT_EQ(1, ref.metaStartCpGs.size());
    ASSERT_EQ(0, ref.metaStartCpGs[0].start);
    ASSERT_EQ(0, ref.metaStartCpGs[0].end);
    ASSERT_EQ(0, ref.metaCpGs.size());

    // test if hashTable is initialized correctly
    ASSERT_EQ(42, ref.kmerTable.size());
    ASSERT_EQ(42, ref.strandTable.size());

    // index to handle collisions in test
    // position i safes how many times we already accessed corresponding struct in kmerTable
    std::vector<unsigned int> kmerInd (ref.tabIndex.size(), 0);

    // count the hashes
    for (unsigned int i = 5; i <= (45 - MyConst::KMERLEN); ++i)
    {


        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }
    for (unsigned int i = 0; i <= (40 - MyConst::KMERLEN); ++i)
    {

        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        ++kmerInd[hash % kmerInd.size()];
    }

    unsigned long sum = 0;
    // test if hash offsets are correct
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        sum += kmerInd[i];
        kmerInd[i] = sum;

    }

    // test if everything is placed correctly
    for (unsigned int i = 5; i <= (45 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = 40 - i + 5 - MyConst::KMERLEN;

        // get hash of corresponding reverse kmer
        uint64_t hash = ntHash::NTP64(redRevSeq.data() + i);
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        // lookup if kmer is present in hash table
        bool strand = ref.strandTable[index];
        ASSERT_EQ(0, strand);
        KMER::kmer kRev = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(kRev));
        ASSERT_EQ(1, KMER::isStartCpG(kRev));
        ASSERT_EQ(off, KMER::getOffset(kRev));
    }

    for (unsigned int i = 0; i <= (40 - MyConst::KMERLEN); ++i)
    {

        uint32_t off = i;
        // get hash of corresponding kmer
        uint64_t hash = ntHash::NTP64(redSeq.data() + i);
        // lookup if kmer is present in hash table
        unsigned int index = --kmerInd[hash % ref.tabIndex.size()];
        KMER::kmer k = ref.kmerTable[index];
        ASSERT_EQ(0, KMER::getMetaCpG(k));
        ASSERT_EQ(off, KMER::getOffset(k));
        ASSERT_EQ(1, KMER::isStartCpG(k));
        bool strand = ref.strandTable[index];
        ASSERT_EQ(1, strand);
    }
    for (unsigned int i = 0; i < ref.tabIndex.size(); ++i)
    {

        ASSERT_EQ(kmerInd[i], ref.tabIndex[i]);

    }
}
// TODO atStartN atEndN

TEST(RefGenome_test, simpleGetSeed)
{
    // set up sequence container
    //                    |                           ||                           |  |                           ||                           |
    std::string seq = "ATGTTGCCTAATTTCACTATTCAGGGTTATACGCCTGGAATATTCTAGGATTCCTAGTCAATTTATGCCATTAGATACTAGCTGTGACCCTCGAGCTGCTGGGAGCTATTGCATGGGTAGT";
    // sequence with reduced alphabet
    std::string redSeq = "ATGTTGTTTAATTTTATTATTTAGGGTTATATGTTTGGAATATTTTAGGATTTTTAGTTAATTTATGTTATTAGATATTAGTTGTGATTTTTGAGTTGTTGGGAGTTATTGTATGGGTAGT";
    // reverse sequence
    //                    |                           ||                           |  |                           ||                           |
    std::string revSeq = "ACTACCCATGCAATAGCTCCCAGCAGCTCGAGGGTCACAGCTAGTATCTAATGGCATAAATTGACTAGGAATCCTAGAATATTCCAGGCGTATAACCCTGAATAGTGAAATTAGGCAACAT";
    std::string redRevSeq = "ATTATTTATGTAATAGTTTTTAGTAGTTTGAGGGTTATAGTTAGTATTTAATGGTATAAATTGATTAGGAATTTTAGAATATTTTAGGTGTATAATTTTGAATAGTGAAATTAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    cpgTab.emplace_back(0, 3);
    cpgTab.emplace_back(0, 63);

    std::vector<struct CpG> cpgStart;

    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);

    // meta 0 off 3
    // CpG: 0 off 3
    // fwd
    std::string read1 = "TTTAATTTTATTATTTAGGGTTATATGTTT";
    // meta 0 off 61
    // CpG: 1 off 1
    // fwd
    std::string read2 = "ATGTTATTAGATATTAGTTGTGATTTTTGA";

    // meta 0 off 78
    // CpG 1 off 18
    // rev
    std::string read3 = "TAATAGTTTTTAGTAGTTTGAGGGTTATAG";

    // test matching first read
    std::vector<std::vector<KMER::kmer> > kmers1;
    std::vector<std::vector<bool> > strands1;
    std::vector<char> read1V(read1.begin(), read1.end());
    ref.getSeeds(read1V, kmers1, strands1);
    // check that for each read kmer we have an entry
    ASSERT_EQ(11, kmers1.size());
    ASSERT_EQ(kmers1.size(), strands1.size());
    for (unsigned int i = 0; i < kmers1.size(); ++i)
    {
        // check if each kmer is matched
        ASSERT_LE(1, kmers1[i].size());
        ASSERT_LE(1, strands1[i].size());

        bool matchflag = false;
        // check if there is one match that is the correct one
        for (unsigned int j = 0; j < kmers1[i].size(); ++j)
        {

            // check if kmer in seeds is the same as query read
            if (KMER::getOffset(kmers1[i][j]) == 3 + i && strands1[i][j] == true && KMER::getMetaCpG(kmers1[i][j]) == 0)
            {
                matchflag = true;
            }
        }
        ASSERT_EQ(true, matchflag);
    }

    // test matching second read
    std::vector<std::vector<KMER::kmer> > kmers2;
    std::vector<std::vector<bool> > strands2;
    std::vector<char> read2V(read2.begin(), read2.end());
    ref.getSeeds(read2V, kmers2, strands2);
    // check that for each read kmer we have an entry
    ASSERT_EQ(11, kmers2.size());
    ASSERT_EQ(kmers2.size(), strands2.size());
    for (unsigned int i = 0; i < kmers2.size(); ++i)
    {
        // check if each kmer is matched
        ASSERT_LE(1, kmers2[i].size());
        ASSERT_LE(1, strands2[i].size());

        bool matchflag = false;
        // check if there is one match that is the correct one
        for (unsigned int j = 0; j < kmers2[i].size(); ++j)
        {

            // check if kmer in seeds is the same as query read
            if (KMER::getOffset(kmers2[i][j]) == 61 + i && strands2[i][j] == true && KMER::getMetaCpG(kmers2[i][j]) == 0)
            {
                matchflag = true;
            }
        }
        ASSERT_EQ(true, matchflag);
    }

    // test matching third read
    std::vector<std::vector<KMER::kmer> > kmers3;
    std::vector<std::vector<bool> > strands3;
    std::vector<char> read3V(read3.begin(), read3.end());
    ref.getSeeds(read3V, kmers3, strands3);
    // check that for each read kmer we have an entry
    ASSERT_EQ(11, kmers3.size());
    ASSERT_EQ(kmers3.size(), strands3.size());
    for (unsigned int i = 0; i < kmers3.size(); ++i)
    {
        // check if each kmer is matched
        ASSERT_LE(1, kmers3[i].size());
        ASSERT_LE(1, strands3[i].size());

        bool matchflag = false;
        // check if there is one match that is the correct one
        for (unsigned int j = 0; j < kmers3[i].size(); ++j)
        {

            // check if kmer in seeds is the same as query read
            if (KMER::getOffset(kmers3[i][j]) == 78 + (10 - i) && strands3[i][j] == false && KMER::getMetaCpG(kmers3[i][j]) == 0)
            {
                matchflag = true;
            }
        }
        ASSERT_EQ(true, matchflag);
    }
}


TEST(RefGenome_test, simpleGetSeedStart)
{

}
