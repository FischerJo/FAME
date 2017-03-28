

#include "RefGenome.h"


using namespace std;

RefGenome::RefGenome(vector<struct CpG>& cpgTab, vector<const char*>& genSeq, vector<size_t> genSeqLen) :
        cpgTable(cpgTab)
    ,   genomeSeq(genSeq)
    ,   genomeSeqLen(genSeqLen)
    ,   genomeBit()
    ,   kmerTable()
{
}


void RefGenome::generateHashes()
{


    // size of the hash table is approximately |CpGs| * hashed values per CpG
    kmerTable.resize(cpgTable.size() * (2 * MyConst::READLEN - MyConst::KMERLEN));
    // hash each CpG
    for (vector<struct CpG>::const_iterator it = cpgTable.begin(); it != cpgTable.end(); ++it)
    {
        // hash individual kmers
        //
        // offset where to start reading chars
        const int off = it->pos - MyConst::READLEN + 2;
        // CpG is not too near to the start
        if (off >= 0)
        {
            // how long is the rest of the sequence after start of this cpg
            const unsigned int remainderBps = genomeSeqLen[it->chrom] - it->pos;
            // if we can read the full sequence breadth
            if (remainderBps >= MyConst::READLEN)
            {
                ntHashChunk(genomeSeq[it->chrom] + off, *it);
            } else {

                // substract 2 because seq points to CG which should be excluded of the count
                ntHashLast(genomeSeq[it->chrom] + off, remainderBps - 2, *it);
            }
        // CpG is near to the start of chromosome sequence
        } else {

            // ASSUMPTION: chromosome length is bigger than readlength
            ntHashFirst(genomeSeq[it->chrom] + off, it->pos, *it);
        }

    }


}

void RefGenome::ntHashLast(const char* seq, const unsigned int bpsAfterCpG, const struct CpG& cpg)
{

    // find last occurence of N before CpG and first occurence after
    const char* lastBef;
    for (lastBef = seq; (lastBef = (char*) memchr(lastBef, 'N', (MyConst::READLEN - 1 - (lastBef - seq)))); ++lastBef)
    {}
    // now lastBef points to character one after the last N before the CpG

    char* firstAft = (char*) memchr(seq + MyConst::READLEN, 'N', bpsAfterCpG);

    // No N found before CpG, set char* to start of argument sequence
    if (lastBef == NULL)
    {
        lastBef = seq;
    }
    // how many chars to consider, init as full sequence around CpG
    std::size_t offset = MyConst::READLEN + bpsAfterCpG - (lastBef - seq);
    // set offset until first N after CpG, if present
    if (firstAft != NULL)
    {

        offset = firstAft - lastBef;
    }
    // look if enough sequence is there to hash something
    if (offset < MyConst::KMERLEN)
    {
        return;
    }

    // how often should we apply the rolling
    // if zero, make just initial hash
    unsigned int rolls = offset - MyConst::KMERLEN;
    uint64_t fhVal = 0;
    uint64_t rhVal = 0;
    ntHash::NTPC64(lastBef, fhVal, rhVal);
    // note that READLEN is bounded by a small integer (less then 16 bits) and
    // kmer offset is bounded by readlen so pointer arithmetic is fine here
    uint16_t kPos = lastBef - seq;
    kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    // set forward strand flag
    kPos |= 0x4000;
    kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);


    for (unsigned int i = 0; i < rolls; ++i)
    {

        ntHash::NTPC64(lastBef[i], lastBef[MyConst::KMERLEN + i], fhVal, rhVal);
        kPos = lastBef - seq + i + 1;
        kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
        // set forward strand flag
        kPos |= 0x4000;
        kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    }



}

void RefGenome::ntHashFirst(const char* seq, const unsigned int posOfCpG, const struct CpG& cpg)
{

    // find last occurence of N before CpG and first occurence after
    const char* lastBef;
    for (lastBef = seq; (lastBef = (char*) memchr(lastBef, 'N', posOfCpG - (lastBef - seq))); ++lastBef)
    {}
    // now lastBef points to character one after the last N before the CpG

    char* firstAft = (char*) memchr(seq + posOfCpG + 2, 'N', MyConst::READLEN - 2);

    // No N found before CpG, set char* to start of argument sequence
    if (lastBef == NULL)
    {
        lastBef = seq;
    }
    // how many chars to consider, init as full sequence around CpG
    std::size_t offset = MyConst::READLEN + posOfCpG - (lastBef - seq);
    // set offset until first N after CpG, if present
    if (firstAft != NULL)
    {

        offset = firstAft - lastBef;
    }
    // look if enough sequence is there to hash something
    if (offset < MyConst::KMERLEN)
    {
        return;
    }

    // how often should we apply the rolling
    // if zero, make just initial hash
    unsigned int rolls = offset - MyConst::KMERLEN;
    uint64_t fhVal = 0;
    uint64_t rhVal = 0;
    ntHash::NTPC64(lastBef, fhVal, rhVal);
    // note that READLEN is bounded by a small integer (less then 16 bits) and
    // kmer offset is bounded by readlen so pointer arithmetic is fine here
    uint16_t kPos = lastBef - seq;
    kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    // set forward strand flag
    kPos |= 0x4000;
    kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);


    for (unsigned int i = 0; i < rolls; ++i)
    {

        ntHash::NTPC64(lastBef[i], lastBef[MyConst::KMERLEN + i], fhVal, rhVal);
        kPos = lastBef - seq + i + 1;
        kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
        // set forward strand flag
        kPos |= 0x4000;
        kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    }

}
