

#include "RefGenome.h"


using namespace std;

RefGenome::RefGenome(vector<struct CpG>& cpgTab, vector<struct CpG>& cpgStartTab, vector<vector<char> >& genSeq) :
        cpgTable(cpgTab)
    ,   cpgStartTable(cpgStartTab)
    ,   genomeSeq(genSeq)
    ,   genomeBit()
    ,   kmerTable()
{
}


void RefGenome::generateHashes()
{


    // size of the hash table is approximately |CpGs| * hashed values per CpG
    kmerTable.resize(cpgTable.size() * (2 * MyConst::READLEN - MyConst::KMERLEN));

    // hash CpGs from the start
    for (vector<struct CpG>::const_iterator it = cpgStartTable.begin(); it != cpgStartTable.end(); ++it)
    {

        ntHashFirst(genomeSeq[it->chrom], it->pos, *it);
    }

    // hash each CpG
    for (vector<struct CpG>::const_iterator it = cpgTable.begin(); it != cpgTable.end(); ++it)
    {
        // hash individual kmers
        //
        // how long is the rest of the sequence after start of this cpg
        // last term resembles position directly after CpG
        const unsigned int remainderBps = genomeSeq[it->chrom].size() - 1 - (it->pos - MyConst::READLEN);
        // if we can read the full sequence breadth
        if (remainderBps >= (MyConst::READLEN - 2) )
        {
            ntHashChunk(genomeSeq[it->chrom], it->pos, *it);
        } else {

            // substract 2 because seq points to CG which should be excluded of the count
            ntHashLast(genomeSeq[it->chrom], it->pos, remainderBps, *it);
        }

    }


}

void RefGenome::ntHashLast(const vector<char>& seq,const unsigned int& pos, const unsigned int& bpsAfterCpG, const struct CpG& cpg)
{

    // retrieve the underlying char vector of the sequence and offset it to the cpg context position
    const char* cpgContext = seq.data() + pos;

    // find last N before the CpG
    unsigned int lasN = 0;
    for ( ; lasN < (MyConst::READLEN - 2); ++lasN)
    {

        if ( *(cpgContext + MyConst::READLEN - 3 - lasN) == 'N')
        {
            break;
        }
    }

    // find first N after cpg
    unsigned int off = 0;
    for ( ; off < bpsAfterCpG; ++off)
    {

        if ( *(cpgContext + MyConst::READLEN + off) == 'N')
        {

            break;
        }
    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = lasN + off + 2;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }
    // move char* to one after the last N occuring before a CpG
    unsigned int start = MyConst::READLEN - 2 - lasN;
    cpgContext += start;

    uint64_t fhVal = 0;
    uint64_t rhVal = 0;
    ntHash::NTPC64(cpgContext, fhVal, rhVal);
    // note that READLEN is bounded by a small integer (less then 16 bits) and
    // kmer offset is bounded by readlen so pointer arithmetic is fine here
    uint16_t kPos = start;
    kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    // set forward strand flag
    kPos |= 0x4000;
    kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);


    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {

        ntHash::NTPC64(cpgContext[i], cpgContext[MyConst::KMERLEN + i], fhVal, rhVal);
        kPos = start + i + 1;
        kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
        // set forward strand flag
        kPos |= 0x4000;
        kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    }
}

void RefGenome::ntHashFirst(const vector<char>& seq, const unsigned int& posOfCpG, const struct CpG& cpg)
{

    // retrieve the underlying char vector of the sequence and offset it to the cpg context position
    const char* cpgContext = seq.data();

    // find last N before the CpG
    unsigned int lasN = 0;
    for ( ; lasN < posOfCpG; ++lasN)
    {

        if ( *(cpgContext + posOfCpG - 1 - lasN) == 'N')
        {
            break;
        }
    }

    // find first N after cpg
    unsigned int off = 0;
    for ( ; off < (MyConst::READLEN - 2); ++off)
    {

        if ( *(cpgContext + posOfCpG + 2 + off) == 'N')
        {

            break;
        }
    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = lasN + off + 2;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }
    // move char* to one after the last N occuring before a CpG
    unsigned int start = posOfCpG - lasN;
    cpgContext += start;

    uint64_t fhVal = 0;
    uint64_t rhVal = 0;
    ntHash::NTPC64(cpgContext, fhVal, rhVal);
    // note that READLEN is bounded by a small integer (less then 16 bits) and
    // kmer offset is bounded by readlen so pointer arithmetic is fine here
    uint16_t kPos = start;
    kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    // set forward strand flag
    kPos |= 0x4000;
    kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);


    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {

        ntHash::NTPC64(cpgContext[i], cpgContext[MyConst::KMERLEN + i], fhVal, rhVal);
        kPos = start + i + 1;
        kmerTable[rhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
        // set forward strand flag
        kPos |= 0x4000;
        kmerTable[fhVal % kmerTable.capacity()].emplace_back(&cpg, kPos);
    }
}
