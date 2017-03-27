

#include "RefGenome.h"
#include "ntHash-1.0.2/nthash.hpp"


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


    for (vector<struct CpG>::const_iterator it = cpgTable.begin(); it != cpgTable.end(); ++it)
    {
        ntHashChunk(genomeSeq[it->chrom] + it->pos, genomeSeqLen[it->chrom]);

    }


}
