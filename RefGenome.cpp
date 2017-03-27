

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
        // offset where to start reading chars
        const int off = it->pos - MyConst::READLEN + 2;
        // CpG is not too near to the start
        if (off >= 0))
        {
            // if we can read the full sequence breadth
            if ((genomeSeqLen[it->chrom] - it->pos) >= MyConst::READLEN)
            {
                ntHashChunk(genomeSeq[it->chrom] + off);
            } else {

                // TODO
                ntHashLast(genomeSeq[it->chrom], bpsAfterCpG);
            }
        } else {

            // ASSUMPTION: chromosome length is bigger than readlength
            // TODO
            ntHasFirst(genomeSeq[it->chrom], posOfCpG);
        }

    }


}
