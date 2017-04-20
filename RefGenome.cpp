

#include "RefGenome.h"


using namespace std;

RefGenome::RefGenome(vector<struct CpG>&& cpgTab, vector<struct CpG>&& cpgStartTab, vector<vector<char> >& genomeSeq) :
        cpgTable(cpgTab)
    ,   cpgStartTable(cpgStartTab)
    ,   genomeBit()
    ,   tabIndex(MyConst::HTABSIZE, 0)
    ,   kmerTable()
    ,   strandTable()
    ,   metaCpGs()
    ,   metaStartCpGs()
{
    // find out genome size
    unsigned int gensize = 0;
    for (vector<char> chr : genomeSeq)
    {
        gensize += chr.size();
    }
    // init meta table with upper bound on required windows
    metaCpGs.reserve(gensize/MyConst::WINLEN);
    metaStartCpGs.reserve(MyConst::CHROMNUM);
    // fill meta table
    generateMetaCpGs();
    // generate encoding of genome
    generateBitStrings(genomeSeq);
    cout << "Done generating Genome bit representation" << endl;
    // hash all kmers of reduced alphabet
    generateHashes(genomeSeq);
    cout << "Done hashing CpGs" << endl;
}

void RefGenome::generateMetaCpGs()
{

    uint32_t cpgStartInd = 0;
    // current chromosome we are looking at
    unsigned int currChr = 0;
    // start and end index for CpGs in current window
    uint32_t start = 0;
    // Assuming window is larger then readlength
    while (cpgStartInd < cpgStartTable.size())
    {

        while (cpgStartInd < cpgStartTable.size() && cpgStartTable[cpgStartInd].chrom == currChr)
        {

            ++cpgStartInd;
        }

        metaStartCpGs.emplace_back(start, cpgStartInd - 1);
        ++currChr;
        start = cpgStartInd;
    }


    uint32_t cpgTabInd = 0;
    start = 0;
    currChr = 0;

    // while there are CpGs left to look at
    //
    while (cpgTabInd < cpgTable.size())
    {

        // fill current metaCpG window with CpGs
        const uint32_t windowEnd = cpgTable[start].pos + MyConst::WINLEN;

        for (; cpgTabInd < cpgTable.size(); ++cpgTabInd)
        {
            if (cpgTable[cpgTabInd].pos >= windowEnd || cpgTable[cpgTabInd].chrom != currChr)
            {

                // generate meta CpG
                metaCpGs.emplace_back(start, cpgTabInd - 1);
                // set start for next meta CpG
                start = cpgTabInd;
                // update current chromosome index if necessary
                if (cpgTable[cpgTabInd].chrom != currChr) ++currChr;
                // break out of inner for loop (i.e. filling of current metaCpG
                break;
            }

        }

    }
    // put in meta CpG containging last CpGs
    metaCpGs.emplace_back(start, cpgTabInd - 1);

    metaCpGs.shrink_to_fit();
    metaStartCpGs.shrink_to_fit();
}


void RefGenome::generateBitStrings(vector<vector<char> >& genomeSeq)
{

    unsigned int genSeqIndex = 0;
    // generate genome encodings with bitmasks
    for (vector<char>& seq : genomeSeq)
    {
        const unsigned int segLen = seq.size() / 32;
        genomeBit.emplace_back(seq.size());

        for (unsigned int i = 0; i < segLen; ++i)
        {

            genomeBit[genSeqIndex].setBitStrN( string(seq.data() + (i*32), 32), i);
        }
        if ((seq.size() % 32) != 0)
        {

            genomeBit[genSeqIndex].setBitStrLast( string(seq.data() + (segLen*32), seq.size() % 32));
        }
        ++genSeqIndex;
    }

}

void RefGenome::generateHashes(vector<vector<char> >& genomeSeq)
{


    estimateTablesizes(genomeSeq);

    // cout << "\nKmer table size: " << kmerTable.size() << endl;

    // hash CpGs from the start
    uint32_t cpgCount = 0;
    for (vector<struct metaCpG>::const_iterator it = metaStartCpGs.begin(); it != metaStartCpGs.end(); ++it, ++cpgCount)
    {

        unsigned int lastPos = 0;
        for (uint32_t cpgInd = it->start; cpgInd <= it->end; ++cpgInd)
        {

            ntHashFirst(genomeSeq[cpgStartTable[it->start].chrom], lastPos, cpgStartTable[cpgInd].pos, cpgCount);
        }
    }

    cpgCount = 0;
    // hash each CpG
    for (vector<struct metaCpG>::const_iterator it = metaCpGs.begin(); it != metaCpGs.end(); ++it, ++cpgCount)
    {

        unsigned int lastPos = 0;

        for (uint32_t cpgInd = it->start; cpgInd <= it->end; ++cpgInd)
        {

            // how long is the rest of the sequence after the end of last cpg in meta cpg
            const unsigned int remainderBps = genomeSeq[cpgTable[cpgInd].chrom].size() - (cpgTable[cpgInd].pos + MyConst::READLEN);
            // if we can read the full sequence breadth after the CpG
            if (remainderBps >= (MyConst::READLEN - 2) )
            {

                ntHashChunk(genomeSeq[cpgTable[it->start].chrom], lastPos, cpgTable[cpgInd].pos, cpgCount, cpgTable[cpgInd].pos - cpgTable[it->start].pos);


            } else {

                ntHashLast(genomeSeq[cpgTable[it->start].chrom], lastPos, cpgStartTable[cpgInd].pos, remainderBps, cpgCount, cpgTable[cpgInd].pos - cpgTable[it->start].pos);
            }
        }
    }

}

void RefGenome::ntHashLast(const vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG, const uint32_t& metacpg, uint32_t&& metaOff)
{
    //
    // // construct corresponding sequence with reduced alphabet
    // vector<char> redSeq(MyConst::READLEN + bpsAfterCpG);
    // vector<char> redSeqRev(MyConst::READLEN + bpsAfterCpG);
    //
    // vector<char>::const_iterator start = seq.begin();
    // advance(start, pos);
    // vector<char>::const_iterator end = seq.begin();
    // advance(end, pos + MyConst::READLEN);
    //
    // // position in substring
    // unsigned int j = 0;
    // // last N before the CpG
    // int lasN = -1;
    // // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
    // for ( ; start != end; ++start, ++j)
    // {
    //
    //     const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
    //     switch (*start)
    //     {
    //         case 'N':
    //             lasN = j;
    //             break;
    //
    //         case 'C':
    //             redSeq[j] = 'T';
    //             redSeqRev[revPos] = 'G';
    //             break;
    //
    //         case 'G':
    //             redSeq[j] = 'G';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         case 'T':
    //             redSeq[j] = 'T';
    //             redSeqRev[revPos] = 'A';
    //             break;
    //
    //         case 'A':
    //             redSeq[j] = 'A';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         default:
    //             cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
    //
    //     }
    //
    // }
    // // move to one after last N
    // ++lasN;
    // // move over second half in reverse order
    // // reassign end to the position of G
    // --end;
    // // reassign start to final position
    // start = seq.begin();
    // advance(start, pos + MyConst::READLEN + bpsAfterCpG - 1);
    // j = MyConst::READLEN + bpsAfterCpG - 1;
    // // offset where the first N after the CpG is
    // int off = MyConst::READLEN + bpsAfterCpG;
    // for ( ; start != end; --start, --j)
    // {
    //     const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
    //     switch (*start)
    //     {
    //         case 'N':
    //             off = j;
    //             break;
    //
    //         case 'C':
    //             redSeq[j] = 'T';
    //             redSeqRev[revPos] = 'G';
    //             break;
    //
    //         case 'G':
    //             redSeq[j] = 'G';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         case 'T':
    //             redSeq[j] = 'T';
    //             redSeqRev[revPos] = 'A';
    //             break;
    //
    //         case 'A':
    //             redSeq[j] = 'A';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         default:
    //             cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
    //     }
    //
    //
    // }
    //
    // // length of the context around the CpG without Ns
    // unsigned int contextLen = off - lasN;
    // // if we don't have enough to read, return without hashing
    // if (contextLen < MyConst::KMERLEN)
    // {
    //     return;
    // }
    // char* seqStart = redSeq.data() + lasN;
    // char* seqStartRev = redSeqRev.data() + (MyConst::READLEN + bpsAfterCpG - off);
    //
    // // initial hash backward
    // uint64_t rhVal = ntHash::NTP64(seqStartRev);
    // // first kmer on reverse complement corresponds to last kmer in forward sequence
    // uint16_t kPosRev = off - MyConst::KMERLEN;
    //
    // // update kmer table
    // kmerTable[rhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(0, cpg, kPosRev, 0)));
    //
    // // hash kmers of backward strand
    // for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    // {
    //     --kPosRev;
    //     ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
    //     // update kmer table
    //     kmerTable[rhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(0, cpg, kPosRev, 0)));
    //
    // }
    //
    // // initial hash forward
    // uint64_t fhVal = ntHash::NTP64(seqStart);
    // uint16_t kPos = lasN;
    //
    // // update kmer table
    // kmerTable[fhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(1, cpg, kPos, 0)));
    //
    // for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    // {
    //
    //     ++kPos;
    //     ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
    //     // update kmer table
    //     kmerTable[fhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(1, cpg, kPos, 0)));
    // }
}

void RefGenome::ntHashFirst(const vector<char>& seq, uint32_t& lastPos, const unsigned int& posOfCpG, const uint32_t& metacpg)
{
    //
    // // construct corresponding sequence with reduced alphabet
    // vector<char> redSeq(MyConst::READLEN + posOfCpG);
    // vector<char> redSeqRev(MyConst::READLEN + posOfCpG);
    //
    // vector<char>::const_iterator start = seq.begin();
    // vector<char>::const_iterator end = seq.begin();
    // advance(end, posOfCpG + 2);
    //
    // // position in substring
    // unsigned int pos = 0;
    // // last N before the CpG
    // int lasN = -1;
    // // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
    // for ( ; start != end; ++start, ++pos)
    // {
    //
    //     const int revPos = MyConst::READLEN + posOfCpG - 1 - pos;
    //     switch (*start)
    //     {
    //         case 'N':
    //             lasN = pos;
    //             break;
    //
    //         case 'C':
    //             redSeq[pos] = 'T';
    //             redSeqRev[revPos] = 'G';
    //             break;
    //
    //         case 'G':
    //             redSeq[pos] = 'G';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         case 'T':
    //             redSeq[pos] = 'T';
    //             redSeqRev[revPos] = 'A';
    //             break;
    //
    //         case 'A':
    //             redSeq[pos] = 'A';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         default:
    //             cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
    //
    //     }
    //
    // }
    // // move to one after last N
    // ++lasN;
    // // move over second half in reverse order
    // // reassign end to the position of G
    // --end;
    // // reassign start to final position
    // start = seq.begin();
    // advance(start, MyConst::READLEN + posOfCpG - 1);
    // pos = MyConst::READLEN + posOfCpG - 1;
    // // offset where the first N after the CpG is
    // int off = MyConst::READLEN + posOfCpG;
    // for ( ; start != end; --start, --pos)
    // {
    //     const int revPos = MyConst::READLEN + posOfCpG - 1 - pos;
    //     switch (*start)
    //     {
    //         case 'N':
    //             off = pos;
    //             break;
    //
    //         case 'C':
    //             redSeq[pos] = 'T';
    //             redSeqRev[revPos] = 'G';
    //             break;
    //
    //         case 'G':
    //             redSeq[pos] = 'G';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         case 'T':
    //             redSeq[pos] = 'T';
    //             redSeqRev[revPos] = 'A';
    //             break;
    //
    //         case 'A':
    //             redSeq[pos] = 'A';
    //             redSeqRev[revPos] = 'T';
    //             break;
    //
    //         default:
    //             cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
    //     }
    //
    //
    // }
    //
    // // length of the context around the CpG without Ns
    // unsigned int contextLen = off - lasN;
    // // if we don't have enough to read, return without hashing
    // if (contextLen < MyConst::KMERLEN)
    // {
    //     return;
    // }
    // char* seqStart = redSeq.data() + lasN;
    // char* seqStartRev = redSeqRev.data() + (MyConst::READLEN + posOfCpG - off);
    //
    // // initial hash backward
    // uint64_t rhVal = ntHash::NTP64(seqStartRev);
    // // first kmer on reverse complement corresponds to last kmer in forward sequence
    // uint16_t kPosRev = off - MyConst::KMERLEN;
    //
    // // update kmer table
    // kmerTable[rhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(0, cpg, kPosRev, 1)));
    //
    // // hash kmers of backward strand
    // for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    // {
    //     --kPosRev;
    //     ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
    //     // update kmer table
    //     kmerTable[rhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(0, cpg, kPosRev, 1)));
    //
    // }
    //
    // // initial hash forward
    // uint64_t fhVal = ntHash::NTP64(seqStart);
    // uint16_t kPos = lasN;
    //
    // // update kmer table
    // kmerTable[fhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(1, cpg, kPos, 1)));
    //
    // // hash kmers of forward strand
    // for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    // {
    //
    //
    //     ++kPos;
    //     ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
    //     // update kmer table
    //     kmerTable[fhVal % kmerTable.size()].emplace_back(std::move(KMER::constructKmer(1, cpg, kPos, 1)));
    // }
}

void RefGenome::ntCountFirst(vector<char>& seq, uint32_t& lastPos, const unsigned int& cpgOffset)
{
}
void RefGenome::ntCountLast(vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG)
{
}
