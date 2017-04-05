

#include "RefGenome.h"


using namespace std;

RefGenome::RefGenome(vector<struct CpG>& cpgTab, vector<struct CpG>& cpgStartTab, vector<vector<char> >& genomeSeq) :
        cpgTable(cpgTab)
    ,   cpgStartTable(cpgStartTab)
    ,   genomeBit()
    ,   kmerTable()
{
    generateBitStrings(genomeSeq);
    generateHashes(genomeSeq);
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


    // size of the hash table is approximately 2 * |CpGs| * hashed values per CpG
    kmerTable.resize(2 * cpgTable.size() * (2 * MyConst::READLEN - MyConst::KMERLEN));

    // hash CpGs from the start
    uint32_t cpgCount = 0;
    for (vector<struct CpG>::const_iterator it = cpgStartTable.begin(); it != cpgStartTable.end(); ++it, ++cpgCount)
    {

        ntHashFirst(genomeSeq[it->chrom], it->pos, cpgCount);
    }

    cpgCount = 0;
    // hash each CpG
    for (vector<struct CpG>::const_iterator it = cpgTable.begin(); it != cpgTable.end(); ++it, ++cpgCount)
    {
        // hash individual kmers around a CpG
        //
        // how long is the rest of the sequence after start of this cpg
        // last term resembles position directly after CpG
        const unsigned int remainderBps = genomeSeq[it->chrom].size() - 1 - (it->pos - MyConst::READLEN);
        // if we can read the full sequence breadth
        if (remainderBps >= (MyConst::READLEN - 2) )
        {

            ntHashChunk(genomeSeq[it->chrom], it->pos, cpgCount);

        } else {

            ntHashLast(genomeSeq[it->chrom], it->pos, remainderBps, cpgCount);
        }

    }


}

void RefGenome::ntHashLast(const vector<char>& seq,const unsigned int& pos, const unsigned int& bpsAfterCpG, const uint32_t& cpg)
{

    // construct corresponding sequence with reduced alphabet
    vector<char> redSeq(MyConst::READLEN + bpsAfterCpG);
    vector<char> redSeqRev(MyConst::READLEN + bpsAfterCpG);

    vector<char>::const_iterator start = seq.begin();
    advance(start, pos);
    vector<char>::const_iterator end = seq.begin();
    advance(end, pos + MyConst::READLEN);

    // position in substring
    unsigned int j = 0;
    // last N before the CpG
    int lasN = -1;
    // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
    for ( ; start != end; ++start, ++j)
    {

        const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
        switch (*start)
        {
            case 'N':
                lasN = j;
                break;

            case 'C':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[j] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[j] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

        }

    }
    // move to one after last N
    ++lasN;
    // move over second half in reverse order
    // reassign end to the position of G
    --end;
    // reassign start to final position
    start = seq.begin();
    advance(start, pos + MyConst::READLEN + bpsAfterCpG - 1);
    j = MyConst::READLEN + bpsAfterCpG - 1;
    // offset where the first N after the CpG is
    int off = MyConst::READLEN + bpsAfterCpG;
    for ( ; start != end; --start, --j)
    {
        const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
        switch (*start)
        {
            case 'N':
                off = j;
                break;

            case 'C':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[j] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[j] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
        }


    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = off - lasN;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }
    char* seqStart = redSeq.data() + lasN;
    char* seqStartRev = redSeqRev.data() + (MyConst::READLEN + bpsAfterCpG - off);

    // initial hash backward
    uint64_t rhVal = ntHash::NTP64(seqStartRev);
    uint16_t kPos = lasN;

    // update kmer table
    ++kmerTable[rhVal % kmerTable.size()].len;
    kmerTable[rhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);

    // initial hash forward
    uint64_t fhVal = ntHash::NTP64(seqStart);
    kPos |= STRANDMASK;

    // update kmer table
    ++kmerTable[fhVal % kmerTable.size()].len;
    kmerTable[fhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);



    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {

        rhVal = ntHash::NTP64(rhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        kPos = lasN + i + 1;
        // update kmer table
        ++kmerTable[rhVal % kmerTable.size()].len;
        kmerTable[rhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);

        fhVal = ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        kPos |= STRANDMASK;
        // update kmer table
        ++kmerTable[fhVal % kmerTable.size()].len;
        kmerTable[fhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);
    }

}

void RefGenome::ntHashFirst(const vector<char>& seq, const unsigned int& posOfCpG, const uint32_t& cpg)
{

    // construct corresponding sequence with reduced alphabet
    vector<char> redSeq(MyConst::READLEN + posOfCpG);
    vector<char> redSeqRev(MyConst::READLEN + posOfCpG);

    vector<char>::const_iterator start = seq.begin();
    vector<char>::const_iterator end = seq.begin();
    advance(end, posOfCpG);

    // position in substring
    unsigned int pos = 0;
    // last N before the CpG
    int lasN = -1;
    // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
    for ( ; start != end; ++start, ++pos)
    {

        const int revPos = MyConst::READLEN + posOfCpG - 1 - pos;
        switch (*start)
        {
            case 'N':
                lasN = pos;
                break;

            case 'C':
                redSeq[pos] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[pos] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[pos] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[pos] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

        }

    }
    // move to one after last N
    ++lasN;
    // move over second half in reverse order
    // reassign end to the position of G
    --end;
    // reassign start to final position
    start = seq.begin();
    advance(start, pos + MyConst::READLEN + posOfCpG - 1);
    pos = MyConst::READLEN + posOfCpG - 1;
    // offset where the first N after the CpG is
    int off = MyConst::READLEN + posOfCpG;
    for ( ; start != end; --start, --pos)
    {
        const int revPos = MyConst::READLEN + posOfCpG - 1 - pos;
        switch (*start)
        {
            case 'N':
                off = pos;
                break;

            case 'C':
                redSeq[pos] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[pos] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[pos] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[pos] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
        }


    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = off - lasN;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }
    char* seqStart = redSeq.data() + lasN;
    char* seqStartRev = redSeqRev.data() + (MyConst::READLEN - posOfCpG - off);

    // initial hash backward
    uint64_t rhVal = ntHash::NTP64(seqStartRev);
    uint16_t kPos = lasN;

    // update kmer table
    ++kmerTable[rhVal % kmerTable.size()].len;
    kmerTable[rhVal % kmerTable.size()].collis.emplace_back(cpg | MyConst::INDMASK, kPos);

    // initial hash forward
    uint64_t fhVal = ntHash::NTP64(seqStart);
    kPos |= STRANDMASK;

    // update kmer table
    ++kmerTable[fhVal % kmerTable.size()].len;
    kmerTable[fhVal % kmerTable.size()].collis.emplace_back(cpg | MyConst::INDMASK, kPos);



    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {

        rhVal = ntHash::NTP64(rhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        kPos = lasN + i + 1;
        // update kmer table
        ++kmerTable[rhVal % kmerTable.size()].len;
        kmerTable[rhVal % kmerTable.size()].collis.emplace_back(cpg | MyConst::INDMASK, kPos);

        fhVal = ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        kPos |= STRANDMASK;
        // update kmer table
        ++kmerTable[fhVal % kmerTable.size()].len;
        kmerTable[fhVal % kmerTable.size()].collis.emplace_back(cpg | MyConst::INDMASK, kPos);
    }

}
