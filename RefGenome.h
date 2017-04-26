#ifndef REFGENOME_H
#define REFGENOME_H

#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <list>
#include <cstring>      // memchr()

// Project includes
#include "CONST.h"
#include "structs.h"
#include "DnaBitStr.h"
#include "ntHash-1.0.2/nthash.hpp"


// Class representing the reference genome
class RefGenome
{

    public:

        RefGenome() = delete;

        // Ctor
        // ARGUMENTS: ( see class members for full information )
        //      cpgTab      table of all CpGs in reference genome except for the ones near the start of a sequence (i.e. less then READLEN away from start)
        //      cpgStartTab table of CpGs near the start
        //      genSeq      genomic sequence seperated by chromosome
        RefGenome(std::vector<struct CpG>&& cpgTab, std::vector<struct CpG>&& cpgStartTab, std::vector<std::vector<char> >& genSeq);

        ~RefGenome() = default;


        // Returns the seeds of the reference genome for the given read
        // A seed is a kmer match between a reference genome part and the read
        // both hashed using nthash with the reduced alphabet
        //
        // ARGUMENTS:   seq     read sequence which should be hashed to obtain initial seeds of reference
        //                      NOTE: this sequence should be with reduced alphabet
        //              seedsK  empty vector - contains the kmers of all seeds on return
        //
        //              seedsS  empty vector - contains the strand info of all seeds on return
        //
        //                      outer vector elements i corresponds to kmers with offset i of read
        //                      inner vector gives the seeds in the reference
        //
        inline void getSeeds(std::vector<char>& seq, std::vector<std::vector<KMER::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS)
        {

            seedsK.reserve(seq.size() - MyConst::KMERLEN + 1);
            seedsS.reserve(seq.size() - MyConst::KMERLEN + 1);

            // iterators for kmertable
            std::vector<KMER::kmer>::iterator startK = kmerTable.begin();
            std::vector<KMER::kmer>::iterator endK = kmerTable.begin();
            // iterators for strandtable
            std::vector<bool>::iterator startS = strandTable.begin();
            std::vector<bool>::iterator endS = strandTable.begin();

            // retrieve kmers for first hash
            uint64_t hVal = ntHash::NTP64(seq.data());

            // get iterators framing subset of kmers corresponding to this hash
            std::advance(startK, tabIndex[hVal % tabIndex.size()]);
            std::advance(endK, tabIndex[(hVal % tabIndex.size()) + 1]);
            std::advance(startS, tabIndex[hVal % tabIndex.size()]);
            std::advance(endS, tabIndex[(hVal % tabIndex.size()) + 1]);

            // retrieve seeds
            seedsK.emplace_back(startK, endK);
            seedsS.emplace_back(startS, endS);

            for (unsigned int i = 0; i < (seq.size() - MyConst::KMERLEN); ++i)
            {

                // use rolling hash
                ntHash::NTP64(hVal, seq[i], seq[i + MyConst::KMERLEN]);

                startK = kmerTable.begin();
                std::advance(startK, tabIndex[hVal % tabIndex.size()]);
                endK = kmerTable.begin();
                std::advance(endK, tabIndex[(hVal % tabIndex.size()) + 1]);

                startS = strandTable.begin();
                std::advance(startS, tabIndex[hVal % tabIndex.size()]);
                endS = strandTable.begin();
                std::advance(endS, tabIndex[(hVal % tabIndex.size()) + 1]);

                // retrieve seeds
                seedsK.emplace_back(startK, endK);
                seedsS.emplace_back(startS, endS);
            }

        }



    // TODO: make this private once its tested
    // private:

        // simple max inline functions
        inline int max(int x, int y) {return x > y ? x : y;}

        // produces all struct Meta CpGs
        void generateMetaCpGs();


        // generate Bit representation of whole genome for genomeBit
        // using full alphabet
        void generateBitStrings(std::vector<std::vector<char> >& genomeSeq);


        // hash all kmers in all CpGs to _kmerTable using ntHash
        // the kmers are represented in REDUCED alphabet {A,T,G}
        // mapping all Cs to Ts
        void generateHashes(std::vector<std::vector<char> >& genomeSeq);


        // generates all kmers in seq and hashes them and their reverse complement using nthash into kmerTable
        // using the reduced alphabet {A,G,T}
        // Arguments:
        //          seq         sequence to look at
        //          lastPos     position of the last kmer hashed inside meta cpg + 1
        //                      THIS VALUE WILL BE UPDATED BY THE FUNCTION
        //          pos         position of CpG in sequence (look at CpG struct for details)
        //          metaCpG     index of meta CpG that we are looking at
        //          metaOff     offset of pos in metaCpG
        //
        inline void ntHashChunk(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const uint32_t& metacpg, const uint32_t&& metaOff)
        {

            // kmers to be skipped
            const int skipKmer = max(lastPos - pos, 0);
            // construct corresponding sequence with reduced alphabet
            std::vector<char> redSeq(2*MyConst::READLEN - 2);
            std::vector<char> redSeqRev(2*MyConst::READLEN - 2);

            std::vector<char>::const_iterator start = seq.begin();
            std::vector<char>::const_iterator end = seq.begin();

            // last N before the CpG or last position that was looked at
            int lasN = skipKmer - 1;
            unsigned int j;

            // if we need to hash kmer starting at the left of current CpG
            if (lastPos < pos + MyConst::READLEN - MyConst::KMERLEN)
            {
                std::advance(start, pos + skipKmer);
                std::advance(end, pos + MyConst::READLEN);

                // position in substring
                j = skipKmer;
                // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
                for ( ; start != end; ++start, ++j)
                {

                    const int revPos = 2*MyConst::READLEN - 3 - j;
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
                            std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

                    }

                }
                // move to one after last N
                ++lasN;
                // reassign end to the position of G or one after last hashed kmer
                end = seq.begin();
                std::advance(end, pos + MyConst::READLEN - 1);
            } else {

                ++lasN;
                end = seq.begin();
                std::advance(end, pos + skipKmer - 1);
            }
            // move over second half in reverse order
            // reassign current position
            j = 2*MyConst::READLEN - 3;
            // reassign start to final position
            start = seq.begin();
            std::advance(start, pos + j);

            // offset where the first N after the CpG is
            int off = 2*MyConst::READLEN - 2;
            for ( ; start != end; --start, --j)
            {
                const int revPos = 2*MyConst::READLEN - 3 - j;
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
                        std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
                }


            }

            // length of the context around the CpG without Ns
            unsigned int contextLen = off - lasN;
            // if we don't have enough to read, return without hashing
            if (contextLen < MyConst::KMERLEN)
            {
                return;
            }

            const char* seqStart = redSeq.data() + lasN;
            const char* seqStartRev = redSeqRev.data() + (2*MyConst::READLEN - 2 - off);

            // initial hash backward
            uint64_t rhVal = ntHash::NTP64(seqStartRev);
            // first kmer on reverse complement corresponds to last kmer in forward sequence
            uint64_t kPosRev = off - MyConst::KMERLEN;

            // update kmer table
            kmerTable[--tabIndex[rhVal % tabIndex.size()]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
            strandTable[tabIndex[rhVal % tabIndex.size()]] = false;


            // hash kmers of backward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {
                --kPosRev;
                ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
                // update kmer table
                kmerTable[--tabIndex[rhVal % tabIndex.size()]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
                strandTable[tabIndex[rhVal % tabIndex.size()]] = false;
            }

            // initial hash forward
            uint64_t fhVal = ntHash::NTP64(seqStart);
            uint64_t kPos = lasN;

            // update kmer table
            kmerTable[--tabIndex[fhVal % tabIndex.size()]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
            strandTable[tabIndex[fhVal % tabIndex.size()]] = true;

            // hash kmers of forward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {
                ++kPos;
                ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
                // update kmer table
                kmerTable[--tabIndex[fhVal % tabIndex.size()]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
                strandTable[tabIndex[fhVal % tabIndex.size()]] = true;
            }
            lastPos = pos + off - MyConst::KMERLEN + 1;

        }
        void ntHashLast(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG, const uint32_t& metacpg, uint32_t&& metaOff);
        void ntHashFirst(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& cpgOffset, const uint32_t& metacpg);



        inline void ntCountChunk(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos)
        {

            // kmers to be skipped
            const int skipKmer = max(lastPos - pos, 0);
            // construct corresponding sequence with reduced alphabet
            std::vector<char> redSeq(2*MyConst::READLEN - 2);
            std::vector<char> redSeqRev(2*MyConst::READLEN - 2);

            std::vector<char>::const_iterator start = seq.begin();
            std::vector<char>::const_iterator end = seq.begin();

            // last N before the CpG or last position that was looked at
            int lasN = skipKmer - 1;
            unsigned int j;

            // if we need to hash kmer starting at the left of current CpG
            if (lastPos < pos + MyConst::READLEN - MyConst::KMERLEN)
            {
                std::advance(start, pos + skipKmer);
                std::advance(end, pos + MyConst::READLEN);

                // position in substring
                j = skipKmer;
                // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
                for ( ; start != end; ++start, ++j)
                {

                    const int revPos = 2*MyConst::READLEN - 3 - j;
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
                            std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

                    }

                }
                // reassign end to the position of G or one after last hashed kmer
                end = seq.begin();
                std::advance(end, pos + MyConst::READLEN - 1);
            } else {

                end = seq.begin();
                std::advance(end, pos + skipKmer - 1);
            }
            // move to one after last N/ one after last hashed kmer
            ++lasN;
            // move over second half in reverse order
            // reassign start to final position
            start = seq.begin();
            std::advance(start, pos + 2*MyConst::READLEN - 3);

            // reassign current position
            j = 2*MyConst::READLEN - 3;
            // offset where the first N after the CpG start is
            int off = 2*MyConst::READLEN - 2;
            for ( ; start != end; --start, --j)
            {
                const int revPos = 2*MyConst::READLEN - 3 - j;
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
                        std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
                }


            }

            // length of the context around the CpG that should be hashed
            unsigned int contextLen = off - lasN;
            // if we don't have enough to read, return without hashing
            if (contextLen < MyConst::KMERLEN)
            {
                return;
            }

            const char* seqStart = redSeq.data() + lasN;
            const char* seqStartRev = redSeqRev.data() + (2*MyConst::READLEN - 2 - off);

            // initial hash backward
            uint64_t rhVal = ntHash::NTP64(seqStartRev);

            // update indices
            ++tabIndex[rhVal % tabIndex.size()];

            // hash kmers of backward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {
                ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
                // update indices
                ++tabIndex[rhVal % tabIndex.size()];

            }

            // initial hash forward
            uint64_t fhVal = ntHash::NTP64(seqStart);

            // update indices
            ++tabIndex[fhVal % tabIndex.size()];

            // hash kmers of forward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {

                ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
                // update indices
                ++tabIndex[fhVal % tabIndex.size()];
            }
            // update position of last hashed kmer (+ 1)
            lastPos = pos + off - MyConst::KMERLEN + 1;

        }
        void ntCountLast(std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG);
        void ntCountFirst(std::vector<char>& seq, uint32_t& lastPos, const unsigned int& cpgOffset);

        // estimates the number of collision per entry and number of overall kmers to be hashed
        // to initialize tabIndex and kmerTable
        inline void estimateTablesizes(std::vector<std::vector<char> >& genomeSeq)
        {

            // count start CpG kmers
            for (metaCpG& m : metaStartCpGs)
            {

                // we know that all of the CpGs at start will overlap
                const uint8_t chr = cpgStartTable[m.start].chrom;

                uint32_t lastPos = 0;

                for (uint32_t cpgInd = m.start; cpgInd <= m.end; ++cpgInd)
                {
                    ntCountFirst(genomeSeq[chr], lastPos, cpgStartTable[cpgInd].pos);

                }
            }
            // count normal CpG kmers
            for (metaCpG& m : metaCpGs)
            {

                const uint8_t chr = cpgTable[m.start].chrom;

                uint32_t lastPos = 0;

                // how long is the rest of the sequence after the end of last cpg in meta cpg
                unsigned int remainderBps = genomeSeq[cpgTable[m.start].chrom].size() - (cpgTable[m.start].pos + MyConst::READLEN);

                // kmers of first CpG
                if (remainderBps >= (MyConst::READLEN - 2) )
                {
                    ntCountChunk(genomeSeq[chr], lastPos, cpgTable[m.start].pos);

                } else {

                    ntCountLast(genomeSeq[chr], lastPos, cpgTable[m.start].pos, remainderBps);
                }

                // consecutive CpG kmers
                for (uint32_t cpgInd = m.start + 1; cpgInd <= m.end; ++cpgInd)
                {

                    remainderBps = genomeSeq[cpgTable[cpgInd].chrom].size() - (cpgTable[cpgInd].pos + MyConst::READLEN);
                    // if we can read the full sequence breadth after the CpG
                    if (remainderBps >= (MyConst::READLEN - 2) )
                    {

                        // count the collisions
                        ntCountChunk(genomeSeq[chr], lastPos, cpgTable[cpgInd].pos);

                    } else {

                        // count the collisions and break
                        ntCountLast(genomeSeq[chr], lastPos, cpgTable[cpgInd].pos, remainderBps);
                        break;
                    }
                }


            }

            unsigned long sum = 0;
            // update to sums of previous entrys
            for (unsigned int i = 0; i < tabIndex.size(); ++i)
            {

                sum += tabIndex[i];
                tabIndex[i] = sum;
            }

            // resize table to number of kmers that have to be hashed
            kmerTable.resize(sum);
            strandTable.resize(sum);

        }





        // table of all CpGs in reference genome
        const std::vector<struct CpG> cpgTable;
        const std::vector<struct CpG> cpgStartTable;

        // table of bitstrings* holding a bit representation of genomeSeq
        // used as a perfect hash later on
        // encoding:
        //          A -> 00
        //          C -> 01
        //          T -> 11
        //          G -> 10
        //
        //          N -> 00   DISCARDED for all computations
        //
        // *own implementation of bitstrings
        std::vector<DnaBitStr> genomeBit;

        // hash table
        // tabIndex [i] points into kmerTable where the first entry with hash value i is saved
        // e.g. tabIndex[2] == 3 <=> hash(kmerTable[i]) < 3 forall i < 3
        // kmerTable holds the kmer (i.e. MetaCpg index and offset)
        // strandTable hold the strand orientation of the corresponding kmer
        std::vector<uint64_t> tabIndex;
        std::vector<KMER::kmer> kmerTable;
        std::vector<bool> strandTable;
        //
        // meta CpG table
        std::vector<struct metaCpG> metaCpGs;
        std::vector<struct metaCpG> metaStartCpGs;

};

#endif /* REFGENOME_H */
