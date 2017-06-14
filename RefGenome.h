#ifndef REFGENOME_H
#define REFGENOME_H

#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <unordered_map>

// Project includes
#include "CONST.h"
#include "structs.h"
#include "DnaBitStr.h"
#include "ntHash-1.0.2/nthash.hpp"


// TODO: Maybe save reduced sequence also, to speed up the hashing process?
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


            // retrieve kmers for first hash
            uint64_t hVal = ntHash::NTP64(seq.data());
            uint64_t hashInd = hVal % MyConst::HTABSIZE;

            // iterators for kmertable
            std::vector<KMER::kmer>::iterator startK = kmerTable.begin() + tabIndex[hashInd];
            std::vector<KMER::kmer>::iterator endK = kmerTable.begin() + tabIndex[hashInd + 1];
            // iterators for strandtable
            std::vector<bool>::iterator startS = strandTable.begin() + tabIndex[hashInd];
            std::vector<bool>::iterator endS = strandTable.begin() + tabIndex[hashInd + 1];

            // retrieve seeds
            seedsK.emplace_back(startK, endK);
            seedsS.emplace_back(startS, endS);

            for (unsigned int i = 0; i < (seq.size() - MyConst::KMERLEN); ++i)
            {

                // use rolling hash
                ntHash::NTP64(hVal, seq[i], seq[i + MyConst::KMERLEN]);
                hashInd = hVal % MyConst::HTABSIZE;

                startK = kmerTable.begin() + tabIndex[hashInd];
                endK = kmerTable.begin() + tabIndex[hashInd + 1];

                startS = strandTable.begin() + tabIndex[hashInd];
                endS = strandTable.begin() + tabIndex[hashInd + 1];

                // retrieve seeds
                seedsK.emplace_back(startK, endK);
                seedsS.emplace_back(startS, endS);
            }

        }

        // yields the internal bit representation of the reference genome at the specified chromosome
        // inline const DnaBitStr getBitRepr(unsigned int chr)
        // {
        //
        //     return genomeBit[chr];
        // }



    // TODO: make this private once its tested
    // private:

        // simple maximum inline function
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
            kmerTable[--tabIndex[rhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
            strandTable[tabIndex[rhVal % MyConst::HTABSIZE]] = false;


            // hash kmers of backward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {
                --kPosRev;
                ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
                // update kmer table
                kmerTable[--tabIndex[rhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
                strandTable[tabIndex[rhVal % MyConst::HTABSIZE]] = false;
            }

            // initial hash forward
            uint64_t fhVal = ntHash::NTP64(seqStart);
            uint64_t kPos = lasN;

            // update kmer table
            kmerTable[--tabIndex[fhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
            strandTable[tabIndex[fhVal % MyConst::HTABSIZE]] = true;

            // hash kmers of forward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {
                ++kPos;
                ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
                // update kmer table
                kmerTable[--tabIndex[fhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
                strandTable[tabIndex[fhVal % MyConst::HTABSIZE]] = true;
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
            ++tabIndex[rhVal % MyConst::HTABSIZE];

            // hash kmers of backward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {
                ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
                // update indices
                ++tabIndex[rhVal % MyConst::HTABSIZE];

            }

            // initial hash forward
            uint64_t fhVal = ntHash::NTP64(seqStart);

            // update indices
            ++tabIndex[fhVal % MyConst::HTABSIZE];

            // hash kmers of forward strand
            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {

                ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
                // update indices
                ++tabIndex[fhVal % MyConst::HTABSIZE];
            }
            // update position of last hashed kmer (+ 1)
            lastPos = pos + off - MyConst::KMERLEN + 1;

        }
        void ntCountLast(std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG);
        void ntCountFirst(std::vector<char>& seq, uint32_t& lastPos, const unsigned int& cpgOffset);

        // estimates the number of collision per entry and number of overall kmers to be hashed
        // to initialize tabIndex and kmerTable
        void estimateTablesizes(std::vector<std::vector<char> >& genomeSeq);


        // blacklist all k-mers that appear more then KMERCUTOFF times in the specified kmerTable slice
        // returns a list of k-mer sequences (as bitstrings) that are blacklisted through argument
        //
        // ARGUMENTS:
        //              KSliceStart Start index where to start blacklisting in KmerTable
        //              KSliceEnd   End index where to stop blacklisting (excluding this index)
        //
        //              blt         map of k-mer strings that are blacklisted
        //                          THIS WILL BE FILLED DURING CALL
        inline void blacklist(const unsigned int& KSliceStart, const unsigned int& KSliceEnd, std::unordered_map<uint64_t, unsigned int>& bl)
        {

            // iterate over kmerTable in the specified range
            for (unsigned int i = KSliceStart; i < KSliceEnd; ++i)
            {

                // get the kmer at that position
                KMER::kmer& k = kmerTable[i];
                bool sFlag = strandTable[i];

                const uint64_t kHash = reproduceKmerSeq(k, sFlag);

                // try to put sequence in map - if already in, count up
                // NOTE:    hash function is perfect, hence no implicit collisions before putting it into hashmap
                auto insertion = bl.find(kHash);
                if (insertion == bl.end())
                {

                    bl[kHash] = 1;

                } else {

                    insertion->second += 1;
                }
            }
        }

        // reproduce the k-mer sequence of a given kmer by looking up the position in the reference genome
        // the sequence will be returned as a bitstring
        //
        // ARGUMENTS:
        //              k       k-mer
        //              sFlag   flag stating if kmer is of forward (true) or reverse (false) strand
        //              kSeq    empty fixed length array to be filled with the kmer sequence
        inline uint64_t reproduceKmerSeq(KMER::kmer& k, bool sFlag)
        {

            uint32_t pos = 0;
            uint8_t chrom;
            // reproduce kmer sequence
            if (KMER::isStartCpG(k))
            {

                // get the start position and chromosome of surrounding meta cpg wrapped into a CpG
                const struct CpG& startCpG = cpgStartTable[metaStartCpGs[KMER::getMetaCpG(k)].start];
                chrom = startCpG.chrom;

            } else {

                // get the start position and chromosome of surrounding meta cpg wrapped into a CpG
                const struct CpG& startCpG = cpgTable[metaCpGs[KMER::getMetaCpG(k)].start];
                pos = startCpG.pos;
                chrom = startCpG.chrom;

            }

            // retrieve sequence
            //
            // get iterator to sequence start
            auto seqStart = fullSeq[chrom].begin() + KMER::getOffset(k) + pos;
            uint64_t kSeq = 0;

            // if from forward strand, just translate
            if (sFlag)
            {

                for (unsigned int i = 0; i < MyConst::KMERLEN; ++i)
                {
                    kSeq = kSeq << 2;

                    switch (seqStart[i])
                    {

                        case 'T':

                            kSeq += 3;
                            break;

                        case 'C':

                            kSeq += 1;
                            break;

                        case 'G':

                            kSeq += 2;
                            break;

                    // end switch
                    }

                // end forloop over sequence
                }

            // otherwise, build reverse complement
            } else {

                for (unsigned int i = MyConst::KMERLEN; i > 0; --i)
                {
                    kSeq = kSeq << 2;

                    switch (seqStart[i - 1])
                    {

                        case 'A':

                            kSeq += 3;
                            break;

                        case 'G':

                            kSeq += 1;
                            break;

                        case 'C':

                            kSeq += 2;
                            break;

                    // end switch
                    }

                // end forloop over sequence
                }
            }

            return kSeq;

        }

        // filter the current hash table according to the given blacklist
        // overwrites the internal kmerTable and strandTable structure, as well as tabIndex
        //
        void filterHashTable();



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
        // full sequence
        std::vector<std::vector<char> > fullSeq;

        // hash table
        // tabIndex [i] points into kmerTable where the first entry with hash value i is saved
        // kmerTable holds the kmer (i.e. MetaCpg index and offset)
        // strandTable hold the strand orientation of the corresponding kmer (true iff forward)
        std::vector<uint64_t> tabIndex;
        std::vector<KMER::kmer> kmerTable;
        std::vector<bool> strandTable;
        //
        // meta CpG table
        std::vector<struct metaCpG> metaCpGs;
        std::vector<struct metaCpG> metaStartCpGs;

};

#endif /* REFGENOME_H */
