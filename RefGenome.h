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
        // Arguments: ( see class members for full information )
        //      cpgTab      table of all CpGs in reference genome except for the ones near the start of a sequence (i.e. less then READLEN away from start)
        //      cpgStartTab table of CpGs near the start
        //      genSeq      genomic sequence seperated by chromosome
        RefGenome(std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<std::vector<char> >& genSeq);

        ~RefGenome() = default;





    // TODO: make this private once its tested
    // private:

        // generate Bit representation of whole genome for genomeBit
        // using full alphabet
        void generateBitStrings(std::vector<std::vector<char> >& genomeSeq);


        // hash all kmers in all CpGs to _kmerTable using ntHash
        // the kmers are represented in REDUCED alphabet {A,T,G}
        // mapping all Cs to Ts
        void generateHashes(std::vector<std::vector<char> >& genomeSeq);


        // generates all kmers in seq and hashes them and their reverse complement using nthash into kmerTable
        inline void ntHashChunk(const std::vector<char>& seq, const unsigned int& pos, const uint32_t& cpg)
        {

            // construct corresponding sequence with reduced alphabet
            std::vector<char> redSeq(2*MyConst::READLEN - 2);
            std::vector<char> redSeqRev(2*MyConst::READLEN - 2);

            std::vector<char>::const_iterator start = seq.begin();
            std::advance(start, pos);
            std::vector<char>::const_iterator end = seq.begin();
            std::advance(end, pos + MyConst::READLEN);

            // position in substring
            unsigned int j = 0;
            // last N before the CpG
            int lasN = -1;
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
            // move over second half in reverse order
            // reassign end to the position of G
            --end;
            // reassign start to final position
            start = seq.begin();
            std::advance(start, pos + 2*MyConst::READLEN - 3);
            j = 2*MyConst::READLEN - 3;
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
            char* seqStart = redSeq.data() + lasN;
            char* seqStartRev = redSeqRev.data() + (2*MyConst::READLEN - 2 - off);

            // initial hash backward
            uint64_t fhVal = ntHash::NTP64(seqStartRev);
            uint16_t kPos = lasN;

            // update kmer table
            ++kmerTable[fhVal % kmerTable.size()].len;
            kmerTable[fhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);

            // initial hash forward
            uint64_t rhVal = ntHash::NTP64(seqStart);
            kPos |= 0x4000;

            // update kmer table
            ++kmerTable[rhVal % kmerTable.size()].len;
            kmerTable[rhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);



            for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
            {

                rhVal = ntHash::NTP64(rhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
                kPos = lasN + i + 1;
                // update kmer table
                ++kmerTable[rhVal % kmerTable.size()].len;
                kmerTable[rhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);

                fhVal = ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
                kPos |= 0x4000;
                // update kmer table
                ++kmerTable[fhVal % kmerTable.size()].len;
                kmerTable[fhVal % kmerTable.size()].collis.emplace_back(cpg, kPos);
            }
        }
        void ntHashLast(const std::vector<char>& seq, const unsigned int& pos, const unsigned int& bpsAfterCpG, const uint32_t& cpg);
        void ntHashFirst(const std::vector<char>& seq, const unsigned int& cpgOffset, const uint32_t& cpg);

        // table of all CpGs in reference genome
        const std::vector<struct CpG> cpgTable;
        const std::vector<struct CpG> cpgStartTable;

        // table of strings holding chromosome code
        // and table holding the length of the sequences
        // convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
        // const std::vector<std::vector<char> > genomeSeq;
        // const std::vector<std::size_t> genomeSeqLen;

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

        // hash table of all (reduced alphabet) kmers that surround CpGs in reference,
        // hash values are computed using nthash
        std::vector<colList > kmerTable;


};

#endif /* REFGENOME_H */
