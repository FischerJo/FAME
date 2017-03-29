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


        // hash all kmers in all CpGs to _kmerTable using ntHash
        // the kmers are represented in REDUCED alphabet {A,T,G}
        // mapping all Cs to Ts
        // generate Bit representation of whole genome for genomeBit
        // using full alphabet
        void generateHashes();


        // generate the Bitstrings (i.e. sequence encoding and masks) for sequence
        void generateBitStr();


    private:

        // generates all kmers in seq and hashes them and their reverse complement using nthash into kmerTable
        // CONVENTION:
        inline void ntHashChunk(const std::vector<char>& seq, const unsigned int& pos, const struct CpG& cpg)
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
            for ( ; off < (MyConst::READLEN - 2); ++off)
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
        void ntHashLast(const std::vector<char>& seq, const unsigned int& pos, const unsigned int& bpsAfterCpG, const struct CpG& cpg);
        void ntHashFirst(const std::vector<char>& seq, const unsigned int& cpgOffset, const struct CpG& cpg);

        // table of all CpGs in reference genome
        const std::vector<struct CpG> cpgTable;
        const std::vector<struct CpG> cpgStartTable;

        // table of strings holding chromosome code
        // and table holding the length of the sequences
        // convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
        const std::vector<std::vector<char> > genomeSeq;
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
        std::vector<std::list<struct kmer> > kmerTable;


};

#endif /* REFGENOME_H */
