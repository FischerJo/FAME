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
        //      cpgTab      table of all CpGs in reference genome
        //      genSeq      genomic sequence seperated by chromosome
        RefGenome(std::vector<struct CpG>& cpgTab, std::vector<const char*>& genSeq, std::vector<std::size_t> genSeqLen);

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
        //              seq points to a character sequence containing at least 2*READLEN - 2 many chars
        inline void ntHashChunk(const char* seq, const struct CpG& cpg)
        {

            // find last occurence of N before CpG and first occurence after
            const char* lastBef;
            for (lastBef = seq; (lastBef = (char*) memchr(lastBef, 'N', (MyConst::READLEN - 1 - (lastBef - seq)))); ++lastBef)
            {}
            // now lastBef points to character one after the last N before the CpG

            char* firstAft = (char*) memchr(seq + MyConst::READLEN, 'N', MyConst::READLEN - 2);

            // No N found before CpG, set char* to start of argument sequence
            if (lastBef == NULL)
            {
                lastBef = seq;
            }
            // how many chars to consider, init as full sequence around CpG
            std::size_t offset = 2*MyConst::READLEN - 2 - (lastBef - seq);
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
        void ntHashLast(const char* seq, const unsigned int bpsAfterCpG, const struct CpG& cpg);
        void ntHashFirst(const char* seq, const unsigned int posOfCpG, const struct CpG& cpg);

        // table of all CpGs in reference genome
        const std::vector<struct CpG> cpgTable;

        // table of strings holding chromosome code
        // and table holding the length of the sequences
        // convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
        const std::vector<const char*> genomeSeq;
        const std::vector<std::size_t> genomeSeqLen;

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
