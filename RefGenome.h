#ifndef REFGENOME_H
#define REFGENOME_H

#include <fstream>
#include <istream>
#include <string>
#include <vector>

// Project includes
#include "CONST.h"
#include "structs.h"
#include "DnaBitStr.h"


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
        inline void ntHashChunk(const char* seq, const std::size_t seqLen)
        {

            // find last occurence of N before CpG and first occurence after
            char* lastBef;
            // TODO what if N is last?
            for (lastBef = seq; (lastBef = (char*) memchr(lastBef, 'N', MyConst::READLEN - 1)); ++lastBef)
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
            if (offset < KMERLEN)
            {
                return;
            }

            // how often should we apply the rolling
            // if zero, make just initial hash
            unsigned int rolls = offset - KMERLEN;
            uint64_t fhVal = 0;
            uint64_t rhVal = 0;
            ntHash::NTPC64(lastBef, fhVal, rhVal);
            // TODO make kmer and put it in table

            for (int i = 0; i < rolls; ++i)
            {

                ntHash::NTPC64(lastBef[i], lastBef[KMERLEN + i], fhVal, rhVal);
                // TODO make kmer and put it in table
            }


        }

        // table of all CpGs in reference genome
        // _cpgTable has size _cpgNum
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
