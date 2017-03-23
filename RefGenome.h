#ifndef REFGENOME_H
#define REFGENOME_H

#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <cstdint>

// Project includes
#include "CONST.h"

// yields value for number of kmers in one read at compile time
constexpr kPerRead(unsigned int k, unsigned int rL) {rL - k}


struct CpG {
    // position where this CpG is present in genome
    // convention: pos points to the location of the C
    const uint8_t chrom;
    const unsigned int pos;
};


struct kmer {
    // pointer to related CpG
    const struct CpG * cpg;
    // offset in this CpG, the C of the CpG is offset 0
    // 2nd highest bit is flag for strand; use STRANDMASK to extract information
    //      1 -> forward strand
    //      0 -> reverse strand
    const uint16_t offset;

};

// Use this mask to extract strand information from kmer.offset
// offset & STRANDMASK == TRUE   <=>   kmer on forward strand
constexpr STRANDMASK = 0x4000;

// Class representing the reference genome
class RefGenome
{

    public:

        delete RefGenome();

        // read file in FASTA format
        // Throws out unlocated reads
        RefGenome(std::ifstream&);


        ~RefGenome();


        // hash all kmers in all CpGs to _kmerTable using ntHash
        // the kmers are represented in REDUCED alphabet {A,T,G}
        // mapping all Cs to Ts
        void hashCpgKmers();


    private:


        // number of CpGs in reference genome
        const unsigned int _cpgNum;

        // table of all CpGs in reference genome
        // _cpgTable has size _cpgNum
        const std::vector<const struct CpG> _cpgTable;

        // table of strings holding chromosome code
        // convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
        const std::vector<const std::string> _genomeSeq;

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
        const std::vector<const dnaBitStr> _genomeBit;

        // hash table of all (reduced alphabet) kmers that surround CpGs in reference,
        // hash values are computed using nthash
        std::vector<const struct kmer> _kmerTable;


};

#endif /* REFGENOME_H */
