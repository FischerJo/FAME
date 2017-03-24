#ifndef REFGENOME_H
#define REFGENOME_H

#include <fstream>
#include <istream>
#include <string>
#include <vector>

// Project includes
#include "CONST.h"
#include "structs.h"


// Class representing the reference genome
class RefGenome
{

    public:

        RefGenome() = delete;

        // Ctor
        // Arguments: ( see class members for full information )
        //      cpgs        number of CpGs in reference genome
        //      cpgTab      table of all CpGs in reference genome
        //      genSeq      genomic sequence seperated by chromosome
        RefGenome(unsigned int cpgs, std::vector<struct CpG>& cpgTab, std::vector<std::string>& genSeq);


        ~RefGenome();


        // hash all kmers in all CpGs to _kmerTable using ntHash
        // the kmers are represented in REDUCED alphabet {A,T,G}
        // mapping all Cs to Ts
        // generate Bit representation of whole genome for genomeBit
        // using full alphabet
        void generateHashes();


    private:


        // number of CpGs in reference genome
        const unsigned int cpgNum;

        // table of all CpGs in reference genome
        // _cpgTable has size _cpgNum
        const std::vector<struct CpG> cpgTable;

        // table of strings holding chromosome code
        // convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
        const std::vector<std::string> genomeSeq;

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
        // std::vector<const dnaBitStr> genomeBit;

        // hash table of all (reduced alphabet) kmers that surround CpGs in reference,
        // hash values are computed using nthash
        std::vector<struct kmer> kmerTable;


};

#endif /* REFGENOME_H */
