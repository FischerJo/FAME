#ifndef READ_H
#define READ_H

#include <string>
#include <fstream>
#include <vector>

#include "structs.h"

class Read
{

    public:

        Read(std::string& seq, std::string& id);

        // print the matching info of this read to ofs
        // FORMAT:
        //      Chromosome\tPosition\tKmerOverlap\tReadSeq\n*
        //      , where
        //      Chromosome is the index of the chromosome
        //      Position is the 0 based offset in the chromosome
        //      KmerOverlap is the number of kmers that match the position
        //          (up to ReadSize - MyConst::KMERLEN + 1)
        //      ReadSeq is the sequence of the read
        //
        //      done for each match
        //
        void printMatch(std::ofstream& ofs);


        // print the matching info to ofs in SAM format
        void printMatchSam(std::ofstream& ofs);

        std::string id;

        // (DNA) sequence of the read
        std::string seq;

        struct match {
            struct CpG cpg;
            unsigned int offset;
        };
        // The matching positions of the read
        std::vector<struct match> matches;

        // flag stating if read has an N in sequence or is too small (< Kmerlength)
        // if this is the case, it won't be processed
        bool isInvalid;

};

#endif /* READ_H */
