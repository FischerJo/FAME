#ifndef CONST_H
#define CONST_H

#include <cstdint>

namespace MyConst {



//  -------- SET THESE VARIABLES --------
//
//
// number of chromosomes in organism
constexpr unsigned int CHROMNUM = 24;

// (more than) size in bp of biggest chromosome in organism
constexpr unsigned int CHROMMAX = 260000000;

// (more than) number of CpGs in organism
constexpr unsigned int CPGMAX = 24000000;

// maximum read length of the reads in bp
constexpr unsigned int READLEN = 30;



//  --------------------------------------




// --------------------
// ------INTERNAL------
// --------------------


// Length of a kmer in bp
constexpr unsigned int KMERLEN = 20;

// Kmer bitmask to build reverse complement of 64bit kmer
constexpr uint64_t KMERMASK = 0xffffffffffffffffULL;

// yields value for number of kmers in one read at compile time
// constexpr KPERREAD = READLEN - KMERLEN + 1;



// internal buffer for reading reference file
constexpr unsigned int BUFSIZE = 256 * 1024;


}

#endif /* CONST_H */
