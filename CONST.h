#ifndef CONST_H
#define CONST_H

#include <cstdint>
#include <array>

namespace MyConst {



//  -------- SET THESE VARIABLES --------
//
//
// number of chromosomes in organism
constexpr unsigned int CHROMNUM = 24;

// (more than) size in bp of biggest chromosome in organism
constexpr unsigned int CHROMMAX = 1000000000;

// (more than) number of CpGs in organism
constexpr unsigned int CPGMAX = 100000000;

// maximum read length of the reads in bp
constexpr unsigned int READLEN = 110;

// number of reads that should be read per batch
// (note that very large values may increase the required amount of RAM drastically)
constexpr unsigned int CHUNKSIZE = 100000;

// maximum number of times a k-mer is allowed to occur in the whole genome
constexpr uint64_t KMERCUTOFF = 1000;

// Number of cores that this program is allowed to occupy at any given point
#define CORENUM 1


//  --------------------------------------




// --------------------
// ------INTERNAL------
// --------------------


// Length of a kmer in bp
// recommended is 30
constexpr unsigned int KMERLEN = 30;

// Kmer bitmask to build reverse complement of 64bit kmer
constexpr uint64_t KMERMASK = (KMERLEN == 32 ? 0xffffffffffffffffULL : (1ULL << 2*KMERLEN) - 1);

// yields value for number of kmers in one read at compile time
// constexpr KPERREAD = READLEN - KMERLEN + 1;


// size of hash table
// recommended is 1 << 30
constexpr uint64_t HTABSIZE = 1ULL << 30;


// window length for meta CpGs
constexpr unsigned int WINLEN = 1 << 10;

// number of mismatches we allow
constexpr unsigned int MISCOUNT = 1;

// threshold over which reads are immediately discarded if that many metaCpGs are hit
constexpr unsigned int QUERYTHRESHOLD = 30000;

// Checks the usefulness of the set parameters
void sanityChecks();


}

#endif /* CONST_H */
