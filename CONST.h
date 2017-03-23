#ifndef CONST_H
#define CONST_H

// Length of the reads in bp
constexpr unsigned int READLEN = 100;

// Length of a kmer in bp
constexpr unsigned int KMERLEN = 20;

// yields value for number of kmers in one read at compile time
constexpr kPerRead = READLEN - KMERLEN;

// INTERNAL


// internal buffer for reading reference file
constexpr unsigned int BUFSIZE = 256 * 1024;


#endif /* CONST_H */
