#ifndef HASHSTATS_H
#define HASHSTATS_H

#include <vector>
#include "../CONST.h"

// how many kmers should be drawn
// pow(3,20) same size as the distinct sequences for KMERLEN == 20
constexpr unsigned long SAMPSIZE = 3486784401;
// constexpr unsigned long SAMPSIZE = 15000000;

// size of the hash table
// 20 million
constexpr unsigned int TABSIZE = 20000000;
// constexpr unsigned int TABSIZE = 200000;

// encoding of numbers
constexpr char LETCODE[3] = {'C','G','T'};

// draw SAMPSIZE manz hashes, draw letters uniformly at random
// encode the kmers using nthash
//
// hashTable should be initialized with TABSIZE capacity
// will contain the number of kmers hashed to an element of a vector
void drawRandHashes(std::vector<unsigned int>& hashTable);

// generate all possible (distinct) kmers and generate corresponding hashes
//
// hashTable should be initialized with TABSIZE capacity
// will contain the number of kmers hashed to an element of a vector
void drawDistHashes(std::vector<unsigned int>& hashTable);

// recursive function generating all possible kmers, call with pos == 0
// will recurse until inclusively KMERLEN - 1 (base case)
// will cout up corresponding value in hashtable
void genDistSequence(char* kmer, unsigned int pos, std::vector<unsigned int>& hashTable);

// write out statistic
void writeStats(char* filename, std::vector<unsigned int>& hashTable);


#endif /* HASHSTATS_H */
