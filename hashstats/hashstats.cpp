
#include <random>
#include <fstream>
#include <iostream>

#include "hashstats.h"

#include "../ntHash-1.0.2/nthash.hpp"

int main(int argc, char** argv)
{

    std::vector<unsigned int> hashTable (TABSIZE, 0);
    drawRandHashes(hashTable);
    writeStats(argv[1], hashTable);
    hashTable = std::vector<unsigned int>(TABSIZE, 0);
    drawDistHashes(hashTable);
    writeStats(argv[2], hashTable);

    return 0;
}

void drawRandHashes(std::vector<unsigned int>& hashTable)
{
    std::random_device seedGen;
    std::mt19937 MT(seedGen());
    std::uniform_int_distribution<int> toIndex(0, 2);

    char kmer [MyConst::KMERLEN];
    // sequences that should be constructed
    for (unsigned long i = 1; i <= SAMPSIZE; ++i)
    {

        // draw the individual positions of the kmer
        for (unsigned int k = 0; k < MyConst::KMERLEN; ++k)
        {

            kmer[k] = LETCODE[toIndex(MT)];

        }

        ++hashTable[ntHash::NTP64(kmer) % TABSIZE];

    }
}

void drawDistHashes(std::vector<unsigned int>& hashTable)
{

    char kmer [MyConst::KMERLEN];
    genDistSequence(kmer, 0, hashTable);
}

void genDistSequence(char* kmer, unsigned int pos, std::vector<unsigned int>& hashTable)
{
    // base case
    if (pos == MyConst::KMERLEN - 1)
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            kmer[pos] = LETCODE[i];
            ++hashTable[ntHash::NTP64(kmer) % TABSIZE];

        }
        return;
    }

    // recurse on pos
    for (unsigned int i = 0; i < 3; ++i)
    {
        kmer[pos] = LETCODE[i];
        genDistSequence(kmer, pos + 1, hashTable);

    }
}

void writeStats(char* filename, std::vector<unsigned int>& hashTable)
{

    std::ofstream ofs (filename);

    for (unsigned long i = 0; i < TABSIZE; ++i)
    {

        ofs << i << '\t' << hashTable[i] << '\n';
    }

    ofs.close();
}
