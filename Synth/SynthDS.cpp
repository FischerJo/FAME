
#include <omp.h>
#include <iostream>

#include "SynthDS.h"

SynthDS::SynthDS(const size_t refLen) :
        toIndex(0,3)
    ,   alphabet {{ 'A', 'C', 'G', 'T' }}
{
    std::random_device seedGen;

    // init random number generator for each thread
    for (unsigned int cID = 0; cID < CORENUM; ++cID)
    {
        randGen[cID] = std::mt19937(seedGen());
    }
    initReference(refLen, seedGen());
    std::cout << "Gerated reference sequence\n\n";
}

SynthDS::SynthDS(const size_t refLen, const unsigned int seed) :
        toIndex(0,3)
    ,   alphabet {{ 'A', 'C', 'G', 'T' }}
{
    std::random_device seedGen;

    // init random number generator for each thread
    for (unsigned int cID = 0; cID < CORENUM; ++cID)
    {
        randGen[cID] = std::mt19937(seedGen());
    }
    initReference(refLen, seed);
    std::cout << "Gerated reference sequence\n\n";
}

std::vector<std::string> SynthDS::genReadsFwdFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<size_t>& offsets)
{

    // will hold the generated reads
    std::vector<std::string> readSet(readNum);
    offsets.resize(readNum);

    // range of indices allowed to be drawn for the reference
    std::uniform_int_distribution<int> toOffset (0, refFwd.size() - readLen);

    // range of indices allowed to be drawn for errors in read
    std::uniform_int_distribution<int> toOffRead (0, readLen - 1);

    // range of errors allowed to be drawn
    std::uniform_int_distribution<int> toErr(0, maxErrNum);

    // coin flip for conversion type
    std::uniform_int_distribution<int> coin(0, 1);

#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
    for (size_t i = 0; i < readNum; ++i)
    {
        std::mt19937& MT = randGen[omp_get_thread_num()];
        // draw an offset where to draw the read
        const size_t offset = toOffset(MT);
        // generate sequence
        std::string read(refFwd, offset, readLen);

        // draw number of errors
        const unsigned int err = toErr(MT);
        // introduce errors at random positions
        for (unsigned int e = 0; e < err; ++e)
        {
            read[toOffRead(MT)] = alphabet[toIndex(MT)];

        }
        // introduce C->T  OR  G->A conversions
        if (coin(MT))
        {

            for (size_t j = 0; j < readLen; ++j)
            {
                if (read[j] == 'C')
                {
                    // try conversion
                    if (coin(MT))
                    {
                        read[j] = 'T';
                    }

                }
            }

        } else {

            for (size_t j = 0; j < readLen; ++j)
            {
                if (read[j] == 'G')
                {
                    // try conversion
                    if (coin(MT))
                    {
                        read[j] = 'A';
                    }

                }
            }
        }

        readSet[i] = std::move(read);
        offsets[i] = offset;
    }

    std::cout << "Generated forward strand read set\n\n";
    return readSet;
}

std::vector<std::string> SynthDS::genReadsRevFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum)
{

    // will hold the generated reads
    std::vector<std::string> readSet(readNum);

    // range of indices allowed to be drawn for the reference
    std::uniform_int_distribution<int> toOffset (0, refRev.size() - readLen);

    // range of indices allowed to be drawn for errors in read
    std::uniform_int_distribution<int> toOffRead (0, readLen - 1);

    // range of errors allowed to be drawn
    std::uniform_int_distribution<int> toErr(0, maxErrNum);

    // coin flip for conversion type
    std::uniform_int_distribution<int> coin(0, 1);

#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
    for (size_t i = 0; i < readNum; ++i)
    {
        std::mt19937& MT = randGen[omp_get_thread_num()];
        // draw an offset where to draw the read
        const size_t offset = toOffset(MT);
        // generate sequence
        std::string read(refRev, offset, readLen);

        // draw number of errors
        const unsigned int err = toErr(MT);
        // introduce errors at random positions
        for (unsigned int e = 0; e < err; ++e)
        {
            read[toOffRead(MT)] = alphabet[toIndex(MT)];

        }
        // introduce C->T  OR  G->A conversions
        if (coin(MT))
        {

            for (size_t j = 0; j < readLen; ++j)
            {
                if (read[j] == 'C')
                {
                    // try conversion
                    if (coin(MT))
                    {
                        read[j] = 'T';
                    }

                }
            }

        } else {

            for (size_t j = 0; j < readLen; ++j)
            {
                if (read[j] == 'G')
                {
                    // try conversion
                    if (coin(MT))
                    {
                        read[j] = 'A';
                    }

                }
            }
        }
        readSet[i] = std::move(read);
    }

    std::cout << "Generated reverse complement strand read set\n\n";
    return readSet;
}

// std::vector<std::string> SynthDS::genReadsFwdRange(const size_t readLenMin, const size_t readLenMax, const size_t readNum, const unsigned int maxErrNum)
// {
//
// }

void SynthDS::initReference(const size_t refLen, const unsigned int seed)
{

    refFwd = std::string(refLen, 'N');
    refRev = std::string(refLen, 'N');
    std::mt19937 MT(seed);

    for (size_t i = 0; i < refLen; ++i)
    {

        // draw rnd number
        const unsigned int r = toIndex(MT);
        // generate random letter from the alphabet and append to reference sequence
        refFwd[i] = alphabet[r];
        refRev[refLen - 1 - i] = alphabet[3-r];
    }
}
