
#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>  // reverse

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
    std::cout << "Generated reference sequence\n\n";
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
    std::cout << "Generated reference sequence\n\n";
}

SynthDS::SynthDS(const char* genFile, const double methRate, const double convRate) :
        toIndex(0,3)
    ,   methToss(methRate)
    ,   convToss(convRate)
    ,   alphabet {{ 'A', 'C', 'G', 'T' }}
{

    std::random_device seedGen;

    // init random number generator for each thread
    for (unsigned int cID = 0; cID < CORENUM; ++cID)
    {
        randGen[cID] = std::mt19937(seedGen());
    }
    loadRefSeq(genFile);
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

std::vector<std::string> SynthDS::genReadsFwdRef(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<std::pair<size_t, size_t> >& offsets, std::vector<std::array<int, errNum> >& errOffs)
{

    std::cout << "Start generating forward strand read set\n\n";
    // will hold the generated reads
    std::vector<std::string> readSet(readNum);
    offsets.resize(readNum);
    errOffs.resize(readNum);

    size_t processedReads = 0;

    // size of the genome
    size_t genomeSize = 0;
    for (size_t chr = 0; chr < refSeqFwd.size(); ++chr)
        genomeSize += refSeqFwd[chr].size();

    // go over all chromosomes
    for (size_t chr = 0; chr < refSeqFwd.size(); ++chr)
    {

        // compute proportion of reads that should be generated by this chromosome
        const size_t chrReadNum = ((double)refSeqFwd[chr].size() / (double)genomeSize) * readNum;

        // range of indices allowed to be drawn for the reference
        std::uniform_int_distribution<int> toOffset (0, refSeqFwd[chr].size() - readLen);

        // range of indices allowed to be drawn for errors in read
        std::uniform_int_distribution<int> toOffRead (0, readLen - 1);

        // range of errors allowed to be drawn
        std::uniform_int_distribution<int> toErr(0, maxErrNum);

        // coin flip for conversion type (G to A or C to T)
        // i.e. forward or reverse of PCR of genomic fragment treated with bisulfite
        std::uniform_int_distribution<int> coin(0, 1);

#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
        for (size_t i = processedReads; i < ((chr == refSeqFwd.size() - 1) ? readNum : processedReads + chrReadNum); ++i)
        {

            std::mt19937& MT = randGen[omp_get_thread_num()];
            // draw an offset where to draw the read
            const size_t offset = toOffset(MT);
            // generate sequence
            std::string read(refSeqFwd[chr], offset, readLen);
            // test if CpG is contained
            bool prevC = false;
            bool hasCpG = false;
            // search for CpGs
            for (const char c : read)
            {
                if (c == 'C')
                {
                    prevC = true;
                    continue;
                }
                if (c == 'G')
                {
                    if (prevC)
                        hasCpG = true;
                }
                prevC = false;
            }

            // produce only reads with CpG
            if (!hasCpG)
            {
                --i;
                continue;
            }

            // draw number of errors
            const unsigned int err = toErr(MT);
            // introduce errors at random positions
            size_t errID = 0;
            for (unsigned int e = 0; e < err; ++e)
            {
                int eOff = toOffRead(MT);
                read[eOff] = alphabet[toIndex(MT)];
                errOffs[i][errID] = eOff;
                ++errID;
            }
            for (; errID < errNum; ++errID)
            {
                errOffs[i][errID] = -1;
            }
            // introduce C->T  OR  G->A conversions
            if (coin(MT))
            {

                for (size_t j = 0; j < readLen; ++j)
                {
                    if (read[j] == 'C')
                    {
                        // try conversion
                        if (!methToss(MT))
                        {
                            if (convToss(MT))
                                read[j] = 'T';
                            // if necessary, count up methylation counters
                            if (j < readLen - 1)
                            {
                                if (read[j+1] == 'G')
                                {
                                    cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].first += 1;
                                }
                            }
                        } else {
                            // if necessary, count up methylation counters
                            if (j < readLen - 1)
                            {
                                if (read[j+1] == 'G')
                                {
                                    cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].second += 1;
                                }
                            }
                        }

                    }
                }

            } else {

                for (size_t j = 0; j < readLen; ++j)
                {
                    if (read[j] == 'G')
                    {
                        // try conversion
                        if (!methToss(MT))
                        {
                            if (convToss(MT))
                                read[j] = 'A';
                            // if necessary, count up methylation counters
                            if (j > 0)
                            {
                                if (read[j-1] == 'C')
                                {
                                    cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].first += 1;
                                }
                            }
                        } else {

                            // if necessary, count up methylation counters
                            if (j > 0)
                            {
                                if (read[j-1] == 'C')
                                {
                                    cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].second += 1;
                                }
                            }
                        }
                    }
                }
            }

            readSet[i] = std::move(read);
            offsets[i] = std::move(std::pair<size_t,size_t>(offset, chr));
        }
        processedReads += chrReadNum;

    }
    std::cout << "Generated forward strand read set\n\n";
    return readSet;
}

std::vector<std::string> SynthDS::genReadsRevRef(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<std::pair<size_t, size_t> >& offsets, std::vector<std::array<int, errNum> >& errOffs)
{

    std::cout << "Start generating reverse strand read set\n\n";
    // will hold the generated reads
    std::vector<std::string> readSet(readNum);
    offsets.resize(readNum);
    errOffs.resize(readNum);

    size_t processedReads = 0;

    // size of the genome
    size_t genomeSize = 0;
    for (size_t chr = 0; chr < refSeqRev.size(); ++chr)
        genomeSize += refSeqRev[chr].size();

    // go over all chromosomes
    for (size_t chr = 0; chr < refSeqRev.size(); ++chr)
    {

        // compute proportion of reads that should be generated by this chromosome
        const size_t chrReadNum = ((double)refSeqRev[chr].size() / (double)genomeSize) * readNum;

        // range of indices allowed to be drawn for the reference
        std::uniform_int_distribution<int> toOffset (0, refSeqRev[chr].size() - readLen);

        // range of indices allowed to be drawn for errors in read
        std::uniform_int_distribution<int> toOffRead (0, readLen - 1);

        // range of errors allowed to be drawn
        std::uniform_int_distribution<int> toErr(0, maxErrNum);

        // coin flip for conversion type (G to A or C to T)
        // i.e. forward or reverse of PCR of genomic fragment treated with bisulfite
        std::uniform_int_distribution<int> coin(0, 1);

#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
        for (size_t i = processedReads; i < ((chr == refSeqRev.size() - 1) ? readNum : processedReads + chrReadNum); ++i)
        {

            std::mt19937& MT = randGen[omp_get_thread_num()];
            // draw an offset where to draw the read
            const size_t offset = toOffset(MT);
            // generate sequence
            std::string read(refSeqRev[chr], offset, readLen);
            // test if CpG is contained
            bool prevC = false;
            bool hasCpG = false;
            // search for CpGs
            for (const char c : read)
            {
                if (c == 'C')
                {
                    prevC = true;
                    continue;
                }
                if (c == 'G')
                {
                    if (prevC)
                        hasCpG = true;
                }
                prevC = false;
            }

            // produce only reads with CpG
            if (!hasCpG)
            {
                --i;
                continue;
            }

            // draw number of errors
            const unsigned int err = toErr(MT);
            // introduce errors at random positions
            size_t errID = 0;
            for (unsigned int e = 0; e < err; ++e)
            {
                int eOff = toOffRead(MT);
                read[eOff] = alphabet[toIndex(MT)];
                errOffs[i][errID] = eOff;
                ++errID;
            }
            for (; errID < errNum; ++errID)
            {
                errOffs[i][errID] = -1;
            }
            // introduce C->T  OR  G->A conversions
            if (coin(MT))
            {

                for (size_t j = 0; j < readLen; ++j)
                {
                    if (read[j] == 'C')
                    {
                        // try conversion
                        if (!methToss(MT))
                        {
                            if (convToss(MT))
                                read[j] = 'T';
                            // if necessary, count up methylation counters
                            if (j < readLen - 1)
                            {
                                if (read[j+1] == 'G')
                                {
                                    cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].first += 1;
                                }
                            }
                        } else {
                            // if necessary, count up methylation counters
                            if (j < readLen - 1)
                            {
                                if (read[j+1] == 'G')
                                {
                                    cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].second += 1;
                                }
                            }
                        }
                    }
                }

            } else {

                for (size_t j = 0; j < readLen; ++j)
                {
                    if (read[j] == 'G')
                    {
                        // try conversion
                        if (!methToss(MT))
                        {
                            if (convToss(MT))
                                read[j] = 'A';
                            // if necessary, count up methylation counters
                            if (j > 0)
                            {
                                if (read[j-1] == 'C')
                                {
                                    cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].first += 1;
                                }
                            }

                        } else {

                            // if necessary, count up methylation counters
                            if (j > 0)
                            {
                                if (read[j-1] == 'C')
                                {
                                    cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ].second += 1;
                                }
                            }

                        }
                    }
                }
            }

            readSet[i] = std::move(read);
            offsets[i] = std::move(std::pair<size_t,size_t>(refSeqRev[chr].size() - (offset + 100 - 1), chr));
        }
        processedReads += chrReadNum;

    }
    std::cout << "Generated reverse strand read set\n\n";
    return readSet;
}

void SynthDS::loadRefSeq(const char* genFile)
{

    constexpr size_t chromNum = 24;
    constexpr size_t chromLen = 450000000;

    std::string line;

    // stores the sequence of each chromosome
    std::string seqFwd;
    std::string seqRev;

    // default init to human ref
    refSeqFwd.reserve(chromNum);
    refSeqRev.reserve(chromNum);
    seqFwd.reserve(chromLen);
    seqRev.reserve(chromLen);

    // stores the chromosome index we are currently reading
    uint8_t chrIndex = 0;

    // flag stating if we currently read a real chromosome assembly
    bool contFlag = false;
    // flag stating that last read character was a C
    bool lastC = false;

    std::ifstream ifs (genFile);

    std::cout << "Start reading reference file " << genFile << "\n";

    while (getline(ifs, line)) {

        // Test for id tag line
        if (line.begin() != line.end() && *(line.begin()) == '>' && (line.begin() + 1) != line.end())
        {

            // if we read primary assembly previously, write it to vectors
            if (contFlag)
            {

                seqFwd.shrink_to_fit();
                seqRev.shrink_to_fit();
                refSeqFwd.emplace_back(move(seqFwd));
                std::reverse(seqRev.begin(), seqFwd.end());
                refSeqRev.emplace_back(move(seqRev));
                // reset buffer
                seqFwd = std::string();
                seqFwd.reserve(chromLen);
                seqRev = std::string();
                seqRev.reserve(chromLen);
                lastC = false;

            }
            // check if primary assembly
            if (*(line.begin() + 1) == 'C')
            {

                ++chrIndex;
                contFlag = true;
                continue;

            // throw out unlocalized contigs
            } else {

                contFlag = false;
                continue;
            }

        }

        // check if we are in real primary chromosome assembly
        if (contFlag)
        {

            // parse the line
            for (char& c : line)
            {

                switch (c)
                {

                    case 'g':
                    case 'G':

                        if (lastC)
                        {
                            cpgMethRateFwd[(static_cast<uint64_t>(chrIndex - 1) << 32 | seqFwd.size())] = {0,0};
                            cpgMethRateRev[(static_cast<uint64_t>(chrIndex - 1) << 32 | seqFwd.size())] = {0,0};
                            lastC = false;
                        }
                        seqFwd.push_back('G');
                        seqRev.push_back('C');
                        break;

                    case 'c':
                    case 'C':

                        seqFwd.push_back('C');
                        seqRev.push_back('G');
                        lastC = true;
                        break;

                    case 'a':
                    case 'A':

                        seqFwd.push_back('A');
                        seqRev.push_back('T');
                        lastC = false;
                        break;

                    case 't':
                    case 'T':

                        seqFwd.push_back('T');
                        seqRev.push_back('A');
                        lastC = false;
                        break;

                    default:

                        seqFwd.push_back('N');
                        seqRev.push_back('N');
                        lastC = false;
                        break;


                }
            }

        } else {

            continue;
        }
    }

    // if we read primary assembly previously, write it to vectors
    if (contFlag)
    {

        seqFwd.shrink_to_fit();
        seqRev.shrink_to_fit();
        refSeqFwd.emplace_back(move(seqFwd));
        std::reverse(seqRev.begin(), seqFwd.end());
        refSeqRev.emplace_back(move(seqRev));
    }

    std::cout << "Done reading reference file\n\n";
}



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
