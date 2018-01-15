
#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>  // reverse

#include "SynthDS.h"

SynthDS::SynthDS(const size_t refLen) :
        toIndex(0,3)
    ,   pairedOffDist(pairedMinDist, pairedMaxDist)
    ,   hasZeroErr(0.75)
    ,   hasOneErr(0.7)
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
    ,   pairedOffDist(pairedMinDist, pairedMaxDist)
    ,   hasZeroErr(0.75)
    ,   hasOneErr(0.7)
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
    ,   pairedOffDist(pairedMinDist, pairedMaxDist)
    ,   hasZeroErr(0.75)
    ,   hasOneErr(0.7)
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

// std::vector<std::string> SynthDS::genReadsFwdFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<size_t>& offsets)
// {
//
//     // will hold the generated reads
//     std::vector<std::string> readSet(readNum);
//     offsets.resize(readNum);
//
//     // range of indices allowed to be drawn for the reference
//     std::uniform_int_distribution<int> toOffset (0, refFwd.size() - readLen);
//
//     // range of indices allowed to be drawn for errors in read
//     std::uniform_int_distribution<int> toOffRead (0, readLen - 1);
//
//     // range of errors allowed to be drawn
//     std::uniform_int_distribution<int> toErr(0, maxErrNum);
//
//     // coin flip for conversion type
//     std::uniform_int_distribution<int> coin(0, 1);
//
// #pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
//     for (size_t i = 0; i < readNum; ++i)
//     {
//         std::mt19937& MT = randGen[omp_get_thread_num()];
//         // draw an offset where to draw the read
//         const size_t offset = toOffset(MT);
//         // generate sequence
//         std::string read(refFwd, offset, readLen);
//
//         // draw number of errors
//         const unsigned int err = toErr(MT);
//         // introduce errors at random positions
//         for (unsigned int e = 0; e < err; ++e)
//         {
//             read[toOffRead(MT)] = alphabet[toIndex(MT)];
//
//         }
//         // introduce C->T  OR  G->A conversions
//         if (coin(MT))
//         {
//
//             for (size_t j = 0; j < readLen; ++j)
//             {
//                 if (read[j] == 'C')
//                 {
//                     // try conversion
//                     if (coin(MT))
//                     {
//                         read[j] = 'T';
//                     }
//
//                 }
//             }
//
//         } else {
//
//             for (size_t j = 0; j < readLen; ++j)
//             {
//                 if (read[j] == 'G')
//                 {
//                     // try conversion
//                     if (coin(MT))
//                     {
//                         read[j] = 'A';
//                     }
//
//                 }
//             }
//         }
//
//         readSet[i] = std::move(read);
//         offsets[i] = offset;
//     }
//
//     std::cout << "Generated forward strand read set\n\n";
//     return readSet;
// }
//
// std::vector<std::string> SynthDS::genReadsRevFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum)
// {
//
//     // will hold the generated reads
//     std::vector<std::string> readSet(readNum);
//
//     // range of indices allowed to be drawn for the reference
//     std::uniform_int_distribution<int> toOffset (0, refRev.size() - readLen);
//
//     // range of indices allowed to be drawn for errors in read
//     std::uniform_int_distribution<int> toOffRead (0, readLen - 1);
//
//     // range of errors allowed to be drawn
//     std::uniform_int_distribution<int> toErr(0, maxErrNum);
//
//     // coin flip for conversion type
//     std::uniform_int_distribution<int> coin(0, 1);
//
// #pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
//     for (size_t i = 0; i < readNum; ++i)
//     {
//         std::mt19937& MT = randGen[omp_get_thread_num()];
//         // draw an offset where to draw the read
//         const size_t offset = toOffset(MT);
//         // generate sequence
//         std::string read(refRev, offset, readLen);
//
//         // draw number of errors
//         const unsigned int err = toErr(MT);
//         // introduce errors at random positions
//         for (unsigned int e = 0; e < err; ++e)
//         {
//             read[toOffRead(MT)] = alphabet[toIndex(MT)];
//
//         }
//         // introduce C->T  OR  G->A conversions
//         if (coin(MT))
//         {
//
//             for (size_t j = 0; j < readLen; ++j)
//             {
//                 if (read[j] == 'C')
//                 {
//                     // try conversion
//                     if (coin(MT))
//                     {
//                         read[j] = 'T';
//                     }
//
//                 }
//             }
//
//         } else {
//
//             for (size_t j = 0; j < readLen; ++j)
//             {
//                 if (read[j] == 'G')
//                 {
//                     // try conversion
//                     if (coin(MT))
//                     {
//                         read[j] = 'A';
//                     }
//
//                 }
//             }
//         }
//         readSet[i] = std::move(read);
//     }
//
//     std::cout << "Generated reverse complement strand read set\n\n";
//     return readSet;
// }

unsigned int SynthDS::getErrNum()
{
    std::mt19937& MT = randGen[omp_get_thread_num()];
    if (hasZeroErr(MT))
        return 0;
    else if (hasOneErr(MT))
        return 1;
    else
        return 2;
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

        // if (chr != 21)
        //     continue;
        // const size_t chrReadNum = readNum;
        // compute proportion of reads that should be generated by this chromosome
        const size_t chrReadNum = ((double)refSeqFwd[chr].size() / (double)genomeSize) * readNum;

        // range of indices allowed to be drawn for the reference
        std::uniform_int_distribution<int> toOffset (0, refSeqFwd[chr].size() - readLen);

        // range of indices allowed to be drawn for errors in read
        std::uniform_int_distribution<int> toOffRead (0, readLen - 1);


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
            bool hasN = false;
            // search for CpGs
            for (const char c : read)
            {
                if (c == 'C')
                {
                    prevC = true;
                } else if (c == 'G') {

                    if (prevC)
                    {
                        hasCpG = true;
                        prevC = false;
                    }

                } else if (c == 'N') {

                    hasN = true;
                    break;
                }
            }

            // produce only reads with CpG
            if (!hasCpG || hasN)
            {
                --i;
                continue;
            }

            // introduce C->T  OR  G->A conversions
            if (coin(MT))
            {

                for (size_t j = 0; j < readLen; ++j)
                {
                    if (read[j] == 'C')
                    {
                        if (j < readLen - 1)
                        {
                            // if G follows, update CpG structure and methylate according to sampling rate
                            if (read[j+1] == 'G')
                            {
                                struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ];
                                std::bernoulli_distribution methIt(met.sampleRate);
                                if (!methIt(MT))
                                {
                                    read[j] = 'T';
#pragma omp atomic
                                    ++met.unmethCount;
                                } else {
#pragma omp atomic
                                    ++met.methCount;
                                }


                            } else {
                                // see if methylated
                                if (!methToss(MT))
                                {
                                    // if not, try bisulfite conversion
                                    if (convToss(MT))
                                        read[j] = 'T';
                                }
                            }
                        } else {

                            // see if methylated
                            if (!methToss(MT))
                            {
                                // if not, try bisulfite conversion
                                if (convToss(MT))
                                    read[j] = 'T';
                            }
                        }
                    }
                }

            } else {

                std::reverse(read.begin(), read.end());
                for (size_t j = 0; j < readLen; ++j)
                {
                    switch (read[j])
                    {
                        // C will generate G on reverse strand -> test for G to A conversion
                        case ('C') :

                            if (j > 0)
                            {
                                // letter j-1 was already processed and hence converted from G to C
                                if (read[j-1] == 'C')
                                {
                                    struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + readLen - j - 1) ];
                                    std::bernoulli_distribution methIt(met.sampleRate);
                                    if (!methIt(MT))
                                    {
                                        read[j] = 'A';
#pragma omp atomic
                                        ++met.unmethCount;
                                    } else {
                                        read[j] = 'G';
#pragma omp atomic
                                        ++met.methCount;
                                    }


                                } else {
                                    // see if methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read[j] = 'A';
                                        else
                                            read[j] = 'G';
                                    } else {

                                        read[j] = 'G';
                                    }
                                }
                            } else {
                                // see if methylated
                                if (!methToss(MT))
                                {
                                    // if not, try bisulfite conversion
                                    if (convToss(MT))
                                        read[j] = 'A';
                                    else
                                        read[j] = 'G';
                                } else {

                                    read[j] = 'G';
                                }
                            }
                            break;

                        case ('G') :

                            read[j] = 'C';
                            break;

                        case ('A') :

                            read[j] = 'T';
                            break;

                        case ('T') :

                            read[j] = 'A';
                            break;
                    }
                }
            }
            // draw number of errors
            const unsigned int err = getErrNum();
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

            readSet[i] = std::move(read);
            offsets[i] = std::pair<size_t,size_t>(offset, chr);
        }
        processedReads += chrReadNum;

    }
    std::cout << "Generated forward strand read set\n\n";
    return readSet;
}

// std::vector<std::string> SynthDS::genReadsRevRef(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<std::pair<size_t, size_t> >& offsets, std::vector<std::array<int, errNum> >& errOffs)
// {
//
//     std::cout << "Start generating reverse strand read set\n\n";
//     // will hold the generated reads
//     std::vector<std::string> readSet(readNum);
//     offsets.resize(readNum);
//     errOffs.resize(readNum);
//
//     size_t processedReads = 0;
//
//     // size of the genome
//     size_t genomeSize = 0;
//     for (size_t chr = 0; chr < refSeqFwd.size(); ++chr)
//         genomeSize += refSeqFwd[chr].size();
//
//     // go over all chromosomes
//     for (size_t chr = 0; chr < refSeqFwd.size(); ++chr)
//     {
//
//         if (chr != 21)
//             continue;
//         // compute proportion of reads that should be generated by this chromosome
//         // const size_t chrReadNum = ((double)refSeqFwd[chr].size() / (double)genomeSize) * readNum;
//         const size_t chrReadNum = readNum;
//
//         // range of indices allowed to be drawn for the reference
//         std::uniform_int_distribution<int> toOffset (0, refSeqFwd[chr].size() - readLen);
//
//         // range of indices allowed to be drawn for errors in read
//         std::uniform_int_distribution<int> toOffRead (0, readLen - 1);
//
//         // range of errors allowed to be drawn
//         std::uniform_int_distribution<int> toErr(0, maxErrNum);
//
//         // coin flip for conversion type (G to A or C to T)
//         // i.e. forward or reverse of PCR of genomic fragment treated with bisulfite
//         std::uniform_int_distribution<int> coin(0, 1);
//
// #pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
//         for (size_t i = processedReads; i < ((chr == refSeqFwd.size() - 1) ? readNum : processedReads + chrReadNum); ++i)
//         {
//
//             std::mt19937& MT = randGen[omp_get_thread_num()];
//             // draw an offset where to draw the read
//             const size_t offset = toOffset(MT);
//             // generate sequence
//             std::string read(refSeqFwd[chr], offset, readLen);
//             // test if CpG is contained
//             bool prevC = false;
//             bool hasCpG = false;
//             // search for CpGs
//             for (const char c : read)
//             {
//                 if (c == 'C')
//                 {
//                     prevC = true;
//                     continue;
//                 }
//                 if (c == 'G')
//                 {
//                     if (prevC)
//                         hasCpG = true;
//                 }
//                 prevC = false;
//             }
//
//             // produce only reads with CpG
//             if (!hasCpG)
//             {
//                 --i;
//                 continue;
//             }
//
//             // introduce C->T  OR  G->A conversions
//             if (coin(MT))
//             {
//
//                 for (size_t j = 0; j < readLen; ++j)
//                 {
//                     if (read[j] == 'G')
//                     {
//                         // try conversion
//                         if (!methToss(MT))
//                         {
//                             if (convToss(MT))
//                                 read[j] = 'A';
//                             // if necessary, count up methylation counters
//                             if (j > 0)
//                             {
//                                 if (read[j-1] == 'C')
//                                 {
// #pragma omp atomic
//                                     cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + j - 1) ].first += 1;
//                                 }
//                             }
//                         } else {
//                             // if necessary, count up methylation counters
//                             if (j > 0)
//                             {
//                                 if (read[j-1] == 'C')
//                                 {
// #pragma omp atomic
//                                     cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + j - 1) ].second += 1;
//                                 }
//                             }
//                         }
//                     }
//                 }
//
//             } else {
//
//                 std::reverse(read.begin(), read.end());
//                 for (size_t j = 0; j < readLen; ++j)
//                 {
//                     switch (read[j])
//                     {
//                         // G will generate C on reverse strand -> test for C to T conversion
//                         case ('G') :
//
//                             // try conversion
//                             if (!methToss(MT))
//                             {
//                                 if (convToss(MT))
//                                     read[j] = 'T';
//                                 else
//                                     read[j] = 'C';
//
//                                 // if necessary, count up methylation counters
//                                 if (j > 0)
//                                 {
//                                     if (read[j-1] == 'C')
//                                     {
// #pragma omp atomic
//                                         cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + readLen - j - 1) ].first += 1;
//                                     }
//                                 }
//
//                             } else {
//
//                                 // if necessary, count up methylation counters
//                                 read[j] = 'C';
//                                 if (j > 0)
//                                 {
//                                     if (read[j-1] == 'C')
//                                     {
// #pragma omp atomic
//                                         cpgMethRateRev[ (static_cast<uint64_t>(chr) << 32) | (offset + readLen - j - 1) ].second += 1;
//                                     }
//                                 }
//
//                             }
//                             break;
//
//                         case ('C') :
//
//                             read[j] = 'G';
//                             break;
//
//                         case ('A') :
//
//                             read[j] = 'T';
//                             break;
//
//                         case ('T') :
//
//                             read[j] = 'A';
//                             break;
//                     }
//                 }
//             }
//
//             // draw number of errors
//             const unsigned int err = toErr(MT);
//             // introduce errors at random positions
//             size_t errID = 0;
//             for (unsigned int e = 0; e < err; ++e)
//             {
//                 int eOff = toOffRead(MT);
//                 read[eOff] = alphabet[toIndex(MT)];
//                 errOffs[i][errID] = eOff;
//                 ++errID;
//             }
//             for (; errID < errNum; ++errID)
//             {
//                 errOffs[i][errID] = -1;
//             }
//
//             readSet[i] = std::move(read);
//             offsets[i] = std::move(std::pair<size_t,size_t>(offset, chr));
//         }
//         processedReads += chrReadNum;
//
//     }
//     std::cout << "Generated reverse strand read set\n\n";
//     return readSet;
// }
//

std::pair<std::vector<std::string> > SynthDS::genReadsPairedRef(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::pair<std::vector<std::pair<size_t, size_t> > >& offsets, std::pair<std::vector<std::array<int, errNum> > >& errOffs)
{

    std::cout << "Start generating paired read set\n\n";
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

        // for single chromosome data set
        // if (chr != 21)
        //     continue;
        // const size_t chrReadNum = readNum;
        //
        // compute proportion of reads that should be generated by this chromosome
        const size_t chrReadNum = ((double)refSeqFwd[chr].size() / (double)genomeSize) * readNum;

        // range of indices allowed to be drawn for the reference
        std::uniform_int_distribution<int> toOffset (0, refSeqFwd[chr].size() - readLen - pairedMaxDist);

        // range of indices allowed to be drawn for errors in read
        std::uniform_int_distribution<int> toOffRead (0, readLen - 1);


        // coin flip for conversion type (G to A or C to T)
        // i.e. forward or reverse of PCR of genomic fragment treated with bisulfite
        std::uniform_int_distribution<int> coin(0, 1);

#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
        for (size_t i = processedReads; i < ((chr == refSeqFwd.size() - 1) ? readNum : processedReads + chrReadNum); ++i)
        {

            std::mt19937& MT = randGen[omp_get_thread_num()];
            // draw an offset where to draw the first read
            const size_t offset = toOffset(MT);
            // draw a distance between the two reads
            const size_t pDist = pairedOffDist(MT);
            std::string read1(refSeqFwd[chr], offset, readLen);
            std::string read2(refSeqFwd[chr], offset + pDist, readLen);
            // make read2 reverse complementary
            size_t readPos = 0;
            bool read1isFwd;
            // flip a coin if read 1 is from main strand
            if (coin(MT))
            {

                read1isFwd = true;

            // read 2 is on main strand
            } else {

                read1isFwd = false;
            }


            // test if CpG is contained
            bool prevC = false;
            bool hasCpG = false;
            // test if N is contained
            bool hasN = false;
            // search for CpGs in both reads
            for (const char c : read1)
            {
                if (c == 'C')
                {
                    prevC = true;
                    continue;
                } else if (c == 'N') {

                    hasN = true;
                    break;
                } else if (c == 'G') {

                    if (prevC)
                        hasCpG = true;
                }
                prevC = false;
            }
            prevC = false;
            for (const char c : read2)
            {
                if (c == 'C')
                {
                    prevC = true;
                    continue;
                } else if (c == 'N') {

                    hasN = true;
                    break;
                } else if (c == 'G') {

                    if (prevC)
                        hasCpG = true;
                }
                prevC = false;
            }

            // produce only reads where at least one of the two has a CpG
            if (!hasCpG || hasN)
            {
                --i;
                continue;
            }

            // draw if main strand or second strand is originally methylated
            // i.e. if read on main strand is C->T or G->A
            if (coin(MT))
            {

                if (read1isFwd)
                {
                    // go through read 1 and introduce methylation
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        if (read1[j] == 'C')
                        {
                            if (j < readLen - 1)
                            {
                                // if G follows, update CpG structure and methylate according to sampling rate
                                if (read1[j+1] == 'G')
                                {
                                    struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + j) ];
                                    std::bernoulli_distribution methIt(met.sampleRate);
                                    if (!methIt(MT))
                                    {
                                        read1[j] = 'T';
#pragma omp atomic
                                        ++met.unmethCount;
                                    } else {
#pragma omp atomic
                                        ++met.methCount;
                                    }


                                } else {

                                    // see if methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read1[j] = 'T';
                                    }
                                }
                            } else {

                                // see if methylated
                                if (!methToss(MT))
                                {
                                    // if not, try bisulfite conversion
                                    if (convToss(MT))
                                        read1[j] = 'T';
                                }
                            }
                        }
                    }
                    // go through read 2 and introduce methylation
                    std::reverse(read2.begin(), read2.end());
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        switch (read2[j])
                        {
                            // C will generate G on reverse strand -> test for G to A conversion
                            case ('C') :

                                if (j > 0)
                                {
                                    // letter j-1 was already processed and hence converted from G to C
                                    if (read2[j-1] == 'C')
                                    {
                                        struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + readLen + pDist - j - 1) ];
                                        std::bernoulli_distribution methIt(met.sampleRate);
                                        if (!methIt(MT))
                                        {
                                            read2[j] = 'A';
#pragma omp atomic
                                            ++met.unmethCount;
                                        } else {
                                            read2[j] = 'G';
#pragma omp atomic
                                            ++met.methCount;
                                        }


                                    } else {
                                        // see if methylated
                                        if (!methToss(MT))
                                        {
                                            // if not, try bisulfite conversion
                                            if (convToss(MT))
                                                read2[j] = 'A';
                                            else
                                                read2[j] = 'G';
                                        } else {

                                            read2[j] = 'G';
                                        }
                                    }
                                } else {
                                    // see if methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read2[j] = 'A';
                                        else
                                            read2[j] = 'G';
                                    } else {

                                        read2[j] = 'G';
                                    }
                                }
                                break;

                            case ('G') :

                                read2[j] = 'C';
                                break;

                            case ('A') :

                                read2[j] = 'T';
                                break;

                            case ('T') :

                                read2[j] = 'A';
                                break;
                        }
                    }


                // read 2 is on main strand
                } else {

                    // go through read 2 and introduce methylation
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        if (read2[j] == 'C')
                        {
                            if (j < readLen - 1)
                            {
                                // if G follows, update CpG structure and methylate according to sampling rate
                                if (read2[j+1] == 'G')
                                {
                                    struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + pDist + j) ];
                                    std::bernoulli_distribution methIt(met.sampleRate);
                                    if (!methIt(MT))
                                    {
                                        read2[j] = 'T';
#pragma omp atomic
                                        ++met.unmethCount;
                                    } else {
#pragma omp atomic
                                        ++met.methCount;
                                    }


                                } else {

                                    // see if methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read2[j] = 'T';
                                    }
                                }
                            } else {

                                // see if methylated
                                if (!methToss(MT))
                                {
                                    // if not, try bisulfite conversion
                                    if (convToss(MT))
                                        read2[j] = 'T';
                                }
                            }
                        }
                    }
                    // go through read 2 and introduce methylation
                    std::reverse(read1.begin(), read1.end());
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        switch (read1[j])
                        {
                            // C will generate G on reverse strand -> test for G to A conversion
                            case ('C') :

                                if (j > 0)
                                {
                                    // letter j-1 was already processed and hence converted from G to C
                                    if (read1[j-1] == 'C')
                                    {
                                        struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + readLen - j - 1) ];
                                        std::bernoulli_distribution methIt(met.sampleRate);
                                        if (!methIt(MT))
                                        {
                                            read1[j] = 'A';
#pragma omp atomic
                                            ++met.unmethCount;
                                        } else {
                                            read1[j] = 'G';
#pragma omp atomic
                                            ++met.methCount;
                                        }


                                    } else {
                                        // see if methylated
                                        if (!methToss(MT))
                                        {
                                            // if not, try bisulfite conversion
                                            if (convToss(MT))
                                                read1[j] = 'A';
                                            else
                                                read1[j] = 'G';
                                        } else {

                                            read1[j] = 'G';
                                        }
                                    }
                                } else {
                                    // see if methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read1[j] = 'A';
                                        else
                                            read1[j] = 'G';
                                    } else {

                                        read1[j] = 'G';
                                    }
                                }
                                break;

                            case ('G') :

                                read1[j] = 'C';
                                break;

                            case ('A') :

                                read1[j] = 'T';
                                break;

                            case ('T') :

                                read1[j] = 'A';
                                break;
                        }
                    }
                }

            // (!coin(MT))
            } else {

                if (read1isFwd)
                {
                    // go through read 1 and introduce methylation
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        if (read1[j] == 'G')
                        {
                            if (j > 0)
                            {
                                // if G follows, update CpG structure and methylate according to sampling rate
                                if (read1[j - 1] == 'C')
                                {
                                    struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + j - 1) ];
                                    std::bernoulli_distribution methIt(met.sampleRate);
                                    if (!methIt(MT))
                                    {
                                        read1[j] = 'A';
#pragma omp atomic
                                        ++met.unmethCount;
                                    } else {
#pragma omp atomic
                                        ++met.methCount;
                                    }

                                } else {

                                    // see if not methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read1[j] = 'A';
                                    }
                                }
                            } else {

                                // see if not methylated
                                if (!methToss(MT))
                                {
                                    // if not, try bisulfite conversion
                                    if (convToss(MT))
                                        read1[j] = 'A';
                                }
                            }
                        }
                    }
                    // go through read 2 and introduce methylation
                    std::reverse(read2.begin(), read2.end());
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        switch (read2[j])
                        {
                            case ('G') :

                                if (j < readLen - 1)
                                {
                                    if (read2[j + 1] == 'C')
                                    {
                                        struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + readLen + pDist - j) ];
                                        std::bernoulli_distribution methIt(met.sampleRate);
                                        if (!methIt(MT))
                                        {
                                            read2[j] = 'T';
#pragma omp atomic
                                            ++met.unmethCount;
                                        } else {
                                            read2[j] = 'C';
#pragma omp atomic
                                            ++met.methCount;
                                        }


                                    } else {
                                        // see if methylated
                                        if (!methToss(MT))
                                        {
                                            // if not, try bisulfite conversion
                                            if (convToss(MT))
                                                read2[j] = 'T';
                                            else
                                                read2[j] = 'C';
                                        } else {

                                            read2[j] = 'C';
                                        }
                                    }
                                } else {
                                    // see if methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read2[j] = 'T';
                                        else
                                            read2[j] = 'C';
                                    } else {

                                        read2[j] = 'C';
                                    }
                                }
                                break;

                            case ('C') :

                                read2[j] = 'G';
                                break;

                            case ('A') :

                                read2[j] = 'T';
                                break;

                            case ('T') :

                                read2[j] = 'A';
                                break;
                        }
                    }


                // read 2 is on main strand
                } else {

                    // go through read 2 and introduce methylation
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        if (read2[j] == 'G')
                        {
                            if (j > 0)
                            {
                                // if G follows, update CpG structure and methylate according to sampling rate
                                if (read2[j - 1] == 'C')
                                {
                                    struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + pDist + j - 1) ];
                                    std::bernoulli_distribution methIt(met.sampleRate);
                                    if (!methIt(MT))
                                    {
                                        read2[j] = 'A';
#pragma omp atomic
                                        ++met.unmethCount;
                                    } else {
#pragma omp atomic
                                        ++met.methCount;
                                    }

                                } else {

                                    // see if not methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read2[j] = 'A';
                                    }
                                }
                            } else {

                                // see if not methylated
                                if (!methToss(MT))
                                {
                                    // if not, try bisulfite conversion
                                    if (convToss(MT))
                                        read2[j] = 'A';
                                }
                            }
                        }
                    }
                    // go through read 2 and introduce methylation
                    std::reverse(read1.begin(), read1.end());
                    for (size_t j = 0; j < readLen; ++j)
                    {
                        switch (read1[j])
                        {
                            case ('G') :

                                if (j < readLen - 1)
                                {
                                    if (read1[j + 1] == 'C')
                                    {
                                        struct MethInfo& met = cpgMethRateFwd[ (static_cast<uint64_t>(chr) << 32) | (offset + readLen - j) ];
                                        std::bernoulli_distribution methIt(met.sampleRate);
                                        if (!methIt(MT))
                                        {
                                            read1[j] = 'T';
#pragma omp atomic
                                            ++met.unmethCount;
                                        } else {
                                            read1[j] = 'C';
#pragma omp atomic
                                            ++met.methCount;
                                        }


                                    } else {
                                        // see if methylated
                                        if (!methToss(MT))
                                        {
                                            // if not, try bisulfite conversion
                                            if (convToss(MT))
                                                read1[j] = 'T';
                                            else
                                                read1[j] = 'C';
                                        } else {

                                            read1[j] = 'C';
                                        }
                                    }
                                } else {
                                    // see if methylated
                                    if (!methToss(MT))
                                    {
                                        // if not, try bisulfite conversion
                                        if (convToss(MT))
                                            read1[j] = 'T';
                                        else
                                            read1[j] = 'C';
                                    } else {

                                        read1[j] = 'C';
                                    }
                                }
                                break;

                            case ('C') :

                                read1[j] = 'G';
                                break;

                            case ('A') :

                                read1[j] = 'T';
                                break;

                            case ('T') :

                                read1[j] = 'A';
                                break;
                        }
                    }
                }
            }

            // draw number of errors
            const unsigned int err = getErrNum();
            // introduce errors at random positions
            size_t errID = 0;
            // for read 1
            for (unsigned int e = 0; e < err; ++e)
            {
                int eOff = toOffRead(MT);
                read[eOff] = alphabet[toIndex(MT)];
                errOffs.first[i][errID] = eOff;
                ++errID;
            }
            for (; errID < errNum; ++errID)
            {
                errOffs.first[i][errID] = -1;
            }

            err = getErrNum();
            errID = 0;
            for (unsigned int e = 0; e < err; ++e)
            {
                int eOff = toOffRead(MT);
                read[eOff] = alphabet[toIndex(MT)];
                errOffs.second[i][errID] = eOff;
                ++errID;
            }
            for (; errID < errNum; ++errID)
            {
                errOffs.second[i][errID] = -1;
            }

            // shuffle pairs
            if (coin(MT))
            {
                readSet.first[i] = std::move(read1);
                offsets.first[i] = std::pair<size_t,size_t>(offset, chr);

                readSet.second[i] = std::move(read2);
                offsets.second[i] = std::pair<size_t,size_t>(offset + pDist, chr);

            } else {

                readSet.second[i] = std::move(read1);
                offsets.second[i] = std::pair<size_t,size_t>(offset, chr);

                readSet.first[i] = std::move(read2);
                offsets.first[i] = std::pair<size_t,size_t>(offset + pDist, chr);
            }
        }
        processedReads += chrReadNum;

    }
    std::cout << "Generated forward strand read set\n\n";
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

    // fair coin flip
    std::uniform_int_distribution<int> coin(0, 1);
    std::normal_distribution<double> unmethDist(0.2,0.08);
    std::normal_distribution<double> methDist(0.8,0.08);

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
                            // sample the methylation rates from bimodal distribution
                            double methRate = coin(randGen[0]) ? methDist(randGen[0]) : unmethDist(randGen[0]);
                            if (methRate < 0)
                                methRate = 0;
                            else if (methRate > 1)
                                methRate = 1;
                            cpgMethRateFwd[(static_cast<uint64_t>(chrIndex - 1) << 32 | (seqFwd.size() - 1))] = {0,0,methRate};
                            cpgMethRateRev[(static_cast<uint64_t>(chrIndex - 1) << 32 | (seqFwd.size() - 1))] = {0,0,methRate};
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
