#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ReadQueue.h"

ReadQueue::ReadQueue(std::string& filePath, RefGenome& reference) :
        file(filePath)
    ,   ref(reference)
    ,   readBuffer()
{
    // reserve space for the buffer object according to user parameter
    readBuffer.reserve(MyConst::CHUNKSIZE);
}

bool ReadQueue::parseChunk()
{

    std::string id;

    // counter on how many reads have been read so far
    unsigned int readCounter = 0;

    // read first line of read (aka @'SEQID')
    while (std::getline(file, id))
    {

        // read the next line (aka raw sequence)
        std::string seq;
        std::getline(file, seq);
        // construct read and push it to buffer
        readBuffer.emplace_back(seq, id);
        // read the rest of read (aka +'SEQID' and quality score sequence)
        std::getline(file,id);
        std::getline(file,seq);

        ++readCounter;

        // if buffer is full, return
        if (readCounter >= MyConst::CHUNKSIZE)
        {
            return true;

        }
    }

    return false;
}

bool ReadQueue::matchReads()
{

#ifdef _OPENMP
#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,10)
#endif
    for (unsigned int i = 0; i < MyConst::CHUNKSIZE; ++i)
    {

        Read& r = readBuffer[i];

        unsigned int readSize = r.seq.size();

        // RETRIEVE SEEDS
        //
        std::vector<char> redSeq (readSize);
        std::vector<char> redRevSeq (readSize);

        // construct reduced alphabet sequence for forward and reverse strand
        for (unsigned int pos = 0; pos < readSize; ++pos)
        {
            // get correct offset for reverse strand (strand orientation must be correct)
            unsigned int revPos = readSize - 1 - pos;

            switch (r.seq[pos])
            {
                case 'A':

                    redSeq[pos] = 'A';
                    redRevSeq[revPos] = 'T';
                    break;

                case 'C':

                    redSeq[pos] = 'T';
                    redRevSeq[revPos] = 'G';
                    break;

                case 'G':

                    redSeq[pos] = 'G';
                    redRevSeq[revPos] = 'T';
                    break;

                case 'T':

                    redSeq[pos] = 'T';
                    redRevSeq[revPos] = 'A';
                    break;

                default:

                    std::cerr << "Unknown character '" << r.seq[pos] << "' in read with sequence id " << r.id << std::endl;
            }
        }

        std::vector<std::vector<KMER::kmer> > fwdSeedsK;
        std::vector<std::vector<bool> > fwdSeedsS;
        ref.getSeeds(redSeq, fwdSeedsK, fwdSeedsS);
        std::vector<std::vector<KMER::kmer> > revSeedsK;
        std::vector<std::vector<bool> > revSeedsS;
        ref.getSeeds(redRevSeq, revSeedsK, revSeedsS);


        // FILTER SEEDS
        bool isFwd = filterSeeds(fwdSeedsK, fwdSeedsS, readSize);
        bool isRev = filterSeeds(revSeedsK, revSeedsS, readSize);

    }
    return true;
}
