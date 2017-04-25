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
#pragma omp parallel for num_threads(CORENUM)
#endif
    for (unsigned int i = 0; i < MyConst::CHUNKSIZE; ++i)
    {

        Read& r = readBuffer[i];

        // RETRIEVE SEEDS
        //
        std::vector<char> redSeq (r.seq.size());
        std::vector<char> redRevSeq (r.seq.size());

        ref.getSeeds(redSeq);
        ref.getSeeds(redRevSeq);

    }
    return true;
}
