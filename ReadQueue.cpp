#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ReadQueue.h"

ReadQueue::ReadQueue(const char* filePath, RefGenome& reference) :
        statFile("seedStats.tsv")
    ,   countFile("seedCounts.tsv")
    ,   file(filePath)
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

        // if buffer is read completely, return
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

        const unsigned int readSize = r.seq.size();

        if (readSize < MyConst::KMERLEN)
        {

            r.isInvalid = true;
            continue;
        }
        // RETRIEVE SEEDS
        //
        std::vector<char> redSeq (readSize);
        std::vector<char> redRevSeq (readSize);

        // flag stating if read contains N
        // reads with N are ignored
        bool nflag = false;

        // get correct offset for reverse strand (strand orientation must be correct)
        unsigned int revPos = readSize - 1;

        // construct reduced alphabet sequence for forward and reverse strand
        for (unsigned int pos = 0; pos < readSize; ++pos, --revPos)
        {

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

                case 'N':

                    nflag = true;
                    break;

                default:

                    std::cerr << "Unknown character '" << r.seq[pos] << "' in read with sequence id " << r.id << std::endl;
            }
        }

        // if N is present in read, skip
        if (nflag)
        {
            r.isInvalid = true;
            continue;
        }


        std::vector<std::vector<KMER::kmer> > fwdSeedsK;
        std::vector<std::vector<bool> > fwdSeedsS;
        ref.getSeeds(redSeq, fwdSeedsK, fwdSeedsS);
        std::vector<std::vector<KMER::kmer> > revSeedsK;
        std::vector<std::vector<bool> > revSeedsS;
        ref.getSeeds(redRevSeq, revSeedsK, revSeedsS);

        // printStatistics(fwdSeedsK);

        // FILTER SEEDS BY COUNTING LEMMA
        filterHeuSeeds(fwdSeedsK, fwdSeedsS, readSize);
        filterHeuSeeds(revSeedsK, revSeedsS, readSize);

        // printStatistics(fwdSeedsK);

        // BITMASK COMPARISON ON FULL ALPHABET FOR REMAINING SEEDS
        bitMatching(r, fwdSeedsK, fwdSeedsS);
        bitMatchingRev(r, revSeedsK, revSeedsS);

        // printStatistics(fwdSeedsK);
        // FILTER MATCHING KMERS BY COUNTING LEMMA
        filterHeuSeeds(fwdSeedsK, fwdSeedsS, readSize);
        filterHeuSeeds(revSeedsK, revSeedsS, readSize);

        // printStatistics(fwdSeedsK);

        // finish line for this read in counter file
        countFile << "\n";
    }
    return true;
}


void ReadQueue::printStatistics(std::vector<std::vector<KMER::kmer> > seedsK)
{

    unsigned int count = 0;

    // statistics on metaCpG repetitions
    std::vector<unsigned int> metaRep (n, 0);

    // count occurences of meta CpGs
    for (unsigned int i = 0; i < seedsK.size(); ++i)
    {

        std::unordered_map<uint32_t, unsigned int> metaCounts;

        // count up general seed counter
        count += seedsK[i].size();
        for (unsigned int j = 0; j < seedsK[i].size(); ++j)
        {

            uint64_t metaId = KMER::getMetaCpG(seedsK[i][j]);

            auto insert = metaCounts.emplace(metaId, 1);

            // if insertion fails because metaCpG is already inserted, count up everything
            if (!insert.second)
            {
                ++((insert.first)->second);
            }
        }
        for (auto& c : metaCounts)
        {
            if (c.second <= n)
            {

                ++metaRep[c.second - 1];

            }
        }

    }

    // print everything
    for (unsigned int i = 0; i <= metaRep.size(); ++i)
    {
        statFile << metaRep[i] << "\t";
    }
    statFile << "\n";

    countFile << count << "\t";

}
