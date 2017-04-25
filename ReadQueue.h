#ifndef READQUEUE_H
#define READQUEUE_H

#include <string>
#include <fstream>
#include <unordered_map>

#include "CONST.h"
#include "RefGenome.h"
#include "Read.h"

class ReadQueue
{

    public:

        ReadQueue(std::string& filePath, RefGenome& ref);


        // Parses a chunk of the ifstream file
        // Reads MyConst::CHUNKSIZE many reads and saves them
        // returns true if neither read error nor EOF occured, false otherwise
        bool parseChunk();

        // Match all reads in readBuffer to reference genome
        // First retrieve seeds using getSeeds(...)
        // Filter seeds using filterHeuSeeds(...) according to simple heuristic
        // Extend remaining seeds with BitMatch(...)
        // Filter result with filterHeuRes(...) according to cumulative sum heuristic
        bool matchReads();


    private:

        // filters seeds according to simple counting criteria
        // #kmers of one metaCpG should be > READLEN - (KMERLEN * MISCOUNT)
        inline bool filterSeeds(std::vector<std::vector<KMER::kmer> >& seeds, unsigned int readSize)
        {
            // containing counts (second template param) on how often a specific metaCpG (index is first template param aka key)
            // is found across all kmers of this read
            std::unordered_map<uint32_t, unsigned int> counts;

            // count occurences of meta CpGs
            for (unsigned int i = 0; i < seeds.size(); ++i)
            {

                for (unsigned int j = 0; j < seeds[0].size(); ++j)
                {

                    auto insert = counts.emplace(KMER::getMetaCpG(seeds[i][j]), 1);

                    // if insertion fails because metaCpG is already inserted, count one up
                    if (!insert.second)
                    {
                        ++((*(insert.first)).second);
                    }
                }
            }

            // More than cutoff many kmers are required per metaCpG
            const unsigned int cutoff = readSize - (MyConst::KMERLEN * MyConst::MISCOUNT);

            bool nonEmpty = false;
            // throw out rare metaCpGs
            for (unsigned int i = 0; i < seeds.size(); ++i)
            {

                std::vector<KMER::kmer> filteredSeeds;
                filteredSeeds.reserve(seeds[i].size());
                for (unsigned int j = 0; j < seeds[0].size(); ++j)
                {

                    if (counts[KMER::getMetaCpG(seeds[i][j])] > cutoff)
                    {

                        filteredSeeds.push_back(seeds[i][j]);

                    }
                }
                filteredSeeds.shrink_to_fit();
                nonEmpty = nonEmpty || !filteredSeeds.empty();
                seeds[i] = std::move(filteredSeeds);
            }

            return nonEmpty;
        }


        // input stream of file given as path to Ctor
        std::ifstream file;

        // representation of the reference genome
        RefGenome& ref;

        // buffer holding MyConst::CHUNKSIZE many reads
        std::vector<Read> readBuffer;

};

#endif /* READQUEUE_H */
