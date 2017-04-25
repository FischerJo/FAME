#ifndef READQUEUE_H
#define READQUEUE_H

#include <string>
#include <fstream>

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

        // input stream of file given as path to Ctor
        std::ifstream file;

        // representation of the reference genome
        RefGenome& ref;

        // buffer holding MyConst::CHUNKSIZE many reads
        std::vector<Read> readBuffer;

};

#endif /* READQUEUE_H */
