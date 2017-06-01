#include "RefReader_istr.h"
#include "RefGenome.h"
#include "ReadQueue.h"

// --------------- MAIN -----------------
//
//
int main(int argc, char** argv) {


    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char>> genSeq;

    readReference(argv[1], cpgTab, cpgStartTab, genSeq);
    RefGenome ref(std::move(cpgTab), std::move(cpgStartTab), genSeq);
    ReadQueue rQue(argv[2], ref);

    unsigned int readCounter = 0;
    unsigned int i = 0;

    while(rQue.parseChunk(readCounter))
    {
        rQue.matchReads(readCounter);
        std::cout << "Processed " << MyConst::CHUNKSIZE * (++i) << " many reads\n";
    }
    // match remaining reads
    rQue.matchReads(readCounter);

    return 0;
}
//
//
//---------------------------------------
