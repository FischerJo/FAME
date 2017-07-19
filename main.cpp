
#include <chrono>


#include "RefReader_istr.h"
#include "RefGenome.h"
#include "ReadQueue.h"

// --------------- MAIN -----------------
//
//
int main(int argc, char** argv)
{


    std::string binFile (argv[3]);
    // std::vector<struct CpG> cpgTab;
    // std::vector<struct CpG> cpgStartTab;
    // std::vector<std::vector<char>> genSeq;
    //
    // readReference(argv[1], cpgTab, cpgStartTab, genSeq);
    // RefGenome ref(std::move(cpgTab), std::move(cpgStartTab), genSeq);
    // ref.save(binFile);
    RefGenome refCopy(binFile);
    ReadQueue rQue(argv[2], refCopy, false);
    // ReadQueue rQue(argv[2], refCopy, true);
    // ReadQueue rQue(argv[2], ref, false);

    unsigned int readCounter = 0;
    unsigned int i = 0;

    rQue.parseChunk(readCounter);
    rQue.matchReads(readCounter);
    // rQue.parseChunkGZ(readCounter);
    // rQue.matchReads(readCounter);
    // std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    // while(rQue.parseChunk(readCounter))
    // {
    //     rQue.matchReads(readCounter);
    //     std::cout << "Processed " << MyConst::CHUNKSIZE * (++i) << " many reads\n";
    // }
    // // match remaining reads
    // rQue.matchReads(readCounter);
    //
    // std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    // auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    //
    // std::cout << "Done processing in " << runtime << "s\n";

    return 0;
}
//
//
//---------------------------------------
