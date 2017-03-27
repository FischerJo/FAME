#include "RefReader.h"
#include "RefGenome.h"

// --------------- MAIN -----------------
//
//
int main(int argc, char** argv) {


    std::vector<struct CpG> cpgTab;
    std::vector<const char*> genSeq;
    std::vector<std::size_t> genSeq;

    readReference(argv[1], cpgTab, genSeq, genSeqLen);
    RefGenome ref(cpgTab, genSeq, genSeqLen);
}
//
//
//---------------------------------------
