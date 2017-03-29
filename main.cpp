#include "RefReader_istr.h"
#include "RefGenome.h"

// --------------- MAIN -----------------
//
//
int main(int argc, char** argv) {


    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStartTab;
    std::vector<std::vector<char>> genSeq;

    readReference(argv[1], cpgTab, cpgStartTab, genSeq);
    RefGenome ref(cpgTab, cpgStartTab, genSeq);
}
//
//
//---------------------------------------
