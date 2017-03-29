#include "RefReader_istr.h"
#include "RefGenome.h"

// --------------- MAIN -----------------
//
//
int main(int argc, char** argv) {


    std::vector<struct CpG> cpgTab;
    std::vector<std::vector<char>> genSeq;

    readReference(argv[1], cpgTab, genSeq);
    RefGenome ref(cpgTab, genSeq);
}
//
//
//---------------------------------------
