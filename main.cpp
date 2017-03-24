#include "RefReader.h"
#include "RefGenome.h"

// --------------- MAIN -----------------
//
//
int main(int argc, char** argv) {


    std::vector<struct CpG> cpgTab;
    std::vector<std::string> genSeq;

    readReference(argv[1], cpgTab, genSeq);
    RefGenome ref(cpgTab, genSeq);
}
//
//
//---------------------------------------
