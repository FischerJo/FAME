

#include "RefGenome.h"


using namespace std;

RefGenome::RefGenome(unsigned int cpgs, vector<struct CpG>& cpgTab, vector<string>& genSeq) :
        cpgNum(cpgs)
    ,   cpgTable(cpgTab)
    ,   genomeSeq(genSeq)
{
}
