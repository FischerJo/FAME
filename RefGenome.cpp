

#include "RefGenome.H"


using namespace std;

RefGenome::RefGenome(unsigned int cpgs, const vector<const struct CpG>&& cpgTab,
        const vector<const string>&& genSeq) :
    ,   cpgNum(cpgs)
    ,   cpgTable(cpgTab)
    ,   genomeSeq(genSeq)
{
}
