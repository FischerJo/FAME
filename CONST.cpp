
#include <iostream>

#include "CONST.h"


void MyConst::sanityChecks()
{

    if (MyConst::READLEN <= MyConst::KMERLEN)
    {

        std::cerr << "The chosen maximum read length is smaller then the k (length of a kmer). Please chose another value for either of the two.\n";
        exit(1);

    }

    if (MyConst::MISCOUNT >= 5)
    {

        std::cout << "The chosen number of allowed mismatches (" << MyConst::MISCOUNT << ") is quite large.\n";

    }

    if (MyConst::WINLEN > (1 << 8))
    {

        std::cout << "The chosen Meta CpG window length (" << MyConst::WINLEN << ") is quite small. You should consider redefining it (default is 2^14).\n";

    }

}
