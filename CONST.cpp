
#include <iostream>

#include "CONST.h"


void MyConst::sanityChecks()
{

    if (MyConst::KMERLEN >= 33)
    {
        std::cerr << "The chosen length for kmers is too large. Maximum value allowed is 32. Please redefine KMERLEN.\n\n";
        exit(1);
    }

    if (MyConst::READLEN <= MyConst::KMERLEN)
    {

        std::cerr << "The chosen maximum read length is smaller then the k (length of a kmer). Please chose another value for either of the two (READLEN, KMERLEN).\n\n";
        exit(1);
    }

    if (MyConst::READLEN - MyConst::KMERLEN - (MyConst::KMERLEN*MyConst::MISCOUNT))
    {
        std::cerr << "The maximum read length is too small to allow for so many errors. Please redefine the number of allowed mismatches (MISCOUNT).\n\n";
        exit(1);
    }

    if (MyConst::MISCOUNT >= 3)
    {

        std::cout << "The chosen number of allowed mismatches (" << MyConst::MISCOUNT << ") is quite large. Most of the heuristics won't work.\n\n";

    }

    if (MyConst::WINLEN > (1 << 8))
    {

        std::cout << "The chosen Meta CpG window length (" << MyConst::WINLEN << ") is quite small. You should consider redefining WINLEN (default is 2^12).\n\n";

    }

    if (MyConst::HTABSIZE < 2 << 16)
    {
        std::cout << "The chosen size of your hash table (" << MyConst::HTABSIZE << ") is quite small. You should consider redefining HTABSIZE (default is 2^30 for the human genome).\n\n";
    }
}
