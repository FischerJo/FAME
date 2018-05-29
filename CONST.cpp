//	Metal - A fast methylation alignment and calling tool for WGBS data.
//	Copyright (C) 2017  Jonas Fischer
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	Jonas Fischer	jonaspost@web.de

#include <iostream>

#include "CONST.h"


void MyConst::sanityChecks()
{

    if (MyConst::KMERLEN >= 33)
    {
        std::cerr << "The chosen length for kmers is too large. Maximum value allowed is 32. Please redefine KMERLEN.\n\n";
        exit(1);
    }

	if (MyConst::KMERLEN != MyConst::SEED.size())
	{
		std::cerr << "Seed size is different from k-mer length! Aborting.\n\n";
		exit(1);
	}

    if (MyConst::READLEN <= MyConst::KMERLEN)
    {

        std::cerr << "The chosen maximum read length is smaller then the k (length of a kmer). Please chose another value for either of the two (READLEN, KMERLEN).\n\n";
        exit(1);
    }

    // if (MyConst::READLEN - MyConst::KMERLEN - (MyConst::KMERLEN*MyConst::MISCOUNT))
    // {
    //     std::cerr << "The maximum read length is too small to allow for so many errors. Please redefine the number of allowed mismatches (MISCOUNT).\n\n";
    //     exit(1);
    // }
    //
    if (MyConst::MISCOUNT >= 3)
    {

        std::cout << "The chosen number of allowed mismatches (" << MyConst::MISCOUNT << ") is quite large. Most of the heuristics won't work.\n\n";

    }

    if (MyConst::WINLEN <= (1 << 9))
    {

        std::cout << "The chosen Meta CpG window length (" << MyConst::WINLEN << ") is quite small. You should consider redefining WINLEN (default is 2^12).\n\n";

    }

    if (MyConst::HTABSIZE < 2 << 16)
    {
        std::cout << "The chosen size of your hash table (" << MyConst::HTABSIZE << ") is quite small. You should consider redefining HTABSIZE (default is 2^30 for the human genome).\n\n";
    }
}
