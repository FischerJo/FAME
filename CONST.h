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

#ifndef CONST_H
#define CONST_H

#include <cstdint>
#include <array>

namespace MyConst {



//  -------- SET THESE VARIABLES --------
//
//

// maximum read length of the reads in bp
constexpr unsigned int READLEN = 110;

// maximum number of times a k-mer is allowed to occur in the whole genome
constexpr uint64_t KMERCUTOFF = 500;

// Number of cores that this program is allowed to occupy at any given point
#define CORENUM 16


//  --------------------------------------




// --------------------
// ------INTERNAL------
// --------------------

// number of chromosomes in organism
constexpr unsigned int CHROMNUM = 24;

// (more than) size in bp of biggest chromosome in organism
constexpr unsigned int CHROMMAX = 1000000000;

// (more than) number of CpGs in organism
constexpr unsigned int CPGMAX = 100000000;

// number of reads that should be read per batch
// (note that very large values may increase the required amount of RAM drastically)
// recommended is 100 000
constexpr unsigned int CHUNKSIZE = 100000;

// Length of a kmer in bp
// recommended is 25
constexpr unsigned int KMERLEN = 25;

// Kmer bitmask to build reverse complement of 64bit kmer
constexpr uint64_t KMERMASK = (KMERLEN == 32 ? 0xffffffffffffffffULL : (1ULL << 2*KMERLEN) - 1);

// yields value for number of kmers in one read at compile time
// constexpr KPERREAD = READLEN - KMERLEN + 1;


// size of hash table
// recommended is 1 << 30
constexpr uint64_t HTABSIZE = 1ULL << 30;


// window length for meta CpGs
// recommended is 1024
constexpr unsigned int WINLEN = 2048;

// number of mismatches we allow
// recommended is 2
constexpr unsigned int MISCOUNT = 2;

// threshold over which reads are immediately discarded if that many metaCpGs are hit
// constexpr unsigned int QUERYTHRESHOLD = 30000;

// Checks the usefulness of the set parameters
void sanityChecks();


}

#endif /* CONST_H */
