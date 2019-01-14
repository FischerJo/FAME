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
#include <vector>
#include <array>

namespace MyConst {



//  -------- SET THESE VARIABLES --------
//
//

// maximum read length of the reads in bp
constexpr unsigned int READLEN = 101;

// Number of cores that this program is allowed to occupy at any given point
#define CORENUM 32

// closed interval borders for distances allowed between paired reads
constexpr uint32_t MINPDIST = 50;
constexpr uint32_t MAXPDIST = 450;

// number of chromosomes in organism
constexpr unsigned int CHROMNUM = 24;

//  --------------------------------------




// --------------------
// ------INTERNAL------
// --------------------

// seed used for spaced k-mer hashing
const std::vector<bool> SEED = {1,1,1,1,0,1,1,1,0,1,1,0,1,1,1,1,0,0,1,0,1,1,1,0,1,1,1,0,1,1,1,1};
constexpr uint32_t SEEDBITS = 0b11110111011011110010111011101111;

// (more than) size in bp of biggest chromosome in organism
constexpr unsigned int CHROMMAX = 1000000000;

// (more than) number of CpGs in organism
constexpr unsigned int CPGMAX = 100000000;

// number of reads that should be read per batch
// (note that very large values may increase the required amount of RAM drastically)
// recommended is 300 000
constexpr unsigned int CHUNKSIZE = 300000;

// Length of a kmer in bp
// WARNING: THIS PARAMETER MUST EQUAL SEED.size() !!!
constexpr unsigned int KMERLEN = 32;

// Kmer bitmask
constexpr uint32_t KMERMASK = (KMERLEN == 32 ? 0xffffffff : ((uint64_t)1 << KMERLEN) - 1);

// minimum number of k-mers required to test for match
// recommended is 10
constexpr uint16_t QTHRESH = 10;


// size of hash table
// recommended is 1 << 30
constexpr uint64_t HTABSIZE = 1ULL << 30;


// window length for meta CpGs
// recommended is 2048
constexpr unsigned int WINLEN = 2048;

// number of mismatches we allow initially (can be extended by ADDMIS constant)
// recommended is 2
constexpr uint8_t MISCOUNT = 2;
// number of mismatches we allow additionally for shift-and and alignment
// recommended is 4
constexpr uint8_t ADDMIS = 4;

// maximum number of times a k-mer is allowed to occur in the whole genome
// recommended is 1500
constexpr uint64_t KMERCUTOFF = 1500;

// Checks the usefulness of the set parameters
void sanityChecks();


}

#endif /* CONST_H */
