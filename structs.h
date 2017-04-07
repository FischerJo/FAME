#ifndef STRUCTS_H
#define STRUCTS_H

#include <cstdint>

// struct that holds information about id tags in buffer
// id tags start with '>'
struct idPos {

    // pointer to id line
    char* const id;
    // pointer to next line after the id
    char* const sec;
    // flag if tag is primary assembly
    const bool imp;
    // Ctor for emplace_back
    idPos(char* const idC, char* const secC, const bool impC) : id(idC), sec(secC), imp(impC) {}
};


struct CpG {
    // position where this CpG is present in genome
    const uint8_t chrom;
    // convention: pos points to position of start of context of CpG (C position - READLEN + 2)
    //              or to position of C if the former is negative
    const uint32_t pos;

    // Ctor for emplace_back
    CpG(uint8_t chromC, unsigned int posC) : chrom(chromC), pos(posC) {}
};


namespace KMER {

    // KMER DEFINITION
    //
    // most signficant bit hold if cpg is start cpg (1) or normal (0)
    // next 31 bit state if the index of the kmer in the cpg vector
    // 33th most significant bit states if kmer is of forward (1) or reverse strand (0)
    // remaining bits define the offset (starting at C position - READLEN + 2 as 0) in the CpG context
    typedef uint64_t kmer;

    // Use this mask to extract strand information from kmer.offset
    // offset & STRANDMASK == TRUE   <=>   kmer on forward strand
    constexpr uint64_t STRANDMASK = 0x0000000080000000ULL;
    // cpg mask to identify start (1) and normal (0) cpgs
    constexpr uint64_t CPGMASK = 0x8000000000000000ULL;

    // returns 1 if kmer lies on forwards strand, 0 on backward
    inline bool isForward(KMER::kmer& k)
    {
        return (k & KMER::STRANDMASK) >> 31;
    }

    // returns whether kmer is of start cpg (1) or normal (0)
    inline bool isStartCpG(KMER::kmer& k)
    {
        return (k & KMER::CPGMASK) >> 63;
    }

    // return the offset of given kmer inside its cpg
    inline uint64_t getOffset(KMER::kmer& k)
    {
        return k & 0x000000007fffffffULL;
    }

    inline uint64_t getCpG(KMER::kmer& k)
    {
        return (k & ~KMER::CPGMASK) >> 32;
    }

    // constructs a kmer according to its definition
    inline KMER::kmer constructKmer(uint64_t strand, uint64_t CpG, uint64_t off, uint64_t isStart)
    {

        return (isStart << 63) | (CpG << 32) | (strand << 31) | off;
    }
}

// struct kmer {
//
//     // element index of CpG in RefGenome
//     // first (most significant) bit is reserved to state if cpg is start cpg (1) or normal (0)
//     const uint32_t cpg;
//     // offset in this CpG, the C of the CpG position - READLEN + 2 is offset 0
//     // 2nd highest bit is flag for strand; use STRANDMASK to extract information
//     //      1 -> forward strand
//     //      0 -> reverse strand
//     const uint16_t offset;
//     // Ctor for emplace_back
//     kmer(const uint32_t cpgC, uint16_t offsetC) : cpg(cpgC), offset(offsetC) {}
// };

// struct colList {
//
//     // how many elements are in this collision list
//     uint16_t len;
//
//     // structure saving elements
//     // Guarantee: collis.size() == len
//     std::vector<struct kmer> collis;
//
//     // Ctor for emplace_back
//     colList() : len(0), collis() {}
// };


#endif /* STRUCTS_H */
