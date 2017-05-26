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


struct metaCpG {
    // index of first CpG of this meta CpG
    uint32_t start;
    // start of this meta CpG is cpgTable[start].pos
    // end of this meta CpG is cpgTable[start].pos + WINLEN
    // if startmetaCpG .start is 0 and end is WINLEN

    // index of last CpG in this meta CpG
    uint32_t end;

    // Ctor for emplace_back
    metaCpG(uint32_t startC, uint32_t endC) : start(startC), end(endC) {}
};


namespace KMER {

    // KMER DEFINITION
    //
    // 32 upper (more significant) bit hold metaCpG index, lower 32 bit hold offset inside metaCpG
    // most significant bit holds flag that is 1 <=> points to start meta CpG
    typedef uint64_t kmer;


    // return the offset of given kmer inside its cpg
    inline uint64_t getOffset(KMER::kmer& k)
    {
        return k & 0x00000000ffffffffULL;
    }

    inline uint64_t getMetaCpG(KMER::kmer& k)
    {
        return (k & 0x7fffffff00000000ULL) >> 32;
    }

    inline bool isStartCpG(KMER::kmer& k)
    {
        return (k & 0x8000000000000000ULL);
    }

    // constructs a kmer according to its definition
    inline KMER::kmer constructKmer(uint64_t isStart, uint64_t metacpg, uint64_t off)
    {

        return  (isStart << 63) | (metacpg << 32) | off;
    }
} // end namespace KMER



#endif /* STRUCTS_H */
