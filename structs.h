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
    // convention: pos points to the location of the C
    const uint8_t chrom;
    const unsigned int pos;
    // Ctor for emplace_back
    CpG(uint8_t chromC, unsigned int posC) : chrom(chromC), pos(posC) {}
};


struct kmer {

    // pointer to related CpG
    const struct CpG * cpg;
    // offset in this CpG, the C of the CpG position - READLEN + 2 is offset 0
    // 2nd highest bit is flag for strand; use STRANDMASK to extract information
    //      1 -> forward strand
    //      0 -> reverse strand
    const uint16_t offset;
    // Ctor for emplace_back
    kmer(const struct CpG* cpgC, uint16_t offsetC) : cpg(cpgC), offset(offsetC) {}
};

// Use this mask to extract strand information from kmer.offset
// offset & STRANDMASK == TRUE   <=>   kmer on forward strand
constexpr uint16_t STRANDMASK = 0x4000;

#endif /* STRUCTS_H */
