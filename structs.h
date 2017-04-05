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


struct kmer {

    // element index of CpG in RefGenome
    // first (most significant) bit is reserved to state if cpg is start cpg (1) or normal (0)
    const uint32_t cpg;
    // offset in this CpG, the C of the CpG position - READLEN + 2 is offset 0
    // 2nd highest bit is flag for strand; use STRANDMASK to extract information
    //      1 -> forward strand
    //      0 -> reverse strand
    const uint16_t offset;
    // Ctor for emplace_back
    kmer(const uint32_t cpgC, uint16_t offsetC) : cpg(cpgC), offset(offsetC) {}
};

struct colList {

    // how many elements are in this collision list
    uint16_t len;

    // structure saving elements
    // Guarantee: collis.size() == len
    std::vector<struct kmer> collis;

    // Ctor for emplace_back
    colList() : len(0), collis() {}
};

// Use this mask to extract strand information from kmer.offset
// offset & STRANDMASK == TRUE   <=>   kmer on forward strand
constexpr uint16_t STRANDMASK = 0x4000;

#endif /* STRUCTS_H */
