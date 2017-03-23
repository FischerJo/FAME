#ifndef STRUCTS_H
#define STRUCTS_H

#include <cstdint>

// struct that holds information about id tags in buffer
// id tags start with '>'
struct idPos {

    // pointer to id line
    const char* id;
    // pointer to next line after the id
    const char* sec;
    // flag if tag is primary assembly
    const bool imp;
};


struct CpG {
    // position where this CpG is present in genome
    // convention: pos points to the location of the C
    const uint8_t chrom;
    const unsigned int pos;
};


struct kmer {
    // pointer to related CpG
    const struct CpG * cpg;
    // offset in this CpG, the C of the CpG is offset 0
    // 2nd highest bit is flag for strand; use STRANDMASK to extract information
    //      1 -> forward strand
    //      0 -> reverse strand
    const uint16_t offset;

};

// Use this mask to extract strand information from kmer.offset
// offset & STRANDMASK == TRUE   <=>   kmer on forward strand
constexpr STRANDMASK = 0x4000;

#endif /* STRUCTS_H */
