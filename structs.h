#ifndef STRUCTS_H
#define STRUCTS_H

#include <cstdint>
#include <type_traits>

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
    uint8_t chrom;
    // convention: pos points to position of start of context of CpG (C position - READLEN + 2)
    //              or to position of C if the former is negative
    uint32_t pos;
};


struct metaCpG {

    // index of first CpG of this meta CpG
    uint32_t start;
    // start of this meta CpG is cpgTable[start].pos
    // end of this meta CpG is cpgTable[start].pos + WINLEN
    // if startmetaCpG .start is 0 and end is WINLEN

    // index of last CpG in this meta CpG
    uint32_t end;
};


// test for POD types of structs
inline constexpr bool testPODs()
{

    return (std::is_pod<struct CpG>::value && std::is_pod<struct metaCpG>::value);
}

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
// // TODO
// namespace KMER {
//
//     // KMER DEFINITION
//     //
//     // metaID holds index of corresponding meta CpG
//     // offset the offset inside the meta CpG
//     // most significant bit of offset holds flag that is 1 <=> points to start meta CpG
//     struct kmer {
//
//         uint32_t metaID;
//         uint16_t offset;
//     };
//
//     // return the offset of given kmer inside its cpg
//     inline uint64_t getOffset(KMER::kmer& k)
//     {
//         return k.offset & static_cast<uint16_t>(0x7fff);
//     }
//
//     inline uint32_t getMetaCpG(KMER::kmer& k)
//     {
//         return k.metaID;
//     }
//
//     inline bool isStartCpG(KMER::kmer& k)
//     {
//         return k.offset & static_cast<uint16_t>(0x8000);
//     }
//
//     // constructs a kmer according to its definition
//     inline KMER::kmer constructKmer(uint16_t isStart, uint32_t metacpg, uint16_t off)
//     {
//         return  {metacpg, (isStart << 15 | off)};
//     }
// } // end namespace KMER


namespace MATCH {

    // MATCH DEFINITION
    //
    // 16 lower (least significant) bit hold offset where the match ends in the sequence
    // 8 higher bits hold number of errors produced by this match
    // next higher bit is strand flag
    // even next higher bit is start flag
    // 32 most significant bits hold meta CpGID
    typedef uint64_t match;

    // return the offset of given match
    inline uint64_t getOffset(const MATCH::match m)
    {
        return m & 0x000000000000ffffULL;
    }
    inline uint64_t getErrNum(const MATCH::match m)
    {
        return (m & 0x0000000000ff0000ULL) >> 16;
    }
    inline bool isFwd(const MATCH::match m)
    {
        return (m & 0x0000000001000000ULL);
    }
    inline bool isStart(const MATCH::match m)
    {
        return (m & 0x0000000002000000ULL);
    }
    inline uint32_t getMetaID(const MATCH::match m)
    {
        return static_cast<uint32_t>(m >> 32);
    }

    inline MATCH::match constructMatch(uint16_t off, uint8_t errNum, uint64_t isFwd, uint64_t isStart, uint64_t metaID)
    {
        return (static_cast<uint64_t>(off)) | (static_cast<uint64_t>(errNum) << 16) | (isFwd << 24) | (isStart << 25) | (metaID << 32);
    }
} // end namespace MATCH

#endif /* STRUCTS_H */
