#ifndef DNABITSTR_H
#define DNABITSTR_H





// bit string representation of the DNA optimized for slice access
// and Methylation procedures using an additional bitmask
// Encoding for strand:
//          A -> 00
//          C -> 01
//          G -> 10
//          T -> 11
//          N -> 00   DISCARDED IN ANY COMPUTATIONS
//
// Bitmask:
//          C -> 01
//          * -> 00 where * \in {A,G,T,N}
class DnaBitStr
{

    public:


        DnaBitStr() = delete;

        // Ctor
        // Arguments:
        //      size    the length (in bp) of the represented sequence
        DnaBitStr(const unsigned int size);


        // Compute the bit string representation and bitmasks of seq
        // IMPORTANT: seq should be of length "size" (member variable) o/w undefined behaviour
        void computeBitStr(std::string& seq);


        // Set the n-th 64 bit element of bitSeq and bitMask according to sequence part
        // CONVENTION:  n*64 should NOT exceed "size" (member variable) o/w undefined behaviour
        //              seq should be of length == 32 for n <  \gaussup size DIV 64 \gaussup
        //                               length <= 32 for n == \gaussup offset DIV 64 \gaussup
        void setBitStrN(std::string& seq, const unsigned int n)


        // get a bit slice of the sequence's bitstring,
        // offset many characters as bit representation starting at pos
        // returns the slice as a vector of 64bit slices
        // length of the vector is \gaussup offset DIV 64 \gaussup
        //
        // CONVENTION: The first character is pos = 0
        //             The value of the last 64 - (offset MOD 64) bits in the last 64bit element are undefined
        std::vector<uint64_t> getSeqSlice(const unsigned int pos, const unsigned int offset);


        // return bitmask slice
        // see getSeqSlice for more info
        std::vector<uint64_t> getMaskSlice(const unsigned int pos, const unsigned int offset);

    private:

        // length of the represented sequence
        const unsigned int size;

        // represented sequence as bitstring
        // last 64 might not be filled completely
        std::vector<uint64_t> bitSeq;

        // corresponding bitmask
        std::vector<uint64_t> bitMask;

};

#endif /* DNABITSTR_H */
