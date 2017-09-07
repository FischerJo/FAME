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

#ifndef DNABITSTR_H
#define DNABITSTR_H


#include <vector>
#include <string>
#include <cstdint>

#include <iostream>  // for tests
#include <bitset>   // for tests

#include "BitFunctions.h"

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
//          * -> 11 where * \in {A,G,T,N}
class DnaBitStr
{

    public:


        DnaBitStr() = delete;

        // Ctor
        // Arguments:
        //      size    the length (in bp) of the represented sequence
        DnaBitStr(const unsigned int size);


        // Get n letters of sequence in bit representation starting at pos
        // Returns n/32 many vectors, comprising the sequence.
        // if n%32 != 0 then the last element of return vector will hold n%32 * 2 bits
        // in the less significant bits for the remainder-sequence
        // inline std::vector<uint64_t> getBitStr(const unsigned int pos, const unsigned int n)
        // {
        //
        //     // get position of first part of slice in vector
        //     unsigned int index = pos / 32;
        //     const unsigned int shift = 2*(pos % 32);
        //
        //     // reserve space for output
        //     std::vector<uint64_t> result;
        //     result.reserve(n/32 + 1);
        //
        //     // retrieve sequence
        //     for (unsigned int i = 0; i < n/32; ++i, ++index)
        //     {
        //
        //         result.emplace_back((bitSeq[index] << shift) | (bitSeq[index + 1] >> (64 - shift)));
        //
        //     }
        //     // Retrieve last part of sequence
        //     if (n % 32 > 0)
        //     {
        //
        //         if (2*(n % 32) > shift)
        //         {
        //
        //             result.emplace_back((bitSeq[index] << shift) >> (64 - 2*(n % 32)));
        //
        //         }
        //     }
        //
        //
        //     return result;
        // }


        // Set the n-th 64 bit element of bitSeq and bitMask according to sequence part
        // CONVENTION:  we start counting here at 0 (i.e. n>= 0)
        //              (n + 1)*64 should NOT exceed "size" - 1 (member variable) o/w undefined behaviour
        //              seq should be of length == 32 for n <  \gaussdown size DIV 64 \gaussdown
        // for last element use setBitStrLast
        inline void setBitStrN(std::string&& seq, const unsigned int n)
        {

            uint64_t bitStr = 0;
            uint64_t bitM = 0xffffffffffffffffULL;
            uint64_t bitRevM = 0xffffffffffffffffULL;
            for (unsigned int i = 1; i <= 32; ++i)
            {

                // we don't need A here
                switch (seq[i - 1]){

                    case 'C':
                        {
                            const unsigned int shift = (64 - 2*i);
                            bitStr |= (1ULL << shift);
                            bitM ^= (2ULL << shift);
                        }
                        break;

                    case 'G':
                        {
                            const unsigned int shift = (64 - 2*i);
                            bitStr |= (2ULL << shift);
                            bitRevM ^= (2ULL << shift);
                        }
                        break;

                    case 'T':
                        bitStr |= (3ULL << (64 - 2*i));
                        break;


                    // unknowns, A and N will be encoded as zero; nothing to do
                    default:
                        break;
                }
            }
            bitSeq[n] = bitStr;
            bitMask[n] = bitM;
            bitRevMask[n] = bitRevM;
        }

        void setBitStrLast(std::string&& seq);


        // get a bit slice (kmer) of the sequence's bitstring (or of the reverse complement sequence),
        // starting at pos (pos == 0 is first letter of sequence)
        // returns the kmer as 64bit slice
        // if kmer length < 32 then the last (least significant) kmer length*2 bits will hold the bitstring
        inline uint64_t getSeqKmer(const unsigned int pos)
        {

            // get position of first part of kmer in vector
            const unsigned int  k1 = pos / 32;

            // maximum position that kmer start can have in 64bit word without exceeding the 64 bit
            constexpr unsigned int maxBitPos = 64 - (2 * MyConst::KMERLEN);
            // offset in word
            const unsigned int offBitPos = 2 * (pos % 32);
            if ( offBitPos <= maxBitPos )
            {

                return ((bitSeq[k1] << offBitPos) >> maxBitPos);

            // if necessary get second part of kmer
            } else {

                uint64_t tmp = ((bitSeq[k1] << offBitPos) >> offBitPos) << (offBitPos - maxBitPos);
                // right operand of shift is < 64 so we will be fine
                return tmp | (bitSeq[k1 + 1] >> (64 - (offBitPos - maxBitPos)));
            }
        }

        inline uint64_t getSeqKmerRev(const unsigned int pos)
        {

            // get position of first part of kmer in vector
            const unsigned int  k1 = pos / 32;

            // maximum position that kmer start can have in 64bit word without exceeding the 64 bit
            constexpr unsigned int maxBitPos = 64 - (2 * MyConst::KMERLEN);
            // offset in word
            const unsigned int offBitPos = 2 * (pos % 32);
            if ( offBitPos <= maxBitPos )
            {

                uint64_t tmp = ((bitSeq[k1] << offBitPos) >> maxBitPos);
                return BitFun::revKmer(tmp);

            // if necessary get second part of kmer
            } else {

                uint64_t tmp = ((bitSeq[k1] << offBitPos) >> offBitPos) << (offBitPos - maxBitPos);
                // right operand of shift is < 64 so we will be fine
                tmp = tmp | (bitSeq[k1 + 1] >> (64 - (offBitPos - maxBitPos)));
                return BitFun::revKmer(tmp);
            }
        }



        // return bitmask slice of sequence (or of the reverse complement sequence)
        // see getSeqKmer for more info
        inline uint64_t getMaskKmer(const unsigned int pos)
        {

            // get position of first part of kmer in vector
            const unsigned int  k1 = pos / 32;


            // maximum position that kmer start can have in 64bit word without exceeding the 64 bit
            constexpr unsigned int maxBitPos = 64 - 2*MyConst::KMERLEN;
            // offset in word
            const unsigned int offBitPos = 2 * (pos % 32);

            if ( offBitPos <= maxBitPos )
            {

                return ((bitMask[k1] << offBitPos) >> maxBitPos);

            // if necessary get second part of kmer
            } else {

                uint64_t tmp = ((bitMask[k1] << offBitPos) >> offBitPos) << (offBitPos - maxBitPos);
                // right operand of shift is < 64 so we will be fine
                return tmp | (bitMask[k1 + 1] >> (64 - (offBitPos - maxBitPos)));
            }
        }
        uint64_t getMaskKmerRev(const unsigned int pos)
        {

            // get position of first part of kmer in vector
            const unsigned int  k1 = pos / 32;

            // maximum position that kmer start can have in 64bit word without exceeding the 64 bit
            constexpr unsigned int maxBitPos = 64 - 2*MyConst::KMERLEN;
            // offset in word
            const unsigned int offBitPos = 2 * (pos % 32);
            if ( offBitPos <= maxBitPos )
            {

                return BitFun::rev64((bitRevMask[k1] << offBitPos) >> maxBitPos) >> maxBitPos;

            // if necessary get second part of kmer
            } else {

                uint64_t tmp = ((bitRevMask[k1] << offBitPos) >> offBitPos) << (offBitPos - maxBitPos);
                // right operand of shift is < 64 so we will be fine
                tmp |= (bitRevMask[k1 + 1] >> (64 - (offBitPos - maxBitPos)));
                return BitFun::rev64(tmp) >> maxBitPos;
            }
        }

    private:

        // length of the represented sequence
        const unsigned int size;

        // represented sequence as bitstring
        // last 64 might not be filled completely
        std::vector<uint64_t> bitSeq;

        // corresponding bitmask
        std::vector<uint64_t> bitMask;
        // corresponding reverse complement bitmask
        // in order of the FORWARD strand, so reverse bitmask before use
        std::vector<uint64_t> bitRevMask;

};

#endif /* DNABITSTR_H */
