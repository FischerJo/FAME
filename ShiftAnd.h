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

#ifndef SHIFTAND_H
#define SHIFTAND_H

#include <vector>
#include <array>
#include <cstdint>
#include <iostream>

#include "CONST.h"

// Implementation of the approximate shift and algorithm with maximum E errors allowed
template<size_t E>
class ShiftAnd
{

    // Internal representation of states
    // one layer of states encoded by 128 bits (Note that our shift-and automata
    // therefore represents at most 128-1 letters, first state is always active dummy state)
    struct states {

        uint64_t B_0;
        uint64_t B_1;

    };
    using bitMasks = struct states;
    using bitStates = struct states;

    public:

        // C-tor -------------

        ShiftAnd() = delete;

        // ARGUMENTS:
        //              seq     sequence of characters for which the bitmasks should be
        //                      initialized
        //              lMap    array that maps letters to mask indices
        ShiftAnd(std::string& seq, std::array<uint8_t, 16>& lMap);

        // -------------------

        // Query the sequence slice framed by [start,end) to the internal automata
        // Query forward sequence to queryRevSeq to start a query of the reverse complement of the argument to the automaton
        //
        // ARGUMENTS:
        //              start       iterator to start of sequence to query / reverse seq: iterator to last letter of seq to query
        //              end         iterator to end of sequence to query / reverse seq: iterator one before first letter of seq to query
        //              matches     will contain the matchings as offset relative to the start iterator (the end of the match)
        //              errors      same size as matches; will contain number of errors for each match
        //
        inline void querySeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<uint64_t>& matches, std::vector<uint8_t>& errors);
        inline void queryRevSeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<uint64_t>& matches, std::vector<uint8_t>& errors);

        // returns the size of the represented pattern sequence
        inline uint64_t size() { return pLen; }


        // TODO make this private once its tested
    // private:

        // Reset to initial state
        // is called by querySeq
        inline void reset();

        // query a single letter to the automaton
        inline void queryLetter(const char& c);

		// get index for letter without overhead in memory
		// ASSUMING ASCII ENCODING!!!
		// inline uint8_t getCharID(const char c)
		// {
		// 	return ((c % 16) >> 1) ^ 1;
		// }

        // test if we have reached an accepting state
        //
        // ARGUMENTS:
        //          errNum      will contain number of errors produced by current match
        //
        // RETURN:
        //          true iff we have a match
        inline bool isMatch(uint8_t& errNum);

        // load the bitmasks for the given sequences
        inline void loadBitmasks(std::string& seq);

        // Bitvector data structure holding the set of active states indicated by a 1
        // initial state is least significant bit of active[0].B_0 (dummy state),
        // different layers indicating the number of errors are indexed by the vector
        std::array<bitStates, E + 1> active;

        // Set of bitmasks for the given sequence
        // Generally, Masks[m] give the bitmask for letter m
        std::array<bitMasks, 4> masks;
        // masks to filter out accepting states
        bitMasks accepted;

        // length of the sequence that is represented by the automata
        const uint64_t pLen;


        // maps characters 'A', 'C', 'G', 'T' to their index
        std::array<uint8_t, 16>& lmap;


};


template<size_t E>
ShiftAnd<E>::ShiftAnd(std::string& seq, std::array<uint8_t, 16>& lMap) :
        pLen(seq.size())
    ,   lmap(lMap)
{
    loadBitmasks(seq);
}


template<size_t E>
inline void ShiftAnd<E>::querySeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<uint64_t>& matches, std::vector<uint8_t>& errors)
{

    reset();

    bool wasMatch = false;
    uint8_t prevErrs = MyConst::MISCOUNT + 1;
    size_t numCompLets = 0;
    for (auto it = start; it < end; ++it)
    {

        // we do not consider Ns for matches - restart whole automaton for next letter
        if (*it == 'N')
        {
            reset();
            continue;
        }
        queryLetter(*it);

        ++numCompLets;

        // There can only be a match after at least this->size() - E many chars are queried so only then compare
        if (numCompLets < (pLen - E))
            continue;

        uint8_t errNum;
        if (isMatch(errNum))
        {

            // if we matched in previous round, overwrite that match
            if (wasMatch)
            {
                if (errNum < prevErrs)
                {
                    matches.back() = it - start;
                    errors.back() = errNum;
                    prevErrs = errNum;
                }

            } else {

                matches.emplace_back(it - start);
                errors.emplace_back(errNum);
                wasMatch = true;
                prevErrs = errNum;
            }

        } else {

            wasMatch = false;
        }
    }
}


template<size_t E>
inline void ShiftAnd<E>::queryRevSeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<uint64_t>& matches, std::vector<uint8_t>& errors)
{

    reset();

    bool wasMatch = false;
    uint8_t prevErrs = MyConst::MISCOUNT + 1;
    size_t numCompLets = 0;
    for (auto it = start; it != end; --it)
    {

        char c;
        switch (*it)
        {
            case 'A':
                c = 'T';
                break;
            case 'C':
                c = 'G';
                break;
            case 'G':
                c = 'C';
                break;
            case 'T':
                c = 'A';
                break;
            // we do not consider Ns for matches - restart whole automaton for next letter
            case 'N':
                reset();
                continue;
            default:
                std::cerr << "[ShiftAnd] Could not parse letter " << *it << " at offset " << it - end + 1 << "\nStopping program...\n\n";
                exit(1);
        }
        queryLetter(c);

        ++numCompLets;

        // There can only be a match after at least this->size() - E many chars are queried so only then compare
        if (numCompLets < (pLen - E))
            continue;

        uint8_t errNum;
        if (isMatch(errNum))
        {

            // if we matched in previous round, overwrite that match
            if (wasMatch)
            {
                if (errNum < prevErrs)
                {
                    matches.back() = it - end + pLen - 2;
                    errors.back() = errNum;
                    prevErrs = errNum;
                }

            } else {

                matches.emplace_back(it - end + pLen - 2);
                errors.emplace_back(errNum);
                wasMatch = true;
                prevErrs = errNum;
            }

        } else {

            wasMatch = false;
        }
    }
}


template<size_t E>
inline void ShiftAnd<E>::reset()
{
    for (size_t i = 0; i <= E; ++i)
    {

        // set all initial states (with respect to epsilon transitions) active
        active[i].B_0 = (static_cast<uint64_t>(1) << (i+1)) - 1;
        active[i].B_1 = 0;

    }
}
template <>
inline void ShiftAnd<0>::reset()
{
    active[0].B_0 = 1;
    active[0].B_1 = 0;

}
template <>
inline void ShiftAnd<1>::reset()
{
    active[0].B_0 = 1;
    active[0].B_1 = 0;
    active[1].B_0 = 3;
    active[1].B_1 = 0;

}
template <>
inline void ShiftAnd<2>::reset()
{
    active[0].B_0 = 1;
    active[0].B_1 = 0;
    active[1].B_0 = 3;
    active[1].B_1 = 0;
    active[2].B_0 = 7;
    active[2].B_1 = 0;

}
template <>
inline void ShiftAnd<3>::reset()
{
    active[0].B_0 = 1;
    active[0].B_1 = 0;
    active[1].B_0 = 3;
    active[1].B_1 = 0;
    active[2].B_0 = 7;
    active[2].B_1 = 0;
    active[3].B_0 = 15;
    active[3].B_1 = 0;
}
template <>
inline void ShiftAnd<4>::reset()
{
    active[0].B_0 = 1;
    active[0].B_1 = 0;
    active[1].B_0 = 3;
    active[1].B_1 = 0;
    active[2].B_0 = 7;
    active[2].B_1 = 0;
    active[3].B_0 = 15;
    active[3].B_1 = 0;
    active[4].B_0 = 31;
    active[4].B_1 = 0;
}
template <>
inline void ShiftAnd<5>::reset()
{
    active[0].B_0 = 1;
    active[0].B_1 = 0;
    active[1].B_0 = 3;
    active[1].B_1 = 0;
    active[2].B_0 = 7;
    active[2].B_1 = 0;
    active[3].B_0 = 15;
    active[3].B_1 = 0;
    active[4].B_0 = 31;
    active[4].B_1 = 0;
    active[5].B_0 = 63;
    active[5].B_1 = 0;
}
template <>
inline void ShiftAnd<6>::reset()
{
    active[0].B_0 = 1;
    active[0].B_1 = 0;
    active[1].B_0 = 3;
    active[1].B_1 = 0;
    active[2].B_0 = 7;
    active[2].B_1 = 0;
    active[3].B_0 = 15;
    active[3].B_1 = 0;
    active[4].B_0 = 31;
    active[4].B_1 = 0;
    active[5].B_0 = 63;
    active[5].B_1 = 0;
    active[6].B_0 = 127;
    active[6].B_1 = 0;
}


template<size_t E>
inline void ShiftAnd<E>::queryLetter(const char& c)
{

    const bitMasks& mask = masks[lmap[c%16]];

    // Bottom up update part for old values of previous iteration
    for (size_t i = E; i > 0; --i)
    {

        // Update second part of pattern states
        //
        //                                  Match                                       Insertion                       Substitution
        active[i].B_1 = ((active[i].B_1 << 1 | active[i].B_0 >> 63) & mask.B_1) | (active[i-1].B_1) | (active[i-1].B_1 << 1 | active[i-1].B_0 >> 63);

        // Update first part of pattern states
        //
        //                          Match                           Insertion               Substitution
        active[i].B_0 = ((active[i].B_0 << 1 | 1) & mask.B_0) | (active[i-1].B_0) | (active[i-1].B_0 << 1);

    }
    // update zero error layer (at the top)
    active[0].B_1 = ((active[0].B_1 << 1 | active[0].B_0 >> 63) & mask.B_1);
    active[0].B_0 = ((active[0].B_0 << 1 | 1) & mask.B_0);
    //
    // Top down update for values of this iteration
    for (size_t i = 1; i <= E; ++i)
    {

        active[i].B_1 |= active[i-1].B_1 << 1 | active[i-1].B_0 >> 63;
        active[i].B_0 |= active[i-1].B_0 << 1;

    }
}
// template spec for 0 errors
template<>
inline void ShiftAnd<0>::queryLetter(const char& c)
{

    const bitMasks& mask = masks[lmap[c%16]];

    // update zero error layer (at the top)
    active[0].B_1 = ((active[0].B_1 << 1 | active[0].B_0 >> 63) & mask.B_1);
    active[0].B_0 = ((active[0].B_0 << 1 | 1) & mask.B_0);
}
// template spec for 1 error
template<>
inline void ShiftAnd<1>::queryLetter(const char& c)
{

    const bitMasks& mask = masks[lmap[c%16]];

    // Bottom up update part for old values of previous iteration for bottom layer
    // Update second part of pattern states
    //
    //                                  Match                                       Insertion                       Substitution
    active[1].B_1 = ((active[1].B_1 << 1 | active[1].B_0 >> 63) & mask.B_1) | (active[0].B_1) | (active[0].B_1 << 1 | active[0].B_0 >> 63);

    // Update first part of pattern states
    //
    //                          Match                           Insertion               Substitution
    active[1].B_0 = ((active[1].B_0 << 1 | 1) & mask.B_0) | (active[0].B_0) | (active[0].B_0 << 1);

    // update zero error layer (at the top)
    active[0].B_1 = ((active[0].B_1 << 1 | active[0].B_0 >> 63) & mask.B_1);
    active[0].B_0 = ((active[0].B_0 << 1 | 1) & mask.B_0);
    //
    // Top down update for values of this iteration
    active[1].B_1 |= active[0].B_1 << 1 | active[0].B_0 >> 63;
    active[1].B_0 |= active[0].B_0 << 1;
}
// template spec for 2 errors
template<>
inline void ShiftAnd<2>::queryLetter(const char& c)
{

    const bitMasks& mask = masks[lmap[c%16]];

    // Bottom up update part for old values of previous iteration for bottom layer
    // Update second part of pattern states
    //
    //                                  Match                                       Insertion                       Substitution
    active[2].B_1 = ((active[2].B_1 << 1 | active[2].B_0 >> 63) & mask.B_1) | (active[1].B_1) | (active[1].B_1 << 1 | active[1].B_0 >> 63);

    // Update first part of pattern states
    //
    //                          Match                           Insertion               Substitution
    active[2].B_0 = ((active[2].B_0 << 1 | 1) & mask.B_0) | (active[1].B_0) | (active[1].B_0 << 1);
    // Update second part of pattern states
    //
    //                                  Match                                       Insertion                       Substitution
    active[1].B_1 = ((active[1].B_1 << 1 | active[1].B_0 >> 63) & mask.B_1) | (active[0].B_1) | (active[0].B_1 << 1 | active[0].B_0 >> 63);

    // Update first part of pattern states
    //
    //                          Match                           Insertion               Substitution
    active[1].B_0 = ((active[1].B_0 << 1 | 1) & mask.B_0) | (active[0].B_0) | (active[0].B_0 << 1);

    // update zero error layer (at the top)
    active[0].B_1 = ((active[0].B_1 << 1 | active[0].B_0 >> 63) & mask.B_1);
    active[0].B_0 = ((active[0].B_0 << 1 | 1) & mask.B_0);
    //
    // Top down update for values of this iteration
    active[1].B_1 |= active[0].B_1 << 1 | active[0].B_0 >> 63;
    active[1].B_0 |= active[0].B_0 << 1;
    active[2].B_1 |= active[1].B_1 << 1 | active[1].B_0 >> 63;
    active[2].B_0 |= active[1].B_0 << 1;
}



template<size_t E>
inline bool ShiftAnd<E>::isMatch(uint8_t& errNum)
{
    // go through layers
    for (size_t i = 0; i <= E; ++i)
    {

        // test if this layer is match
        if ((active[i].B_1 & accepted.B_1) || (active[i].B_0 & accepted.B_0))
        {
            errNum = i;
            return true;
        }

    }
    return false;
}
// template spec for 0 errors
template<>
inline bool ShiftAnd<0>::isMatch(uint8_t& errNum)
{

    if ((active[0].B_1 & accepted.B_1) || (active[0].B_0 & accepted.B_0))
    {
        errNum = 0;
        return true;
    }
    return false;
}
// template spec for 1 error
template<>
inline bool ShiftAnd<1>::isMatch(uint8_t& errNum)
{

    if ((active[0].B_1 & accepted.B_1) || (active[0].B_0 & accepted.B_0) )
    {
        errNum = 0;
        return true;
    }
    if ((active[1].B_1 & accepted.B_1) || (active[1].B_0 & accepted.B_0) )
    {
        errNum = 1;
        return true;
    }
    return false;
}
// template spec for 2 errors
template<>
inline bool ShiftAnd<2>::isMatch(uint8_t& errNum)
{

    if ((active[0].B_1 & accepted.B_1) || (active[0].B_0 & accepted.B_0) )
    {
        errNum = 0;
        return true;
    }
    if ((active[1].B_1 & accepted.B_1) || (active[1].B_0 & accepted.B_0) )
    {
        errNum = 1;
        return true;
    }
    if ((active[2].B_1 & accepted.B_1) || (active[2].B_0 & accepted.B_0) )
    {
        errNum = 2;
        return true;
    }
    return false;
}


template<size_t E>
inline void ShiftAnd<E>::loadBitmasks(std::string& seq)
{

    // retrieve length of the sequence
    const size_t seqLen = seq.size();

    // load full mask for trailing 1s for sequences that are shorter than pLen
    uint64_t maskA = 0xffffffffffffffffULL;
    uint64_t maskC = 0xffffffffffffffffULL;
    uint64_t maskG = 0xffffffffffffffffULL;
    uint64_t maskT = 0xffffffffffffffffULL;

    if (seqLen < 64)
    {

        accepted.B_1 = 0;
        accepted.B_0 = static_cast<uint64_t>(1) << seqLen;
        // save masks for second part of sequence
        // dummy bitmask to propagate states to the end
        masks[0].B_1 = maskA;
        masks[1].B_1 = maskA;
        masks[2].B_1 = maskA;
        masks[3].B_1 = maskA;

        for (size_t i = seqLen; i > 0; )
        {

            // shift to make space for next letter
            maskA = maskA << 1;
            maskC = maskC << 1;
            maskG = maskG << 1;
            maskT = maskT << 1;

            switch (seq[--i])
            {
                case ('A'):

                    maskA |= 1;
                    break;

                case ('C'):

                    maskC |= 1;
                    break;

                case ('G'):

                    maskG |= 1;
                    break;

                case ('T'):

                    maskT |= 1;
                    maskC |= 1;

            }

        }
        // load initial state
        maskA = (maskA << 1) | 1;
        maskC = (maskC << 1) | 1;
        maskG = (maskG << 1) | 1;
        maskT = (maskT << 1) | 1;
        // save masks
        masks[lmap['A'%16]].B_0 = maskA;
        masks[lmap['C'%16]].B_0 = maskC;
        masks[lmap['G'%16]].B_0 = maskG;
        masks[lmap['T'%16]].B_0 = maskT;



    } else if (seqLen <= 127) {

        accepted.B_0 = 0;
        accepted.B_1 = static_cast<uint64_t>(1) << (seqLen - 64);

        // load second part of sequence
        for (size_t i = seqLen - 1; i > 62; --i)
        {

            // shift to make space for next letter
            maskA = maskA << 1;
            maskC = maskC << 1;
            maskG = maskG << 1;
            maskT = maskT << 1;

            switch (seq[i])
            {
                case ('A'):

                    maskA |= 1;
                    break;

                case ('C'):

                    maskC |= 1;
                    break;

                case ('G'):

                    maskG |= 1;
                    break;

                case ('T'):

                    maskT |= 1;
                    maskC |= 1;

            }

        }
        // save masks
        masks[lmap['A'%16]].B_1 = maskA;
        masks[lmap['C'%16]].B_1 = maskC;
        masks[lmap['G'%16]].B_1 = maskG;
        masks[lmap['T'%16]].B_1 = maskT;

        // reset mask buffer
        maskA = 0xffffffffffffffffULL;
        maskC = 0xffffffffffffffffULL;
        maskG = 0xffffffffffffffffULL;
        maskT = 0xffffffffffffffffULL;

        // load first part of sequence
        for (size_t i = 63; i > 0; )
        {
            // shift to make space for next letter
            maskA = maskA << 1;
            maskC = maskC << 1;
            maskG = maskG << 1;
            maskT = maskT << 1;

            switch (seq[--i])
            {
                case ('A'):

                    maskA |= 1;
                    break;

                case ('C'):

                    maskC |= 1;
                    break;

                case ('G'):

                    maskG |= 1;
                    break;

                case ('T'):

                    maskT |= 1;
                    maskC |= 1;

            }

        }
        // load initial state
        maskA = (maskA << 1) | 1;
        maskC = (maskC << 1) | 1;
        maskG = (maskG << 1) | 1;
        maskT = (maskT << 1) | 1;
        // save masks
        masks[lmap['A'%16]].B_0 = maskA;
        masks[lmap['C'%16]].B_0 = maskC;
        masks[lmap['G'%16]].B_0 = maskG;
        masks[lmap['T'%16]].B_0 = maskT;

    } else {

        accepted.B_0 = 0;
        accepted.B_1 = 0x8000000000000000ULL;
        // load second part of sequence
        for (size_t i = 127; i > 62; --i)
        {

            // shift to make space for next letter
            maskA = maskA << 1;
            maskC = maskC << 1;
            maskG = maskG << 1;
            maskT = maskT << 1;

            switch (seq[i])
            {
                case ('A'):

                    maskA |= 1;
                    break;

                case ('C'):

                    maskC |= 1;
                    break;

                case ('G'):

                    maskG |= 1;
                    break;

                case ('T'):

                    maskT |= 1;
                    maskC |= 1;

            }

        }
        // save masks
        masks[lmap['A'%16]].B_1 = maskA;
        masks[lmap['C'%16]].B_1 = maskC;
        masks[lmap['G'%16]].B_1 = maskG;
        masks[lmap['T'%16]].B_1 = maskT;

        // reset mask buffer
        maskA = 0xffffffffffffffffULL;
        maskC = 0xffffffffffffffffULL;
        maskG = 0xffffffffffffffffULL;
        maskT = 0xffffffffffffffffULL;

        // load first part of sequence
        for (size_t i = 63; i > 0; )
        {

            // shift to make space for next letter
            maskA = maskA << 1;
            maskC = maskC << 1;
            maskG = maskG << 1;
            maskT = maskT << 1;
            switch (seq[--i])
            {
                case ('A'):

                    maskA |= 1;
                    break;

                case ('C'):

                    maskC |= 1;
                    break;

                case ('G'):

                    maskG |= 1;
                    break;

                case ('T'):

                    maskT |= 1;
                    maskC |= 1;

            }

        }
        // load initial state
        maskA = (maskA << 1) | 1;
        maskC = (maskC << 1) | 1;
        maskG = (maskG << 1) | 1;
        maskT = (maskT << 1) | 1;
        // save masks
        masks[lmap['A'%16]].B_0 = maskA;
        masks[lmap['C'%16]].B_0 = maskC;
        masks[lmap['G'%16]].B_0 = maskG;
        masks[lmap['T'%16]].B_0 = maskT;
    }

}


#endif /* SHIFTAND_H */
