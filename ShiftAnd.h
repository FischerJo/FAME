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
        ShiftAnd(std::string& seq, std::array<uint8_t, 256>& lMap);

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
        inline void querySeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<size_t>& matches, std::vector<uint16_t>& errors);
        inline void queryRevSeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<uint64_t>& matches, std::vector<uint16_t>& errors);

        // returns the size of the represented pattern sequence
        inline uint64_t size() { return pLen; }


        // TODO make this private once its tested
    // private:

        // Reset to initial state
        // is called by querySeq
        inline void reset();

        // query a single letter to the automaton
        inline void queryLetter(const char& c);

        // test if we have reached an accepting state
        //
        // ARGUMENTS:
        //          errNum      will contain number of errors produced by current match
        //
        // RETURN:
        //          true iff we have a match
        inline bool isMatch(uint16_t& errNum);

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
        std::array<uint8_t, 256>& lmap;


};


template<size_t E>
ShiftAnd<E>::ShiftAnd(std::string& seq, std::array<uint8_t, 256>& lMap) :
        pLen(seq.size())
    ,   lmap(lMap)
{
    loadBitmasks(seq);
}


template<size_t E>
inline void ShiftAnd<E>::querySeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<uint64_t>& matches, std::vector<uint16_t>& errors)
{

    reset();

    for (auto it = start; it != end; ++it)
    {

        queryLetter(*it);

        uint16_t errNum;
        if (isMatch(errNum))
        {

            matches.emplace_back(it - start);
            errors.emplace_back(errNum);

        }
    }
}


template<size_t E>
inline void ShiftAnd<E>::queryRevSeq(std::vector<char>::iterator start, std::vector<char>::iterator end, std::vector<uint64_t>& matches, std::vector<uint16_t>& errors)
{

    reset();

    for (auto it = start; it != end; --it)
    {

        char c;
        // change letter to reverse complement
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
            default:
                std::cerr << "[ShiftAnd] Could not parse letter " << *it << "\n\n";
                exit(1);
        }
        queryLetter(c);

        uint16_t errNum;
        if (isMatch(errNum))
        {

            matches.emplace_back(it - end);
            errors.emplace_back(errNum);

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


// TODO: make template spezialisations for this for E=1,2
//          OR use loop unrolling with templates (see bottom of file)
template<size_t E>
inline void ShiftAnd<E>::queryLetter(const char& c)
{

    const bitMasks& mask = masks[lmap[c]];

    // Bottom up update part for old values of previous iteration
    for (size_t i = E; i > 0; --i)
    {

        // Update first part of pattern states
        //
        //                          Match                           Insertion               Substitution
        active[i].B_0 = ((active[i].B_0 << 1 | 1) & mask.B_0) | (active[i-1].B_0) | (active[i-1].B_0 << 1);

        // Update second part of pattern states
        //
        //                                  Match                                       Insertion                       Substitution
        active[i].B_1 = ((active[i].B_1 << 1 | active[i].B_0 >> 63) & mask.B_1) | (active[i-1].B_1) | (active[i-1].B_1 << 1 | active[i-1].B_0 >> 63);


    }
    // update zero error layer (at the top)
    active[0].B_0 = ((active[0].B_0 << 1 | 1) & mask.B_0);
    active[0].B_1 = ((active[0].B_1 << 1 | active[0].B_0 >> 63) & mask.B_1);
    //
    // Top down update for values of this iteration
    for (size_t i = 1; i <= E; ++i)
    {

        active[i].B_0 |= active[i-1].B_0 << 1;
        active[i].B_1 |= active[i-1].B_1 << 1 | active[i-1].B_0 >> 63;

    }
}


template<size_t E>
inline bool ShiftAnd<E>::isMatch(uint16_t& errNum)
{
    // go through layers
    for (size_t i = 0; i <= E; ++i)
    {

        // test if this layer is match
        if ((active[i].B_0 & accepted.B_0) | (active[i].B_1 & accepted.B_1))
        {
            errNum = i;
            return true;
        }

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
        masks[lmap['A']].B_0 = maskA;
        masks[lmap['C']].B_0 = maskC;
        masks[lmap['G']].B_0 = maskG;
        masks[lmap['T']].B_0 = maskT;



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
        masks[lmap['A']].B_1 = maskA;
        masks[lmap['C']].B_1 = maskC;
        masks[lmap['G']].B_1 = maskG;
        masks[lmap['T']].B_1 = maskT;

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
        masks[lmap['A']].B_0 = maskA;
        masks[lmap['C']].B_0 = maskC;
        masks[lmap['G']].B_0 = maskG;
        masks[lmap['T']].B_0 = maskT;

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
        masks[lmap['A']].B_1 = maskA;
        masks[lmap['C']].B_1 = maskC;
        masks[lmap['G']].B_1 = maskG;
        masks[lmap['T']].B_1 = maskT;

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
        masks[lmap['A']].B_0 = maskA;
        masks[lmap['C']].B_0 = maskC;
        masks[lmap['G']].B_0 = maskG;
        masks[lmap['T']].B_0 = maskT;
    }

}


// COMPILE TIME loop structure to unroll loops for bit updates


// Calling example:
//
// ForLoop<MyConst::MISCOUNT>::iterate<BottomUpChanges>()
// template<size_t loopCounter>
// struct ForLoop
// {
//     template<template <size_t> class function>
//     static void iterate(ShiftAnd<MyConst::MISCOUNT> sa)
//     {
//         function<loopCounter>::apply(sa);
//         ForLoop<loopCounter - 1>::template iterate<function>(sa);
//     }
// };
//
// template<>
// struct ForLoop<1>
// {
//     template<template <size_t> class function>
//     static void iterate(ShiftAnd<MyConst::MISCOUNT> sa)
//     {
//         function<1>::apply(sa);
//     }
// };
//
// // function should only be applied to error leayers - this is edge case
// // for zero mismatch program executions
// template<>
// struct ForLoop<0>
// {
//
//     template<template <size_t> class function>
//     static void iterate(ShiftAnd<MyConst::MISCOUNT> sa)
//     {
//     }
// };
//
// template <size_t layer>
// struct BotomUpChanges
// {
//     // TODO: pass mask
//     static void apply(ShiftAnd<MyConst::MISCOUNT> sa, struct states mask)
//     {
//
//         // Update first part of pattern states
//         //
//         //                                      Match                           Insertion           Substitution
//         sa.active[layer].B_0 = ((sa.active[layer].B_0 << 1 | 1) & mask.B_0) | (sa.active[layer-1].B_0) | (sa.active[layer-1].B_0 << 1);
//
//         // Update second part of pattern states
//         //
//         //                                      Match                                             Insertion                       Substitution
//         sa.active[layer].B_1 = ((sa.active[layer].B_1 << 1 | sa.active[layer].B_0 >> 63) & mask.B_1) | (sa.active[layer-1].B_1) | (sa.active[layer-1].B_1 << 1 | sa.active[layer-1].B_0 >> 63);
//     }
// };
//
// template <size_t layer>
// struct TopDownChanges
// {
//
//     static void apply(ShiftAnd<MyConst::MISCOUNT> sa)
//     {
//
//     }
// };


#endif /* SHIFTAND_H */
