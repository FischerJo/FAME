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
	// one layer of states is split into 64bit int that represents 64 consecutive states
    using bitMasks = std::vector<uint64_t>;
    using bitStates = std::vector<uint64_t>;

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
                if (errNum <= prevErrs)
                {
                    matches.back() = it - start + MyConst::MISCOUNT;
                    errors.back() = errNum;
                    prevErrs = errNum;
                }

            } else {

                matches.push_back(it - start + MyConst::MISCOUNT);
                errors.push_back(errNum);
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

        switch (*it)
        {
            case 'A':
				queryLetter('T');
                break;
            case 'C':
				queryLetter('G');
                break;
            case 'G':
				queryLetter('C');
                break;
            case 'T':
				queryLetter('A');
                break;
            // we do not consider Ns for matches - restart whole automaton for next letter
            case 'N':
                reset();
                continue;
            default:
                std::cerr << "[ShiftAnd] Could not parse letter " << *it << " at offset " << it - end + 1 << "\nStopping program...\n\n";
                exit(1);
        }

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
                if (errNum <= prevErrs)
                {
					// TODO make this safe for insertions!!
                    matches.back() = it - end + pLen - 1 + MyConst::MISCOUNT;
                    errors.back() = errNum;
                    prevErrs = errNum;
                }

            } else {

				// TODO make this safe for insertions!!
                matches.push_back(it - end + pLen - 1 + MyConst::MISCOUNT);
                errors.push_back(errNum);
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
		for (uint64_t j = active[i].size() - 1; j > 0; --j)
		{
			active[i][j] = 0;
		}
        active[i][0] = (static_cast<uint64_t>(1) << (i+1)) - 1;

    }
}


template<size_t E>
inline void ShiftAnd<E>::queryLetter(const char& c)
{

	const bitMasks& mask = masks[lmap[c%16]];
	// Bottom up update part for old values of previous iteration
	for (size_t i = E; i > 0; --i)
	{

		// from right to left through states
		for (uint64_t j = active[i].size() - 1; j > 0; --j)
		{
			// Update j-th part of pattern states
			//
			//                                  Match                                       Insertion                       Substitution
			active[i][j] = ((active[i][j] << 1 | active[i][j-1] >> 63) & mask[j]) | (active[i-1][j]) | (active[i-1][j] << 1 | active[i-1][j-1] >> 63);
		}

		// Update first part of pattern states
		//
		//                          Match                           Insertion               Substitution
		active[i][0] = ((active[i][0] << 1 | 1) & mask[0]) | (active[i-1][0]) | (active[i-1][0] << 1);

	}
	// update zero error layer (at the top)
	for (uint64_t j = active[0].size() - 1; j > 0; --j)
	{
		active[0][j] = ((active[0][j] << 1 | active[0][j-1] >> 63) & mask[j]);
	}
	active[0][0] = ((active[0][0] << 1 | 1) & mask[0]);
	//
	// Top down update for values of this iteration
	for (size_t i = 1; i <= E; ++i)
	{

		// from right to left through states
		for (uint64_t j = active[i].size() - 1; j > 0; --j)
		{
			active[i][j] |= active[i-1][j] << 1 | active[i-1][j-1] >> 63;
		}
		active[i][0] |= active[i-1][0] << 1;

	}
}



template<size_t E>
inline bool ShiftAnd<E>::isMatch(uint8_t& errNum)
{
    // go through layers
    for (size_t i = 0; i <= E; ++i)
    {

		for (size_t j = 0; j < active[i].size(); ++j)
		{
			// test if this layer is match
			if (active[i][j] & accepted[j])
			{
				errNum = i;
				return true;
			}
		}

    }
    return false;
}


template<size_t E>
inline void ShiftAnd<E>::loadBitmasks(std::string& seq)
{

    // retrieve length of the sequence
    const size_t seqLen = seq.size();
	// initializing
	accepted.resize(seqLen/64 + 1);
	masks[lmap['A'%16]].resize(seqLen/64 + 1);
	masks[lmap['C'%16]].resize(seqLen/64 + 1);
	masks[lmap['T'%16]].resize(seqLen/64 + 1);
	masks[lmap['G'%16]].resize(seqLen/64 + 1);
	for (auto& ac : active)
	{
		ac.resize(seqLen/64 + 1);
	}

	for (size_t i = 0; i < seqLen / 64; ++i)
	{
		accepted[i] = 0;
	}
	accepted[seqLen/64] = static_cast<uint64_t>(1) << (seqLen%64);

	// load full mask for trailing 1s for sequences that are shorter than pLen
	uint64_t maskA = 0xffffffffffffffffULL;
	uint64_t maskC = 0xffffffffffffffffULL;
	uint64_t maskG = 0xffffffffffffffffULL;
	uint64_t maskT = 0xffffffffffffffffULL;

	for (uint64_t i = seqLen; i > 0; )
	{

		// shift to make space for next letter
		maskA = maskA << 1;
		maskC = maskC << 1;
		maskG = maskG << 1;
		maskT = maskT << 1;

		switch (seq[i-1])
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
		if (i%64 == 0)
		{
			// save masks
			masks[lmap['A'%16]][i/64] = maskA;
			masks[lmap['C'%16]][i/64] = maskC;
			masks[lmap['G'%16]][i/64] = maskG;
			masks[lmap['T'%16]][i/64] = maskT;
			// load full mask for trailing 1s for sequences that are shorter than pLen
			maskA = 0xffffffffffffffffULL;
			maskC = 0xffffffffffffffffULL;
			maskG = 0xffffffffffffffffULL;
			maskT = 0xffffffffffffffffULL;
		}
		--i;

	}
	// load initial state
	maskA = (maskA << 1) | 1;
	maskC = (maskC << 1) | 1;
	maskG = (maskG << 1) | 1;
	maskT = (maskT << 1) | 1;
	// save masks
	masks[lmap['A'%16]][0] = maskA;
	masks[lmap['C'%16]][0] = maskC;
	masks[lmap['G'%16]][0] = maskG;
	masks[lmap['T'%16]][0] = maskT;

}


#endif /* SHIFTAND_H */
