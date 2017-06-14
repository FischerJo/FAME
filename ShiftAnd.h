#ifndef SHIFTAND_H
#define SHIFTAND_H

#include <vector>
#include <array>
#include <cstdint>

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
        ShiftAnd(std::vector<char>& seq, std::array<uint8_t, 256>& lMap);

        // -------------------

        // Query the sequence slice comprised by [start,end) to the internal automata
        //
        // ARGUMENTS:
        //              start   iterator to start of sequence to query
        //              end     iterator to end of sequence to query
        //
        // RETURN:
        //              vector of offsets (relative to start iterator) were matchings start
        std::vector<unsigned int> querySeq(std::vector<char>::iterator& start, std::vector<char>::iterator& end);


    private:

        // Reset to initial state
        // is called by querySeq
        void reset();

        // load the bitmasks for the given sequences
        void loadBitmasks(std::vector<char>& seq);

        // Bitvector data structure holding the set of active states indicated by a 1
        // initial state is least significant bit of active[0].B_0 (dummy state),
        // different layers indicating the number of errors are indexed by the vector
        std::array<bitStates, E + 1> active;

        // Set of bitmasks for the given sequence
        // Generally, Masks[m] give the bitmask for letter m
        std::array<bitMasks, 4> masks;
        // masks to filter out accepting states
        bitMasks accepted;


        // maps characters 'A', 'C', 'G', 'T' to their index
        // // TODO: init this in ReadQueue
        std::array<uint8_t, 256>& lmap;


};

#endif /* SHIFTAND_H */
