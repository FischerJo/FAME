
#include "ShiftAnd.h"

template<size_t E>
ShiftAnd<E>::ShiftAnd(std::vector<char>& seq, std::array<uint8_t, 256>& lMap) :
        lmap(lMap)
{
    loadBitmasks(seq);
}

template<size_t E>
std::vector<unsigned int> ShiftAnd<E>::querySeq(std::vector<char>::iterator& start, std::vector<char>::iterator& end)
{
    std::vector<unsigned int> matches;

    return matches;
}


template<size_t E>
void ShiftAnd<E>::reset()
{
}

template<size_t E>
void ShiftAnd<E>::loadBitmasks(std::vector<char>& seq)
{

    // retrieve length of the sequence
    const size_t seqLen = seq.size();
    auto seqIt = seq.crbegin();
    // if sequence is too large, only represent first part. move reverse iterator appropriately.
    if (seqLen > (128 - 1))
    {

        seqIt += seqLen - (128 - E - 1);

    }
    int

    // go through sequence from front to back
    for ( ; seqIt != seq.crend(); ++seqIt)
    {

        switch (seq[i])
        {

            case 'A':


        }

    }
}
