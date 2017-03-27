
#include "DnaBitStr.h"


using namespace std;


DnaBitStr::DnaBitStr(const unsigned int s) :
        size(s)
    ,   bitSeq((size / 64) + 1)
    ,   bitMask((size / 64) + 1)
    ,   bitRevMask((size / 64) + 1)
{
}



void DnaBitStr::setBitStrLast(string& seq)
{

    uint64_t bitStr = 0;
    uint64_t bitM = 0xffffffffffffffffULL;
    uint64_t bitRevM = 0xffffffffffffffffULL;
    for (unsigned int i = 1; i <= seq.size(); ++i)
    {

        // we don't need A here
        switch (seq[i - 1]){

            case 'C':
                {
                    const unsigned int shift = (64 - 2*i);
                    bitStr |= (1ULL << shift);
                    bitM ^= (3ULL << shift);
                }
                break;

            case 'G':
                {
                    const unsigned int shift = (64 - 2*i);
                    bitStr |= (2ULL << shift);
                    bitRevM ^= (3ULL << shift);
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
    bitSeq[size-1] = bitStr;
    bitMask[size-1] = bitM;
}
