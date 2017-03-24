
#include "DnaBitStr.h"


using namespace std;


DnaBitStr::DnaBitStr(const unsigned int size) :
        this->size(size)
    ,   bitSeq((size / 64) + 1)
    ,   bitMask((size / 64) + 1)
{
}



void DnaBitStr::setBitStrN(string& seq, const unsigned int n)
{


}
