
#include "Read.h"


Read::Read(std::string& sequence, std::string& identifier) :
        id(identifier)
    ,   seq(sequence)
    ,   matches()
    ,   isInvalid(false)
{
}
