
#include <iostream>

#include "Extract.h"

int main(int args, char** argv)
{

    Extractor ex (argv[1]);
    ex.parseLines();
    std::cout << ex.getCount() << "\n\n";
    return 0;
}
