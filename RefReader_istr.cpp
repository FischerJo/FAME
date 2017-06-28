
#include <fstream>
#include <algorithm>

#include "RefReader_istr.h"


void readReference(const char* const filename, std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<std::vector<char> >& genSeq)
{

    std::string line;

    // stores the sequence of each chromosome
    std::vector<char> seq;

    genSeq.reserve(MyConst::CHROMNUM);
    seq.reserve(MyConst::CHROMMAX);
    cpgTab.reserve(MyConst::CPGMAX);

    // stores the chromosome index we are currently reading
    uint8_t chrIndex = 0;

    // flag stating if we currently read a real chromosome assembly
    bool contFlag = false;
    // flag stating if the last character read in the previous line was a c
    bool lastC = false;

    std::ifstream ifs (filename);

    std::cout << "Start reading file" << std::endl;

    // TODO
    // unsigned int counter = 0;

    while (getline(ifs, line)) {

        // if (counter >= 3121254)
        // {
        //     break;
        // }
        // ++counter;
        // Test for id tag line
        if (line.begin() != line.end() && *(line.begin()) == '>' && (line.begin() + 1) != line.end())
        {

            // if we read primary assembly previously, write it to vectors
            if (contFlag)
            {

                seq.shrink_to_fit();
                genSeq.emplace_back(move(seq));
                // reset buffer
                seq = std::vector<char>();
                seq.reserve(MyConst::CHROMMAX);
                // reset flag for last c
                lastC = false;

            }
            // check if primary assembly
            if (*(line.begin() + 1) == 'C')
            {

                ++chrIndex;
                contFlag = true;
                continue;

            // throw out unlocalized contigs
            } else {

                contFlag = false;
                continue;
            }

        }

        // check if we are in real primary chromosome assembly
        if (contFlag)
        {

            // parse the line
            readLine(line, lastC, chrIndex, cpgTab, cpgStartTab, seq);

        } else {

            continue;
        }



    }

    std::cout << "Done parsing file" << std::endl;

    // if we read primary assembly previously, write it to vectors
    if (contFlag)
    {

        seq.shrink_to_fit();
        genSeq.emplace_back(move(seq));
    }

    cpgTab.shrink_to_fit();
    genSeq.shrink_to_fit();
    ifs.close();
}
