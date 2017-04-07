
#include <fstream>
#include <algorithm>

#include "RefReader_istr.h"

using namespace std;


void readReference(const char* const filename, vector<struct CpG>& cpgTab, vector<struct CpG>& cpgStartTab, vector<vector<char> >& genSeq)
{

    string line;

    // stores the sequence of each chromosome
    vector<char> seq;

    genSeq.reserve(MyConst::CHROMNUM);
    seq.reserve(MyConst::CHROMMAX);
    cpgTab.reserve(MyConst::CPGMAX);

    // stores the chromosome index we are currently reading
    uint8_t chrIndex = 0;

    // flag stating if we currently read a real chromosome assembly
    bool contFlag = false;
    // flag stating if the last character read in the previous line was a c
    bool lastC = false;

    ifstream ifs (filename);

    cout << "Start reading file" << endl;

    while (getline(ifs, line)) {

        // Test for id tag line
        if (line.begin() != line.end() && *(line.begin()) == '>' && (line.begin() + 1) != line.end())
        {

            // if we read primary assembly previously, write it to vectors
            if (contFlag)
            {

                seq.shrink_to_fit();
                genSeq.emplace_back(move(seq));
                // reset buffer
                seq = vector<char>();
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

    cout << "Done parsing file" << endl;

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
