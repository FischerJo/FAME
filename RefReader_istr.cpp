
#include <fstream>

#include "RefReader_istr.h"

using namespace std;


void readReference(const char* const filename, vector<struct CpG>& cpgTab, vector<vector<char> >& genSeq)
{

    string line;

    // stores the sequence of each chromosome
    string seq;
    seq.reserve(MyConst::CHROMMAX);

    // stores the chromosome index we are currently reading
    uint8_t chrIndex = -1;

    // flag stating if we currently read a real chromosome assembly
    bool contFlag = false;
    // flag stating if the last character read in the previous line was a c
    bool lastC = false;

    ifstream ifs (filename);

    while (getline(ifs, line)) {

        // Test for id tag line
        if (line.begin() != line.end() && *(line.begin()) == '>' && (line.begin() + 1) != line.end())
        {

            // if we read primary assembly previously, write it to vectors
            if (contFlag)
            {

                // put a vector of chars of size equal to the stringlength read so far to genSeq
                genSeq.emplace_back(seq.size());
                // copy the content of sequence to the object holding all sequences
                copy(seq.begin(), seq.end(), genSeq[chrIndex].begin());
                // reset buffer
                seq.clear();
                // reset flag for last c
                lastC = false;

            }
            // check if primary assembly
            if (*(++(line.begin())) == 'C')
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

            readLine(line, lastC, chrIndex, cpgTab, seq);

        } else {

            continue;
        }



    }
}
