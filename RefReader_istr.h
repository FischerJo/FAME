#ifndef REFREADER_ISTR_H
#define REFREADER_ISTR_H

#include <vector>
#include <string>

#include "structs.h"
#include "CONST.h"

// read reference FASTA file from filename
// produces a vector of all CpGs present in the genome written to cpgTab
// produces sequence strings seperated by chromosome saved to genSeq, their length to genSeqLen
//      underlying vectors should be empty on calling
//      convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
void readReference(const char* const filename, std::vector<struct CpG>& cpgTab, std::vector<std::vector<char> >& genSeq);

// read a line of the reference, constructing possible CpGs, appending and counting the read characters
//
// Arguments:
//              line        line to parse
//              lastC       true if previous line encoding this sequence ended with a C
//                          will be updated after the function call is finished
//              chrIndex    index of the current chromosome
//              cpgTab      reference to where the CpGs should be saved (will be padded at the end)
//              seq         sequence to which chars should be appended
inline void readLine(std::string& line, bool& lastC, const uint8_t chrIndex, std::vector<struct CpG>& cpgTab, std::string& seq)
{

    // check for G if last letter in previous line was C
    if (lastC)
    {

        std::string::iterator it = line.begin();
        if (it != line.end())
        {

            if (*it == 'G' || *it == 'g')
            {
                // start reading for CpG READLEN - 2 chars before, offsets are started counting at 0 but string.size() starts at 1
                // so we need to start at sequencesize - (READLEN - 2) - 1
                cpgTab.emplace_back(chrIndex, seq.size() - MyConst::READLEN + 1);

            }
        }
    }

    // search for CpGs
    for (std::string::iterator it = line.begin(); it != line.end(); ++it)
    {

        seq += *it;

        // test if letter is C
        if (*it == 'C' || *it == 'c')
        {

            // if there exists next letter, check if g
            if ( (++it) != line.end() )
            {

                seq += *it;
                if (*it == 'G' || *it == 'g')
                {
                    cpgTab.emplace_back(chrIndex, seq.size() - MyConst::READLEN + 1);

                }

            // else update lastC because we are at the end of line
            } else {

                lastC = true;
                return;
            }
        }
    }
    lastC = false;

}


#endif /* REFREADER_ISTR_H */
