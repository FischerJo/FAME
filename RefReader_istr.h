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
void readReference(const char* const filename, std::vector<struct CpG>& cpgTab, std::vector<const char*>& genSeq, std::vector<std::size_t> genSeqLen);

// read a line of the reference, constructing possible CpGs, appending and counting the read characters
//
// Arguments:
//              line        line to parse
//              lastC       true if previous line encoding this sequence ended with a C
//                          will be updated after the function call is finished
//              chrIndex    index of the current chromosome
//              cpgTab      reference to where the CpGs should be saved (will be padded at the end)
//              seq         sequence to which chars should be appended
//              seqLen      current length of this sequence
inline void readLine(std::string& line, bool& lastC, const uint8_t chrIndex, std::vector<struct CpG>& cpgTab, std::string& seq, unsigned int& seqLen)
{
    // append line to sequence
    seq.append(line);
    // search for CpGs
    for (string::iterator it = line.begin(); it != line.end(); ++it)
    {
        ++seqLen;

        // test if letter is C
        if (*it == 'C' || *it == 'c')
        {

            // if there exists next letter, check if g
            if ( (++it) != line.end() )
            {
                ++seqLen;

                if (*it == 'G' || *it == 'g')
                {
                    cpgTab.emplace_back(chrIndex, seqLen - MyConst::READLEN + 2);

                }

            // else update lastC
            } else {

                lastC = true;
                return;
            }

    }
    lastC = false;

}


#endif /* REFREADER_ISTR_H */
