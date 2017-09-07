//	Metal - A fast methylation alignment and calling tool for WGBS data.
//	Copyright (C) 2017  Jonas Fischer
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	Jonas Fischer	jonaspost@web.de

#ifndef REFREADER_ISTR_H
#define REFREADER_ISTR_H

#include <vector>
#include <string>

#include <iostream> // for debugging

#include "structs.h"
#include "CONST.h"

// read reference FASTA file from filename
// produces a vector of all CpGs present in the genome written to cpgTab and cpgStartTab (the latter one is filled with those CpGs that are not READLEN away from the start)
// IMPORTANT:
//              cpg.pos for all CpGs in cpgStartTab will be the offset to the C in CpG!
// produces sequence strings seperated by chromosome saved to genSeq, their length to genSeqLen
//      underlying vectors should be empty on calling
//      convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
void readReference(const std::string& filename, std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<std::vector<char> >& genSeq);

// read a line of the reference, constructing possible CpGs, appending and counting the read characters
//
// Arguments:
//              line        line to parse
//              lastC       true if previous line encoding this sequence ended with a C
//                          will be updated after the function call is finished
//              chrIndex    index of the current chromosome (start counting from 1!)
//              cpgTab      reference to where the CpGs should be saved (will be padded at the end)
//              seq         sequence to which chars should be appended
inline void readLine(std::string& line, bool& lastC, uint8_t chrIndex, std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<char>& seq)
{

    chrIndex = chrIndex - 1;
    // check for G if last letter in previous line was C
    // if (lastC)
    // {
    //
    //     std::string::iterator it = line.begin();
    //     if (it != line.end())
    //     {
    //
    //         if (*it == 'G')
    //         {
    //             // start reading for CpG READLEN - 2 chars before, offsets are started counting at 0 but string.size() starts at 1
    //             // so we need to start at sequencesize - (READLEN - 2) - 1
    //             if (seq.size() < (MyConst::READLEN - 1))
    //             {
    //
    //                 // note that seq.size() > 1 because lastC is only set if we read a C before
    //                 cpgStartTab.emplace_back(chrIndex, seq.size() - 1);
    //
    //             } else {
    //
    //                 cpgTab.emplace_back(chrIndex, seq.size() - MyConst::READLEN + 1);
    //             }
    //
    //         }
    //     }
    // }
    //
    // // search for CpGs
    // for (std::string::iterator it = line.begin(); it != line.end(); ++it)
    // {
    //
    //     seq += *it;
    //
    //     // test if letter is C
    //     if (*it == 'C')
    //     {
    //
    //         // if there exists next letter, check if g
    //         if ( (it + 1) != line.end() )
    //         {
    //
    //             if (*(it + 1) == 'G')
    //             {
    //                 if (seq.size() < (MyConst::READLEN - 1))
    //                 {
    //
    //                     cpgStartTab.emplace_back(chrIndex, seq.size() - 1);
    //
    //                 } else {
    //
    //                     cpgTab.emplace_back(chrIndex, seq.size() - MyConst::READLEN + 1);
    //                 }
    //                 seq += *(++it);
    //             }
    //
    //         // else update lastC because we are at the end of line
    //         } else {
    //
    //             lastC = true;
    //             return;
    //         }
    //     }
    // }
    // lastC = false;
    //

    // parse the line
    // transform lower case to upper case, construct cpg objects and write unknowns as N
    for (char& c : line)
    {

        switch (c)
        {

            case 'g':
            case 'G':

                // if previous letter was C, construct a CpG object
                if (lastC)
                {
                    if (seq.size() >= (MyConst::READLEN - 1))
                    {

                        cpgTab.push_back({chrIndex, static_cast<uint32_t>(seq.size()) - MyConst::READLEN + 1});

                    } else {

                        cpgStartTab.push_back({chrIndex, static_cast<uint32_t>(seq.size()) - 1});
                    }
                    lastC = false;
                }
                seq.emplace_back('G');
                break;

            case 'c':
            case 'C':

                seq.emplace_back('C');
                lastC = true;
                break;

            case 'a':
            case 'A':

                seq.emplace_back('A');
                lastC = false;
                break;

            case 't':
            case 'T':

                seq.emplace_back('T');
                lastC = false;
                break;

            default:

                seq.emplace_back('N');
                lastC = false;
                break;


        }
    }

}


#endif /* REFREADER_ISTR_H */
