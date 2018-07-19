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
#include <unordered_map>

#include <iostream> // for debugging

#include "structs.h"
#include "CONST.h"

// read reference FASTA file from filename
// produces a vector of all CpGs present in the genome written to cpgTab and cpgStartTab (the latter one is filled with those CpGs that are not READLEN away from the start)
// IMPORTANT:
//              cpg.pos for all CpGs in cpgStartTab will be the offset to the C in CpG!
// produces sequence strings seperated by chromosome saved to genSeq, their length to genSeqLen
//      underlying vectors should be empty on calling
void readReference(const std::string& filename, std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<std::vector<char> >& genSeq, std::unordered_map<uint8_t,std::string>& chrMap, const bool humanOptFlag);

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

// Arguments:
// 				chrID	string containing chromosome ID (i.e. content after '>' id line in fasta)
// Return:
// 				true if real primary assembly (i.e. not unlocalized contigs etc)
bool isPrimaryHG(std::string chrID);

#endif /* REFREADER_ISTR_H */
