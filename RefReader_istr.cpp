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

#include <fstream>
#include <algorithm>

#include "RefReader_istr.h"


void readReference(const std::string& filename, std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<std::vector<char> >& genSeq, std::unordered_map<uint8_t, std::string>& chrMap)
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

    std::cout << "Start reading reference file " << filename << std::endl;

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
                seq = std::vector<char>();
                seq.reserve(MyConst::CHROMMAX);
                // reset flag for last c
                lastC = false;

            }
            // check if primary assembly
			// (GRCH versions)
            if (*(line.begin() + 1) == 'C')
            {

                ++chrIndex;
				std::string chrID (line.begin() + 1, line.end());
				chrMap.insert(std::pair<uint8_t,std::string>(chrIndex, chrID));
                contFlag = true;
                continue;

            // throw out unlocalized contigs
			// (hg versions)
			} else if ( *(line.begin() + 1) == 'c')
			{

                ++chrIndex;
				std::string chrID (line.begin() + 1, line.end());
				// test if real primary assembly sequence
				if (!isPrimaryHG(chrID))
				{
					--chrIndex;
					contFlag = false;
					continue;
				}
				chrMap.insert(std::pair<uint8_t,std::string>(chrIndex, chrID));
                contFlag = true;
                continue;

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

    // if we read primary assembly previously, write it to vectors
    if (contFlag)
    {

        seq.shrink_to_fit();
        genSeq.emplace_back(move(seq));
    }

    cpgTab.shrink_to_fit();
    genSeq.shrink_to_fit();
    ifs.close();

    std::cout << "Done reading reference file" << std::endl;

}

bool isPrimaryHG(std::string chrID)
{
	for (unsigned int i = 1; i <= 22; ++i)
	{
		if (chrID == ("chr" + std::to_string(i)))
		{
			return true;
		}
	}
	if (chrID == "chrX")
		return true;
	if (chrID == "chrY")
		return true;

	return false;

}
