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


void readReference(const std::string& filename, std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<std::vector<char> >& genSeq, std::unordered_map<uint8_t, std::string>& chrMap, const bool humanOptFlag)
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
	if(!ifs)
    {
        std::cerr << "Opening genome reference file " << filename << " was unsuccessful! Exiting..." << std::endl;
		exit(1);
    }

    std::cout << "Start reading reference file " << filename << std::endl;

	uint32_t unIDCount = 1;

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
			if (humanOptFlag)
			{
				// check if primary assembly
				// (GRCH versions)
				if (*(line.begin() + 1) == 'C' || (*(line.begin() + 1) == 'N' && *(line.begin() + 2) == 'C'))
				{

					++chrIndex;
					std::string chrHeader (line.begin() + 1, line.end());
					uint32_t pos = chrHeader.find_first_of(" \t");
					std::string chrID(chrHeader, 0, pos);
					for (const auto chrP : chrMap)
					{
						if (chrP.second == chrID)
						{
							std::cout << "WARNING: Chromosome identifier " << chrID << " found in header\n" <<
								chrHeader << "\nis not unique.";
							chrID.append("_");
							chrID.append(std::to_string(unIDCount++));
							std::cout << "Renaming to " << chrID << "\n";
							break;
						}
					}
					chrMap.insert(std::pair<uint8_t,std::string>(chrIndex - 1, chrID));
					contFlag = true;
					continue;

				// throw out unlocalized contigs
				// (hg versions)
				} else if ( *(line.begin() + 1) == 'c')
				{

					++chrIndex;
					std::string chrHeader (line.begin() + 1, line.end());
					uint32_t pos = chrHeader.find_first_of(" \t");
					std::string chrID(chrHeader, 0, pos);
					for (const auto chrP : chrMap)
					{
						if (chrP.second == chrID)
						{
							std::cout << "WARNING: Chromosome identifier " << chrID << " found in header\n" <<
								chrHeader << "\nis not unique.";
							chrID.append("_");
							chrID.append(std::to_string(unIDCount++));
							std::cout << "Renaming to " << chrID << "\n";
							break;
						}
					}
					// test if real primary assembly sequence
					if (!isPrimaryHG(chrID))
					{
						--chrIndex;
						contFlag = false;
						continue;
					}
					chrMap.insert(std::pair<uint8_t,std::string>(chrIndex - 1, chrID));
					contFlag = true;
					continue;

				} else {
					contFlag = false;
				}
            } else {

				++chrIndex;
				if (chrIndex > MyConst::CHROMNUM)
				{
					std::cout << "\nWARNING: Number of read chromosomes and number of specified chromosomes the organism should have do not match!\n"
						<< "Maybe the reference genome contains unlocalized contigs. You should consider removing them from the reference fasta.\n"
						<< "For the human reference genome GRCH and HG versions, use \"--human_opt\" to mitigate this problem.\n";
				}
				std::string chrHeader (line.begin() + 1, line.end());
				uint32_t pos = chrHeader.find_first_of(" \t");
				std::string chrID(chrHeader, 0, pos);
				for (const auto chrP : chrMap)
				{
					if (chrP.second == chrID)
					{
						std::cout << "WARNING: Chromosome identifier " << chrID << " found in header\n" <<
							chrHeader << "\nis not unique.";
						chrID.append("_");
						chrID.append(std::to_string(unIDCount++));
						std::cout << "Renaming to " << chrID << "\n";
						break;
					}
				}
				std::cout << chrHeader << "\n" << chrID << "\n\n";
				chrMap.insert(std::pair<uint8_t,std::string>(chrIndex - 1, chrID));
				contFlag = true;
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
