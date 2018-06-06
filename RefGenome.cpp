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

#include <chrono>
#include <algorithm> // max

#include "RefGenome.h"


RefGenome::RefGenome(std::vector<struct CpG>&& cpgTab, std::vector<struct CpG>&& cpgStartTab, std::vector<std::vector<char> >& genomeSeq, const bool noloss, std::unordered_map<uint8_t, std::string>& chromMap) :
        cpgTable(cpgTab)
    ,   cpgStartTable(cpgStartTab)
    ,   fullSeq(std::move(genomeSeq))
    ,   tabIndex(MyConst::HTABSIZE + 1, 0)
    ,   kmerTable()
    ,   strandTable()
    ,   metaCpGs()
    ,   metaStartCpGs()
	,	chrMap(chromMap)
{

    if (!testPODs())
    {
        std::cout << "\nWARNING: Structures are not POD! Cannot write index to file savely!\n\n";
    }
    // find out genome size
    size_t gensize = 0;
    for (std::vector<char> chr : genomeSeq)
    {
        gensize += chr.size();
    }
    // init meta table with upper bound on required windows
    metaCpGs.reserve(gensize/MyConst::WINLEN);
    metaStartCpGs.reserve(MyConst::CHROMNUM);
    // fill meta table
    generateMetaCpGs();
    // generate encoding of genome
    // generateBitStrings(fullSeq);
    // cout << "Done generating Genome bit representation" << endl;
    std::cout << "\nStart hashing CpGs\n";
    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    // hash all kmers of reduced alphabet
    generateHashes(fullSeq);
    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    std::cout << "\nDone hashing CpGs (" << runtime << "s)\n";
    if (!noloss)
    {
        std::cout << "\nStarting filtering process for Index\nThrowing out highly repetetive kmers...\n";
        // filter out highly repetitive sequences
        filterHashTable();
    }
    std::cout << "\nThrowing out kmers of single metaCpG with same hash...\n";
    filterRedundancyInHashTable();
    std::cout << "\nFinished index processing.\n";
}


RefGenome::RefGenome(std::string filepath)
{
    load(filepath);
}


void RefGenome::save(const std::string& filepath)
{

    std::cout << "Start writing index to file " << filepath << "\n";
    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    // open file in binary mode
    std::ofstream of(filepath, std::ofstream::binary | std::ofstream::trunc);

    // save CONSTANTS
    auto htabs = MyConst::HTABSIZE;
    of.write(reinterpret_cast<char*>(&htabs), sizeof(htabs));
    auto readl = MyConst::READLEN;
    of.write(reinterpret_cast<char*>(&readl), sizeof(readl));
    auto winl = MyConst::WINLEN;
    of.write(reinterpret_cast<char*>(&winl), sizeof(winl));
    auto kmerl = MyConst::KMERLEN;
    of.write(reinterpret_cast<char*>(&kmerl), sizeof(kmerl));
    auto kmerc = MyConst::KMERCUTOFF;
    of.write(reinterpret_cast<char*>(&kmerc), sizeof(kmerc));

    // store NON-start CpGs
    size_t cpgNum = cpgTable.size();
    of.write(reinterpret_cast<char*>(&cpgNum), sizeof(cpgNum));
    of.write(reinterpret_cast<char*>(cpgTable.data()), sizeof(cpgTable[0]) * cpgNum);
    if (of.fail())
    {
        std::cout << "Error writing CpG stuff " << cpgNum << "\n";
    }
    // store start CpGs
    cpgNum = cpgStartTable.size();
    of.write(reinterpret_cast<char*>(&cpgNum), sizeof(cpgNum));
    of.write(reinterpret_cast<char*>(cpgStartTable.data()), sizeof(cpgStartTable[0]) * cpgNum);
    if (of.fail())
    {
        std::cout << "Error writing CpG Start stuff " << cpgNum << "\n";
    }
    // write reference sequence
    size_t chromNum = fullSeq.size();
    of.write(reinterpret_cast<char*>(&chromNum), sizeof(chromNum));
    for (std::vector<char> chromSeq : fullSeq)
    {
        size_t chromLen = chromSeq.size();
        of.write(reinterpret_cast<char*>(&chromLen), sizeof(chromLen));
        of.write(chromSeq.data(), chromLen);
        if (of.fail())
        {
            std::cout << "Error writing sequence " << chromLen << "\n";
        }
    }
    // tabIndex
    size_t tabLen = tabIndex.size();
    of.write(reinterpret_cast<char*>(&tabLen), sizeof(tabLen));
    of.write(reinterpret_cast<char*>(tabIndex.data()), sizeof(tabIndex[0])*tabLen);

    if (of.fail())
    {
        std::cout << "Error writing tabIndex " << tabLen << "\n";
    }
    // store kmers
    size_t kmerNum = kmerTable.size();
    of.write(reinterpret_cast<char*>(&kmerNum), sizeof(kmerNum));
    for (size_t i = 0; i < kmerNum; ++i)
    {
        const uint32_t k = KMER::getCore(kmerTable[i]);
		const uint32_t kTMask = getTMask(kmerTable[i]);
        of.write(reinterpret_cast<const char*>(&k),sizeof(uint32_t));
        of.write(reinterpret_cast<const char*>(&kTMask),sizeof(uint32_t));
    }

    if (of.fail())
    {
        std::cout << "Error writing kmerTable " << kmerNum << "\n";
    }
    // store strands
    std::ofstream strandFile(filepath + "_strands", std::ofstream::binary | std::ofstream::trunc);
    write_strands(strandFile);
    if (strandFile.fail())
    {
        std::cerr << "Error while writing index strand file! Terminating...\n\n";
        exit(1);
    }
    strandFile.close();

    // store meta CpGs
    size_t metaCpGNum = metaCpGs.size();
    of.write(reinterpret_cast<char*>(&metaCpGNum), sizeof(metaCpGNum));
    of.write(reinterpret_cast<char*>(metaCpGs.data()), sizeof(metaCpGs[0]) * metaCpGNum);
    if (of.fail())
    {
        std::cout << "Error while writing metaCpGs " << metaCpGNum << "\n";
    }

    // store start meta CpGs
    metaCpGNum = metaStartCpGs.size();
    of.write(reinterpret_cast<char*>(&metaCpGNum), sizeof(metaCpGNum));
    of.write(reinterpret_cast<char*>(metaStartCpGs.data()), sizeof(metaStartCpGs[0]) * metaCpGNum);
    if (of.fail())
    {
        std::cout << "Error while writing start metaCpGs " << metaCpGNum << "\n";
    }
    // store filtered kmers
    write_filteredKmers(of);
    if (of.fail())
    {
        std::cout << "Error while writing filtered kmer map\n";
    }
    if (of.fail())
    {
        std::cerr << "Error while writing index file... Terminating\n\n";
        exit(1);
    }

	// write chromosome ID mapping
	chromNum = chrMap.size();
	of.write(reinterpret_cast<char*>(&chromNum), sizeof(chromNum));
	for (auto chrMapping : chrMap)
	{
		uint8_t internID = chrMapping.first;
		of.write(reinterpret_cast<char*>(&internID), sizeof(internID));
		size_t strlen = chrMapping.second.size();
		of.write(reinterpret_cast<char*>(&strlen), sizeof(strlen));
		for (char c : chrMapping.second)
		{
			of.write(&c, sizeof(c));
		}
	}


    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    std::cout << "Finished writing index to file " << filepath << " in " << runtime << "s\n\n";
    of.close();
}

void RefGenome::load(const std::string& filepath)
{

    std::cout << "Start reading index file " << filepath << "\n";
    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    // open file in binary mode
    std::ifstream ifs(filepath, std::ifstream::binary);

    // load CONSTANTS
    auto htabs = MyConst::HTABSIZE;
    ifs.read(reinterpret_cast<char*>(&htabs), sizeof(htabs));
    if (htabs != MyConst::HTABSIZE)
    {
        std::cerr << "Hash table size in source code and index file are different!\n\n";
        exit(1);
    }
    auto readl = MyConst::READLEN;
    ifs.read(reinterpret_cast<char*>(&readl), sizeof(readl));
    if (readl != MyConst::READLEN)
    {
        std::cerr << "Read length used in source code and index file are different!\n\n";
        exit(1);
    }
    auto winl = MyConst::WINLEN;
    ifs.read(reinterpret_cast<char*>(&winl), sizeof(winl));
    if (winl != MyConst::WINLEN)
    {
        std::cerr << "Meta CpG length used in source code and index file are different!\n\n";
        exit(1);
    }
    auto kmerl = MyConst::KMERLEN;
    ifs.read(reinterpret_cast<char*>(&kmerl), sizeof(kmerl));
    if (kmerl != MyConst::KMERLEN)
    {
        std::cerr << "k-mer length used in source code and index file are different!\n\n";
        exit(1);
    }
    auto kmerc = MyConst::KMERCUTOFF;
    ifs.read(reinterpret_cast<char*>(&kmerc), sizeof(kmerc));
    if (kmerc != MyConst::KMERCUTOFF)
    {
        std::cerr << "k-mer cutoff used in source code and index file are different!\n\n";
        exit(1);
    }

    // load NON-start CpGs
    size_t cpgNum;
    ifs.read(reinterpret_cast<char*>(&cpgNum), sizeof(cpgNum));
    cpgTable.resize(cpgNum);
    ifs.read(reinterpret_cast<char*>(cpgTable.data()), sizeof(cpgTable[0]) * cpgNum);
    // load start CpGs
    size_t cpgStartNum;
    ifs.read(reinterpret_cast<char*>(&cpgStartNum), sizeof(cpgStartNum));
    cpgStartTable.resize(cpgStartNum);
    ifs.read(reinterpret_cast<char*>(cpgStartTable.data()), sizeof(cpgStartTable[0]) * cpgStartNum);
    // load reference sequence
    size_t chromNum;
    ifs.read(reinterpret_cast<char*>(&chromNum), sizeof(chromNum));
    fullSeq.resize(chromNum);
    for (size_t i = 0; i < chromNum; ++i)
    {
        size_t chromLen;
        ifs.read(reinterpret_cast<char*>(&chromLen), sizeof(chromLen));
        fullSeq[i].resize(chromLen);
        ifs.read(fullSeq[i].data(), chromLen);
    }
    // tabIndex
    size_t tabLen;
    ifs.read(reinterpret_cast<char*>(&tabLen), sizeof(tabLen));
    tabIndex.resize(tabLen);
    ifs.read(reinterpret_cast<char*>(tabIndex.data()), sizeof(tabIndex[0])*tabLen);

    // load kmers
    size_t kmerNum;
    ifs.read(reinterpret_cast<char*>(&kmerNum), sizeof(kmerNum));
    kmerTableSmall.resize(kmerNum);
    for (size_t i = 0; i < kmerNum; ++i)
    {
		uint32_t metaCpG;
		uint32_t tMask;
        ifs.read(reinterpret_cast<char*>(&(metaCpG)),sizeof(uint32_t));
        ifs.read(reinterpret_cast<char*>(&(tMask)),sizeof(uint32_t));
		kmerTableSmall[i] = KMER_S::constructKmerS(metaCpG, tMask);
    }

    // load strands
    std::ifstream ifs_bool(filepath + "_strands", std::ifstream::binary);
    read_strands(ifs_bool);
    if (ifs_bool.fail())
    {
        std::cerr << "Error while reading index strand file! Terminating...\n\n";
        exit(1);
    }

    ifs_bool.close();

    // load meta CpGs
    size_t metaCpGNum;
    ifs.read(reinterpret_cast<char*>(&metaCpGNum), sizeof(metaCpGNum));
    metaCpGs.resize(metaCpGNum);
    ifs.read(reinterpret_cast<char*>(metaCpGs.data()), sizeof(metaCpGs[0]) * metaCpGNum);

    // load start meta CpGs
    size_t metaStartCpGNum;
    ifs.read(reinterpret_cast<char*>(&metaStartCpGNum), sizeof(metaStartCpGNum));
    metaStartCpGs.resize(metaStartCpGNum);
    ifs.read(reinterpret_cast<char*>(metaStartCpGs.data()), sizeof(metaStartCpGs[0]) * metaStartCpGNum);

    // load filtered kmers
    read_filteredKmers(ifs);
    if (ifs.fail())
    {
        std::cerr << "Error while reading index file! Terminating...\n\n";
        exit(1);
    }
	// load chromosome ID mapping
	ifs.read(reinterpret_cast<char*>(&chromNum), sizeof(chromNum));
	for (size_t i = 0; i < chromNum; ++i)
	{
		uint8_t internID;
		ifs.read(reinterpret_cast<char*>(&internID), sizeof(internID));
		std::string externID;
		size_t strlen;
		ifs.read(reinterpret_cast<char*>(&strlen), sizeof(strlen));
		for (size_t strL = 0; strL < strlen; ++strL)
		{
			char c;
			ifs.read(&c, sizeof(c));
			externID += c;
		}
		chrMap.insert(std::pair<uint8_t, std::string>(internID, externID));
	}
    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    std::cout << "Finished reading index file " << filepath << " in " << runtime << "s\n\n";
    ifs.close();
}

// General idea taken from
// https://stackoverflow.com/questions/29623605/how-to-dump-stdvectorbool-in-a-binary-file
// accessed on 21th July 2017
inline void RefGenome::write_strands(std::ofstream& of)
{

    size_t n = strandTable.size();
    of.write(reinterpret_cast<char*>(&n), sizeof(n));

    // go over all bits
    for (size_t i = 0; i < n;)
    {
        // saves 8 flags in one char
        unsigned char mergedBits = 0;
        for (unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1)
        {
            if (strandTable[i])
            {
                mergedBits |= mask;
            }

        }
        of.write(reinterpret_cast<const char*>(&mergedBits), sizeof(mergedBits));
    }
}
inline void RefGenome::read_strands(std::ifstream& ifs)
{
    size_t n;
    ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
    strandTable.resize(n);

    // go over all bits
    for (size_t i = 0; i < n;)
    {
        // will hold 8 bits
        unsigned char mergedBits;
        // read one chunk
        ifs.read(reinterpret_cast<char*>(&mergedBits), sizeof(mergedBits));
        // store bits of chunk
        for(unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1)
        {
            strandTable[i] = mergedBits & mask;
        }
    }
}
inline void RefGenome::write_filteredKmers(std::ofstream& of)
{

    size_t n = filteredKmers.size();
    of.write(reinterpret_cast<char*>(&n), sizeof(n));
    for (const uint64_t k : filteredKmers)
    {

        of.write(reinterpret_cast<const char*>(&k), sizeof(k));
    }
}
inline void RefGenome::read_filteredKmers(std::ifstream& ifs)
{

    size_t n;
    ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
    for (size_t i = 0; i < n; ++i)
    {

        uint64_t k;
        ifs.read(reinterpret_cast<char*>(&k), sizeof(k));
        filteredKmers.emplace(k);
    }

}




void RefGenome::generateMetaCpGs()
{

    uint32_t cpgStartInd = 0;
    // current chromosome we are looking at
    unsigned int currChr = 0;
    // start and end index for CpGs in current window
    uint32_t start = 0;
    // Assuming window is larger then readlength
    while (cpgStartInd < cpgStartTable.size())
    {
		bool isMod = false;

        while (cpgStartInd < cpgStartTable.size() && cpgStartTable[cpgStartInd].chrom == currChr)
        {
			isMod = true;
            ++cpgStartInd;
        }

        if (isMod)
			metaStartCpGs.push_back({start, cpgStartInd - 1});
        ++currChr;
        start = cpgStartInd;
    }


    uint32_t cpgTabInd = 0;
    start = 0;
    currChr = 0;

    // while there are CpGs left to look at
    //
    while (cpgTabInd < cpgTable.size())
    {

        // fill current metaCpG window with CpGs
        const uint32_t windowEnd = cpgTable[start].pos + MyConst::WINLEN;

        for (; cpgTabInd < cpgTable.size(); ++cpgTabInd)
        {
            if (cpgTable[cpgTabInd].pos >= windowEnd || cpgTable[cpgTabInd].chrom != currChr)
            {

                // generate meta CpG
                metaCpGs.push_back({start, cpgTabInd - 1});
                // set start for next meta CpG
                start = cpgTabInd;
                // update current chromosome index if necessary
                if (cpgTable[cpgTabInd].chrom != currChr) ++currChr;
                // break out of inner for loop (i.e. filling of current metaCpG
                break;
            }

        }

    }
    if (cpgTable.size() > 0)
    {
        // put in meta CpG containging last CpGs
        metaCpGs.push_back({start, cpgTabInd - 1});
    }

    metaCpGs.shrink_to_fit();
    metaStartCpGs.shrink_to_fit();
}


void RefGenome::generateHashes(std::vector<std::vector<char> >& genomeSeq)
{


    estimateTablesizes(genomeSeq);

    std::cout << "\nKmer table size: " << kmerTable.size() << std::endl;
    std::cout << "\nMeta CpGs: " << metaCpGs.size() << std::endl;

    // hash CpGs from the start
    uint32_t cpgCount = 0;
    for (std::vector<struct metaCpG>::const_iterator it = metaStartCpGs.begin(); it != metaStartCpGs.end(); ++it, ++cpgCount)
    {

        unsigned int lastPos = 0;
        for (uint32_t cpgInd = it->start; cpgInd <= it->end; ++cpgInd)
        {

            ntHashFirst(genomeSeq[cpgStartTable[it->start].chrom], lastPos, cpgStartTable[cpgInd].pos, cpgCount);
        }
    }

    cpgCount = 0;
    // hash each CpG
    for (std::vector<struct metaCpG>::const_iterator it = metaCpGs.begin(); it != metaCpGs.end(); ++it, ++cpgCount)
    {

        unsigned int lastPos = 0;

        for (uint32_t cpgInd = it->start; cpgInd <= it->end; ++cpgInd)
        {

            // how long is the rest of the sequence after the end of last cpg in meta cpg
            const unsigned int remainderBps = genomeSeq[cpgTable[cpgInd].chrom].size() - (cpgTable[cpgInd].pos + MyConst::READLEN);
            // if we can read the full sequence length after the CpG
            if (remainderBps >= (MyConst::READLEN - 2) )
            {

                ntHashChunk(genomeSeq[cpgTable[it->start].chrom], lastPos, cpgTable[cpgInd].pos, cpgCount, cpgTable[cpgInd].pos - cpgTable[it->start].pos);


            } else {

                ntHashLast(genomeSeq[cpgTable[it->start].chrom], lastPos, cpgTable[cpgInd].pos, remainderBps, cpgCount, cpgTable[cpgInd].pos - cpgTable[it->start].pos);
                // if we have the case that we hit the end of the sequence, we can be sure to hash everything in this run and hence can skip the remaining CpGs
                break;
            }
        }
    }
}






inline void RefGenome::ntHashChunk(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const uint32_t& metacpg, const uint32_t&& metaOff)
{

	// kmers to be skipped
	const int skipKmer = std::max((int)lastPos - (int)pos, (int)0);
	// construct corresponding sequence with reduced alphabet
	std::vector<char> redSeq(2*MyConst::READLEN - 2);
	std::vector<char> redSeqRev(2*MyConst::READLEN - 2);

	std::vector<char>::const_iterator start = seq.begin();
	std::vector<char>::const_iterator end = seq.begin();

	// last N before the CpG or last position that was looked at
	int lasN = skipKmer - 1;
	unsigned int j;

	// if we need to hash kmer starting at the left of current CpG
	if (lastPos < pos + MyConst::READLEN - MyConst::KMERLEN)
	{
		std::advance(start, pos + skipKmer);
		std::advance(end, pos + MyConst::READLEN);

		// position in substring
		j = skipKmer;
		// move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
		for ( ; start != end; ++start, ++j)
		{

			const int revPos = 2*MyConst::READLEN - 3 - j;
			switch (*start)
			{
				case 'N':
					lasN = j;
					break;

				case 'C':
					redSeq[j] = 'T';
					redSeqRev[revPos] = 'G';
					break;

				case 'G':
					redSeq[j] = 'G';
					redSeqRev[revPos] = 'T';
					break;

				case 'T':
					redSeq[j] = 'T';
					redSeqRev[revPos] = 'A';
					break;

				case 'A':
					redSeq[j] = 'A';
					redSeqRev[revPos] = 'T';
					break;

				default:
					std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

			}

		}
		// move to one after last N
		++lasN;
		// reassign end to the position of G or one after last hashed kmer
		end = seq.begin() + pos + MyConst::READLEN - 1;
	} else {

		++lasN;
		end = seq.begin() + pos + skipKmer - 1;
	}
	// move over second half in reverse order
	// reassign current position
	j = 2*MyConst::READLEN - 3;
	// reassign start to final position
	start = seq.begin() + pos + j;

	// offset where the first N after the CpG is
	int off = 2*MyConst::READLEN - 2;
	for ( ; start != end; --start, --j)
	{
		const int revPos = 2*MyConst::READLEN - 3 - j;
		switch (*start)
		{
			case 'N':
				off = j;
				break;

			case 'C':
				redSeq[j] = 'T';
				redSeqRev[revPos] = 'G';
				break;

			case 'G':
				redSeq[j] = 'G';
				redSeqRev[revPos] = 'T';
				break;

			case 'T':
				redSeq[j] = 'T';
				redSeqRev[revPos] = 'A';
				break;

			case 'A':
				redSeq[j] = 'A';
				redSeqRev[revPos] = 'T';
				break;

			default:
				std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
		}


	}

	// length of the context around the CpG without Ns
	unsigned int contextLen = off - lasN;
	// if we don't have enough to read, return without hashing
	if (contextLen < MyConst::KMERLEN)
	{
		return;
	}

	const char* seqStart = redSeq.data() + lasN;
	const char* seqStartRev = redSeqRev.data() + (2*MyConst::READLEN - 2 - off);

	// initial hash backward
	uint64_t rhVal;
	uint64_t srVal = ntHash::NTPS64(seqStartRev, MyConst::SEED, MyConst::KMERLEN, rhVal);
	// first kmer on reverse complement corresponds to last kmer in forward sequence
	uint64_t kPosRev = off - MyConst::KMERLEN + metaOff;

	// update kmer table
	kmerTable[--tabIndex[srVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev));
	strandTable[tabIndex[srVal % MyConst::HTABSIZE]] = false;


	// hash kmers of backward strand
	for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
	{
		--kPosRev;
		srVal = ntHash::NTPS64(seqStartRev+i+1, MyConst::SEED, seqStartRev[i], seqStartRev[i + MyConst::KMERLEN], MyConst::KMERLEN, rhVal);
		// update kmer table
		kmerTable[--tabIndex[srVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev));
		strandTable[tabIndex[srVal % MyConst::HTABSIZE]] = false;
	}

	// initial hash forward
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seqStart, MyConst::SEED, MyConst::KMERLEN, fhVal);
	uint64_t kPos = lasN + metaOff;

	// update kmer table
	kmerTable[--tabIndex[sfVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos));
	strandTable[tabIndex[sfVal % MyConst::HTABSIZE]] = true;

	// hash kmers of forward strand
	for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
	{
		++kPos;
		sfVal = ntHash::NTPS64(seqStart+i+1, MyConst::SEED, seqStart[i], seqStart[i + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
		// update kmer table
		kmerTable[--tabIndex[sfVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos));
		strandTable[tabIndex[sfVal % MyConst::HTABSIZE]] = true;
	}
	lastPos = pos + off - MyConst::KMERLEN + 1;

}





void RefGenome::ntHashLast(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG, const uint32_t& metacpg, uint32_t&& metaOff)
{
    // kmers to be skipped
	const int skipKmer = std::max((int)lastPos - (int)pos, (int)0);
    // construct corresponding sequence with reduced alphabet
    std::vector<char> redSeq(MyConst::READLEN + bpsAfterCpG);
    std::vector<char> redSeqRev(MyConst::READLEN + bpsAfterCpG);

    std::vector<char>::const_iterator start = seq.begin();
    std::vector<char>::const_iterator end = seq.begin();

    // last N before the CpG or last position that was looked at
    int lasN = skipKmer - 1;
    unsigned int j;

    // if we need to hash kmer starting at the left of current CpG
    if (lastPos < pos + MyConst::READLEN - MyConst::KMERLEN)
    {
        std::advance(start, pos + skipKmer);
        std::advance(end, pos + MyConst::READLEN);

        // position in substring
        j = skipKmer;
        // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
        for ( ; start != end; ++start, ++j)
        {

            const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
            switch (*start)
            {
                case 'N':
                    lasN = j;
                    break;

                case 'C':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'G';
                    break;

                case 'G':
                    redSeq[j] = 'G';
                    redSeqRev[revPos] = 'T';
                    break;

                case 'T':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'A';
                    break;

                case 'A':
                    redSeq[j] = 'A';
                    redSeqRev[revPos] = 'T';
                    break;

                default:
                    std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

            }

        }
        // move to one after last N
        ++lasN;
        // reassign end to the position of G or one after last hashed kmer
        end = seq.begin();
        std::advance(end, pos + MyConst::READLEN - 1);
    } else {

        ++lasN;
        end = seq.begin();
        std::advance(end, pos + skipKmer - 1);
    }
    // reassign current position
    j = MyConst::READLEN + bpsAfterCpG - 1;
    // move over second half in reverse order
    // reassign start to final position
    start = seq.begin();
    std::advance(start, pos + j);

    // offset where the first N after the CpG is
    int off = MyConst::READLEN + bpsAfterCpG;
    for ( ; start != end; --start, --j)
    {
        const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
        switch (*start)
        {
            case 'N':
                off = j;
                break;

            case 'C':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[j] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[j] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
        }


    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = off - lasN;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }

    const char* seqStart = redSeq.data() + lasN;
    const char* seqStartRev = redSeqRev.data() + (MyConst::READLEN + bpsAfterCpG - off);

    // initial hash backward
	uint64_t rhVal;
	uint64_t srVal = ntHash::NTPS64(seqStartRev, MyConst::SEED, MyConst::KMERLEN, rhVal);
    // first kmer on reverse complement corresponds to last kmer in forward sequence
    uint64_t kPosRev = off - MyConst::KMERLEN;

    // update kmer table
    kmerTable[--tabIndex[srVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
    strandTable[tabIndex[srVal % MyConst::HTABSIZE]] = false;


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        --kPosRev;
		srVal = ntHash::NTPS64(seqStartRev+i+1, MyConst::SEED, seqStartRev[i], seqStartRev[i + MyConst::KMERLEN], MyConst::KMERLEN, rhVal);
        // update kmer table
        kmerTable[--tabIndex[srVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
        strandTable[tabIndex[srVal % MyConst::HTABSIZE]] = false;
    }

    // initial hash forward
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seqStart, MyConst::SEED, MyConst::KMERLEN, fhVal);
    uint64_t kPos = lasN;

    // update kmer table
    kmerTable[--tabIndex[sfVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
    strandTable[tabIndex[sfVal % MyConst::HTABSIZE]] = true;

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ++kPos;
		sfVal = ntHash::NTPS64(seqStart+i+1, MyConst::SEED, seqStart[i], seqStart[i + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
        // update kmer table
        kmerTable[--tabIndex[sfVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
        strandTable[tabIndex[sfVal % MyConst::HTABSIZE]] = true;
    }
}

void RefGenome::ntHashFirst(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& posOfCpG, const uint32_t& metacpg)
{
    // kmers to be skipped
    const int skipKmer = lastPos;
    // construct corresponding sequence with reduced alphabet
    std::vector<char> redSeq(MyConst::READLEN + posOfCpG);
    std::vector<char> redSeqRev(MyConst::READLEN + posOfCpG);

    std::vector<char>::const_iterator start = seq.begin();
    std::vector<char>::const_iterator end = seq.begin();

    // last N before the CpG or last position that was looked at
    int lasN = skipKmer - 1;
    unsigned int j;

    // if we need to hash kmer starting at the left of current CpG
    if (lastPos < static_cast<long>(posOfCpG) + 2 - MyConst::KMERLEN)
    {
        std::advance(start, skipKmer);
        std::advance(end, posOfCpG + 2);

        // position in substring
        j = skipKmer;
        // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
        for ( ; start != end; ++start, ++j)
        {

            const int revPos = MyConst::READLEN + posOfCpG - 1 - j;
            switch (*start)
            {
                case 'N':
                    lasN = j;
                    break;

                case 'C':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'G';
                    break;

                case 'G':
                    redSeq[j] = 'G';
                    redSeqRev[revPos] = 'T';
                    break;

                case 'T':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'A';
                    break;

                case 'A':
                    redSeq[j] = 'A';
                    redSeqRev[revPos] = 'T';
                    break;

                default:
                    std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

            }

        }
        // move to one after last N
        ++lasN;
        // reassign end to the position of G or one after last hashed kmer
        end = seq.begin();
        std::advance(end, posOfCpG - 1);
    } else {

        ++lasN;
        end = seq.begin();
        std::advance(end, skipKmer - 1);
    }
    // reassign current position
    j = MyConst::READLEN + posOfCpG - 1;
    // move over second half in reverse order
    // reassign start to final position
    start = seq.begin();
    std::advance(start, j);

    // offset where the first N after the CpG is
    int off = MyConst::READLEN + posOfCpG;
    for ( ; start != end; --start, --j)
    {
        const int revPos = MyConst::READLEN + posOfCpG - 1 - j;
        switch (*start)
        {
            case 'N':
                off = j;
                break;

            case 'C':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[j] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[j] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
        }


    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = off - lasN;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }

    const char* seqStart = redSeq.data() + lasN;
    const char* seqStartRev = redSeqRev.data() + (MyConst::READLEN + posOfCpG - off);

    // initial hash backward
	uint64_t rhVal;
	uint64_t srVal = ntHash::NTPS64(seqStartRev, MyConst::SEED, MyConst::KMERLEN, rhVal);
    // first kmer on reverse complement corresponds to last kmer in forward sequence
    uint64_t kPosRev = off - MyConst::KMERLEN;

    // update kmer table
    kmerTable[--tabIndex[srVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPosRev));
    strandTable[tabIndex[srVal % MyConst::HTABSIZE]] = false;


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        --kPosRev;
		srVal = ntHash::NTPS64(seqStartRev+i+1, MyConst::SEED, seqStartRev[i], seqStartRev[i + MyConst::KMERLEN], MyConst::KMERLEN, rhVal);
        // update kmer table
        kmerTable[--tabIndex[srVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPosRev));
        strandTable[tabIndex[srVal % MyConst::HTABSIZE]] = false;
    }

    // initial hash forward
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seqStart, MyConst::SEED, MyConst::KMERLEN, fhVal);
    uint64_t kPos = lasN;

    // update kmer table
    kmerTable[--tabIndex[sfVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPos));
    strandTable[tabIndex[sfVal % MyConst::HTABSIZE]] = true;

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ++kPos;
		sfVal = ntHash::NTPS64(seqStart+i+1, MyConst::SEED, seqStart[i], seqStart[i + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
        // update kmer table
        kmerTable[--tabIndex[sfVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPos));
        strandTable[tabIndex[sfVal % MyConst::HTABSIZE]] = true;
    }
    lastPos = off - MyConst::KMERLEN + 1;

}





inline void RefGenome::ntCountChunk(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos)
{

	// kmers to be skipped
	const int skipKmer = std::max((int)lastPos - (int)pos, (int)0);
	// construct corresponding sequence with reduced alphabet
	std::vector<char> redSeq(2*MyConst::READLEN - 2);
	std::vector<char> redSeqRev(2*MyConst::READLEN - 2);

	std::vector<char>::const_iterator start = seq.begin();
	std::vector<char>::const_iterator end = seq.begin();

	// last N before the CpG or last position that was looked at
	int lasN = skipKmer - 1;
	unsigned int j;

	// if we need to hash kmer starting at the left of current CpG
	if (lastPos < pos + MyConst::READLEN - MyConst::KMERLEN)
	{
		std::advance(start, pos + skipKmer);
		std::advance(end, pos + MyConst::READLEN);

		// position in substring
		j = skipKmer;
		// move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
		for ( ; start != end; ++start, ++j)
		{

			const int revPos = 2*MyConst::READLEN - 3 - j;
			switch (*start)
			{
				case 'N':
					lasN = j;
					break;

				case 'C':
					redSeq[j] = 'T';
					redSeqRev[revPos] = 'G';
					break;

				case 'G':
					redSeq[j] = 'G';
					redSeqRev[revPos] = 'T';
					break;

				case 'T':
					redSeq[j] = 'T';
					redSeqRev[revPos] = 'A';
					break;

				case 'A':
					redSeq[j] = 'A';
					redSeqRev[revPos] = 'T';
					break;

				default:
					std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

			}

		}
		// reassign end to the position of G or one after last hashed kmer
		end = seq.begin();
		std::advance(end, pos + MyConst::READLEN - 1);
	} else {

		end = seq.begin();
		std::advance(end, pos + skipKmer - 1);
	}
	// move to one after last N/ one after last hashed kmer
	++lasN;
	// move over second half in reverse order
	// reassign start to final position
	start = seq.begin();
	std::advance(start, pos + 2*MyConst::READLEN - 3);

	// reassign current position
	j = 2*MyConst::READLEN - 3;
	// offset where the first N after the CpG start is
	int off = 2*MyConst::READLEN - 2;
	for ( ; start != end; --start, --j)
	{
		const int revPos = 2*MyConst::READLEN - 3 - j;
		switch (*start)
		{
			case 'N':
				off = j;
				break;

			case 'C':
				redSeq[j] = 'T';
				redSeqRev[revPos] = 'G';
				break;

			case 'G':
				redSeq[j] = 'G';
				redSeqRev[revPos] = 'T';
				break;

			case 'T':
				redSeq[j] = 'T';
				redSeqRev[revPos] = 'A';
				break;

			case 'A':
				redSeq[j] = 'A';
				redSeqRev[revPos] = 'T';
				break;

			default:
				std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
		}


	}

	// length of the context around the CpG that should be hashed
	unsigned int contextLen = off - lasN;
	// if we don't have enough to read, return without hashing
	if (contextLen < MyConst::KMERLEN)
	{
		return;
	}

	const char* seqStart = redSeq.data() + lasN;
	const char* seqStartRev = redSeqRev.data() + (2*MyConst::READLEN - 2 - off);

	// initial hash backward
	uint64_t rhVal;
	uint64_t srVal = ntHash::NTPS64(seqStartRev, MyConst::SEED, MyConst::KMERLEN, rhVal);

	// update indices
	++tabIndex[srVal % MyConst::HTABSIZE];

	// hash kmers of backward strand
	for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
	{
		srVal = ntHash::NTPS64(seqStartRev+i+1, MyConst::SEED, seqStartRev[i], seqStartRev[i + MyConst::KMERLEN], MyConst::KMERLEN, rhVal);
		// update indices
		++tabIndex[srVal % MyConst::HTABSIZE];

	}

	// initial hash forward
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seqStart, MyConst::SEED, MyConst::KMERLEN, fhVal);

	// update indices
	++tabIndex[sfVal % MyConst::HTABSIZE];

	// hash kmers of forward strand
	for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
	{

		sfVal = ntHash::NTPS64(seqStart+i+1, MyConst::SEED, seqStart[i], seqStart[i + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
		// update indices
		++tabIndex[sfVal % MyConst::HTABSIZE];
	}
	// update position of last hashed kmer (+ 1)
	lastPos = pos + off - MyConst::KMERLEN + 1;

}





void RefGenome::ntCountFirst(std::vector<char>& seq, uint32_t& lastPos, const unsigned int& posOfCpG)
{

    // kmers to be skipped
    const int skipKmer = lastPos;
    // construct corresponding sequence with reduced alphabet
    std::vector<char> redSeq(MyConst::READLEN + posOfCpG);
    std::vector<char> redSeqRev(MyConst::READLEN + posOfCpG);

    std::vector<char>::const_iterator start = seq.begin();
    std::vector<char>::const_iterator end = seq.begin();

    // last N before the CpG or last position that was looked at
    int lasN = skipKmer - 1;
    unsigned int j;

    // if we need to hash kmer starting at the left of current CpG
    if (lastPos < static_cast<long>(posOfCpG) + 2 - MyConst::KMERLEN)
    {
        std::advance(start, skipKmer);
        std::advance(end, posOfCpG + 2);

        // position in substring
        j = skipKmer;
        // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
        for ( ; start != end; ++start, ++j)
        {

            const int revPos = MyConst::READLEN + posOfCpG - 1 - j;
            switch (*start)
            {
                case 'N':
                    lasN = j;
                    break;

                case 'C':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'G';
                    break;

                case 'G':
                    redSeq[j] = 'G';
                    redSeqRev[revPos] = 'T';
                    break;

                case 'T':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'A';
                    break;

                case 'A':
                    redSeq[j] = 'A';
                    redSeqRev[revPos] = 'T';
                    break;

                default:
                    std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

            }

        }
        // move to one after last N
        ++lasN;
        // reassign end to the position of G or one after last hashed kmer
        end = seq.begin();
        std::advance(end, posOfCpG - 1);
    } else {

        ++lasN;
        end = seq.begin();
        std::advance(end, skipKmer - 1);
    }
    // reassign current position
    j = MyConst::READLEN + posOfCpG - 1;
    // move over second half in reverse order
    // reassign start to final position
    start = seq.begin();
    std::advance(start, j);

    // offset where the first N after the CpG is
    int off = MyConst::READLEN + posOfCpG;
    for ( ; start != end; --start, --j)
    {
        const int revPos = MyConst::READLEN + posOfCpG - 1 - j;
        switch (*start)
        {
            case 'N':
                off = j;
                break;

            case 'C':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[j] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[j] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
        }


    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = off - lasN;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }

    const char* seqStart = redSeq.data() + lasN;
    const char* seqStartRev = redSeqRev.data() + (MyConst::READLEN + posOfCpG - off);

    // initial hash backward
	uint64_t rhVal;
	uint64_t srVal = ntHash::NTPS64(seqStartRev, MyConst::SEED, MyConst::KMERLEN, rhVal);

    // update kmer table
    ++tabIndex[srVal % MyConst::HTABSIZE];


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
		srVal = ntHash::NTPS64(seqStartRev+i+1, MyConst::SEED, seqStartRev[i], seqStartRev[i + MyConst::KMERLEN], MyConst::KMERLEN, rhVal);
        ++tabIndex[srVal % MyConst::HTABSIZE];
    }

    // initial hash forward
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seqStart, MyConst::SEED, MyConst::KMERLEN, fhVal);

    ++tabIndex[sfVal % MyConst::HTABSIZE];

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
		sfVal = ntHash::NTPS64(seqStart+i+1, MyConst::SEED, seqStart[i], seqStart[i + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
        ++tabIndex[sfVal % MyConst::HTABSIZE];
    }
    lastPos = off - MyConst::KMERLEN + 1;

}






void RefGenome::ntCountLast(std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG)
{
    // kmers to be skipped
	const int skipKmer = std::max((int)lastPos - (int)pos, (int)0);
    // construct corresponding sequence with reduced alphabet
    std::vector<char> redSeq(MyConst::READLEN + bpsAfterCpG);
    std::vector<char> redSeqRev(MyConst::READLEN + bpsAfterCpG);

    std::vector<char>::const_iterator start = seq.begin();
    std::vector<char>::const_iterator end = seq.begin();

    // last N before the CpG or last position that was looked at
    int lasN = skipKmer - 1;
    unsigned int j;

    // if we need to hash kmer starting at the left of current CpG
    if (lastPos < pos + MyConst::READLEN - MyConst::KMERLEN)
    {
        std::advance(start, pos + skipKmer);
        std::advance(end, pos + MyConst::READLEN);

        // position in substring
        j = skipKmer;
        // move over sequence until CpG, construct reduced alphabet string and retrieve the positions of the last N
        for ( ; start != end; ++start, ++j)
        {

            const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
            switch (*start)
            {
                case 'N':
                    lasN = j;
                    break;

                case 'C':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'G';
                    break;

                case 'G':
                    redSeq[j] = 'G';
                    redSeqRev[revPos] = 'T';
                    break;

                case 'T':
                    redSeq[j] = 'T';
                    redSeqRev[revPos] = 'A';
                    break;

                case 'A':
                    redSeq[j] = 'A';
                    redSeqRev[revPos] = 'T';
                    break;

                default:
                    std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";

            }

        }
        // move to one after last N
        ++lasN;
        // reassign end to the position of G or one after last hashed kmer
        end = seq.begin();
        std::advance(end, pos + MyConst::READLEN - 1);
    } else {

        ++lasN;
        end = seq.begin();
        std::advance(end, pos + skipKmer - 1);
    }
    // reassign current position
    j = MyConst::READLEN + bpsAfterCpG - 1;
    // move over second half in reverse order
    // reassign start to final position
    start = seq.begin();
    std::advance(start, pos + j);

    // offset where the first N after the CpG is
    int off = MyConst::READLEN + bpsAfterCpG;
    for ( ; start != end; --start, --j)
    {
        const int revPos = MyConst::READLEN + bpsAfterCpG - 1 - j;
        switch (*start)
        {
            case 'N':
                off = j;
                break;

            case 'C':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'G';
                break;

            case 'G':
                redSeq[j] = 'G';
                redSeqRev[revPos] = 'T';
                break;

            case 'T':
                redSeq[j] = 'T';
                redSeqRev[revPos] = 'A';
                break;

            case 'A':
                redSeq[j] = 'A';
                redSeqRev[revPos] = 'T';
                break;

            default:
                std::cerr << "Read has unknown character \'" << *start << "\' result will not be reliable.\n\n";
        }


    }

    // length of the context around the CpG without Ns
    unsigned int contextLen = off - lasN;
    // if we don't have enough to read, return without hashing
    if (contextLen < MyConst::KMERLEN)
    {
        return;
    }

    const char* seqStart = redSeq.data() + lasN;
    const char* seqStartRev = redSeqRev.data() + (MyConst::READLEN + bpsAfterCpG - off);

    // initial hash backward
	uint64_t rhVal;
	uint64_t srVal = ntHash::NTPS64(seqStartRev, MyConst::SEED, MyConst::KMERLEN, rhVal);

    // update indices
    ++tabIndex[srVal % MyConst::HTABSIZE];


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
		srVal = ntHash::NTPS64(seqStartRev+i+1, MyConst::SEED, seqStartRev[i], seqStartRev[i + MyConst::KMERLEN], MyConst::KMERLEN, rhVal);
        // update indices
        ++tabIndex[srVal % MyConst::HTABSIZE];
    }

    // initial hash forward
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seqStart, MyConst::SEED, MyConst::KMERLEN, fhVal);

    // update indices
    ++tabIndex[sfVal % MyConst::HTABSIZE];

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
		sfVal = ntHash::NTPS64(seqStart+i+1, MyConst::SEED, seqStart[i], seqStart[i + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
        // update indices
        ++tabIndex[sfVal % MyConst::HTABSIZE];
    }
}


void RefGenome::estimateTablesizes(std::vector<std::vector<char> >& genomeSeq)
{

    // count start CpG kmers
    for (metaCpG& m : metaStartCpGs)
    {

        // we know that all of the CpGs at start will overlap
        const uint8_t chr = cpgStartTable[m.start].chrom;

        uint32_t lastPos = 0;

        for (uint32_t cpgInd = m.start; cpgInd <= m.end; ++cpgInd)
        {
            ntCountFirst(genomeSeq[chr], lastPos, cpgStartTable[cpgInd].pos);

        }
    }
    // count normal CpG kmers
    for (metaCpG& m : metaCpGs)
    {

        const uint8_t chr = cpgTable[m.start].chrom;

        uint32_t lastPos = 0;

        // how long is the rest of the sequence after the end of last cpg in meta cpg
        unsigned int remainderBps = genomeSeq[cpgTable[m.start].chrom].size() - (cpgTable[m.start].pos + MyConst::READLEN);

        // kmers of first CpG
        if (remainderBps >= (MyConst::READLEN - 2) )
        {
            ntCountChunk(genomeSeq[chr], lastPos, cpgTable[m.start].pos);

        } else {

            ntCountLast(genomeSeq[chr], lastPos, cpgTable[m.start].pos, remainderBps);
        }

        // consecutive CpG kmers
        for (uint32_t cpgInd = m.start + 1; cpgInd <= m.end; ++cpgInd)
        {

            remainderBps = genomeSeq[cpgTable[cpgInd].chrom].size() - (cpgTable[cpgInd].pos + MyConst::READLEN);
            // if we can read the full sequence breadth after the CpG
            if (remainderBps >= (MyConst::READLEN - 2) )
            {

                // count the collisions
                ntCountChunk(genomeSeq[chr], lastPos, cpgTable[cpgInd].pos);

            } else {

                // count the collisions and break
                ntCountLast(genomeSeq[chr], lastPos, cpgTable[cpgInd].pos, remainderBps);
                break;
            }
        }


    }

    uint64_t sum = 0;
    // update to sums of previous entrys
    for (unsigned int i = 0; i < MyConst::HTABSIZE; ++i)
    {

        sum += tabIndex[i];
        tabIndex[i] = sum;
    }

    // resize table to number of kmers that have to be hashed
    kmerTable.resize(sum);
    strandTable.resize(sum);

    // fill dummy value
    tabIndex[MyConst::HTABSIZE] = sum;
}





inline void RefGenome::blacklist(const unsigned int& KSliceStart, const unsigned int& KSliceEnd, std::unordered_map<uint64_t, unsigned int>& bl)
{

	// iterate over kmerTable in the specified range
	for (unsigned int i = KSliceStart; i < KSliceEnd; ++i)
	{

		// get the kmer at that position
		KMER::kmer& k = kmerTable[i];
		bool sFlag = strandTable[i];

		const uint64_t kHash = reproduceKmerSeq(k, sFlag);

		// try to put sequence in map - if already in, count up
		// NOTE:    hash function is perfect, hence no implicit collisions before putting it into hashmap
		auto insertion = bl.find(kHash);
		if (insertion == bl.end())
		{

			bl[kHash] = 1;

		} else {

			insertion->second += 1;
		}
	}
}






inline uint64_t RefGenome::reproduceKmerSeq(KMER::kmer& k, bool sFlag)
{

	uint32_t pos = 0;
	uint8_t chrom;
	// reproduce kmer sequence
	if (KMER::isStartCpG(k))
	{

		// get the start position and chromosome of surrounding meta cpg wrapped into a CpG
		const struct CpG& startCpG = cpgStartTable[metaStartCpGs[KMER::getMetaCpG(k)].start];
		chrom = startCpG.chrom;

	} else {

		// get the start position and chromosome of surrounding meta cpg wrapped into a CpG
		const struct CpG& startCpG = cpgTable[metaCpGs[KMER::getMetaCpG(k)].start];
		pos = startCpG.pos;
		chrom = startCpG.chrom;

	}

	// retrieve sequence
	//
	// get iterator to sequence start
	auto seqStart = fullSeq[chrom].begin() + KMER::getOffset(k) + pos;
	uint64_t kSeq = 0;

	// if from forward strand, just translate
	if (sFlag)
	{

		for (unsigned int i = 0; i < MyConst::KMERLEN; ++i)
		{
			kSeq = kSeq << 2;

			if (MyConst::SEED[i])
			{
				switch (seqStart[i])
				{

					case 'C':
					case 'T':

						kSeq += 3;
						break;

					// case 'C':
					//
					//     kSeq += 1;
					//     break;

					case 'G':

						kSeq += 2;
						break;

				// end switch
				}
			}

		// end forloop over sequence
		}

	// otherwise, build reverse complement
	} else {

		for (unsigned int i = MyConst::KMERLEN; i > 0; --i)
		{
			kSeq = kSeq << 2;

			if (MyConst::SEED[MyConst::KMERLEN - i])
			{
				switch (seqStart[i - 1])
				{

					case 'G':
					case 'A':

						kSeq += 3;
						break;

					// case 'G':
					//
					//     kSeq += 1;
					//     break;

					case 'C':

						kSeq += 2;
						break;

				// end switch
				}
			}

		// end forloop over sequence
		}
	}

	return kSeq;

}





void RefGenome::filterRedundancyInHashTable()
{

    std::cout << "\nHash table size before filter: " << kmerTable.size() << std::endl;
    std::chrono::high_resolution_clock::time_point filterStartTime = std::chrono::high_resolution_clock::now();

    // iterator to the element that we process
    auto srcItK = kmerTable.begin();
    auto srcItS = strandTable.begin();
    // iterator to the position one after the last inserted FILTERED element, always at most as far as srcIt
    auto filterItK = kmerTable.begin();
    auto filterItS = strandTable.begin();

    // iterate over hashTable
    for (uint64_t i = 0; i < MyConst::HTABSIZE; ++i)
    {

        // previous IDs and flags
        bool isStart = false;
        bool isFwd = false;
        uint64_t metaID = 0xffffffffffffffffULL;
        const uint64_t prevSizeK = filterItK - kmerTable.begin();

        // iterate through vector elements
        for (uint64_t j = tabIndex[i]; j < tabIndex[i+1]; ++j, ++srcItK, ++srcItS)
        {

            if (KMER::getMetaCpG(*srcItK) != metaID || *srcItS != isFwd || KMER::isStartCpG(*srcItK) != isStart)
            {

                isFwd = *srcItS;
                isStart = KMER::isStartCpG(*srcItK);
                metaID = KMER::getMetaCpG(*srcItK);
                *filterItK = *srcItK;
                *filterItS = *srcItS;
                ++filterItK;
                ++filterItS;
            }

        }
        // update the indexing structure for collisions
        tabIndex[i] = prevSizeK;
    }

    // update dummy value used for efficient indexing
    tabIndex[MyConst::HTABSIZE] = filterItK - kmerTable.begin();

    // shrink to new size
    kmerTable.resize(filterItK - kmerTable.begin());
    strandTable.resize(filterItS - strandTable.begin());

    std::chrono::high_resolution_clock::time_point filterEndTime = std::chrono::high_resolution_clock::now();

    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(filterEndTime - filterStartTime).count();


    std::cout << "Hash table size after filter (running " << runtime << "s): " << tabIndex[MyConst::HTABSIZE] << "\n\n";
}


void RefGenome::filterHashTable()
{

    std::cout << "\nHash table size before filter: " << kmerTable.size() << std::endl;
    std::chrono::high_resolution_clock::time_point filterStartTime = std::chrono::high_resolution_clock::now();

    // will hold the counts of individual kmers for each hash table cell
    std::unordered_map<uint64_t, unsigned int> kmerCount;

    // iterator to the element that we process
    auto srcItK = kmerTable.begin();
    auto srcItS = strandTable.begin();
    // iterator to the position of the last inserted FILTERED element, always at most as far as srcIt
    auto filterItK = kmerTable.begin();
    auto filterItS = strandTable.begin();

    // iterate over hashTable
    for (uint64_t i = 0; i < MyConst::HTABSIZE; ++i)
    {

        const uint64_t prevSizeK = filterItK - kmerTable.begin();

        // check if this slice contains enough elements - if not, do not make blacklist and copy all
        if ( (tabIndex[i+1] - tabIndex[i]) < MyConst::KMERCUTOFF)
        {

            // just copy over the old kmer slice
            for (uint64_t j = tabIndex[i]; j < tabIndex[i+1]; ++j, ++srcItK, ++srcItS, ++filterItK, ++filterItS)
            {

                *filterItK = *srcItK;
                *filterItS = *srcItS;

            }



        // if there are enough elements for potential kick outs, start a blacklisting
        } else {

            // clear container from previous round
            kmerCount.clear();

            // generate blacklist
            blacklist(tabIndex[i], tabIndex[i+1], kmerCount);

            // iterate through vector elements
            for (uint64_t j = tabIndex[i]; j < tabIndex[i+1]; ++j, ++srcItK, ++srcItS)
            {

                const uint64_t kHash = reproduceKmerSeq(*srcItK, *srcItS);

                // retrieve count how often it occurs
                if (kmerCount[kHash] < MyConst::KMERCUTOFF)
                {
                    // if it occurs not too often, copy
                    *filterItK = *srcItK;
                    *filterItS = *srcItS;
                    ++filterItK;
                    ++filterItS;

                } else {

                    // TODO
                    // filteredKmers.emplace(kHash);
                }
            }

        }
        // update the indexing structure for collisions
        tabIndex[i] = prevSizeK;
    }

    // update dummy value used for efficient indexing
    tabIndex[MyConst::HTABSIZE] = filterItK - kmerTable.begin();

    // shrink to new size
    kmerTable.resize(filterItK - kmerTable.begin());
    strandTable.resize(filterItS - strandTable.begin());

    std::chrono::high_resolution_clock::time_point filterEndTime = std::chrono::high_resolution_clock::now();

    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(filterEndTime - filterStartTime).count();


    std::cout << "Hash table size after filter (running " << runtime << "s): " << tabIndex[MyConst::HTABSIZE] << "\n\n";
}
