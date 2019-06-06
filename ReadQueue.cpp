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

#include <iostream>
#include <chrono>
#include <algorithm>

#include "ReadQueue.h"

ReadQueue::ReadQueue(const char* filePath, RefGenome& reference, const bool isGZ, const bool bsFlag) :
        ref(reference)
    ,   readBuffer(MyConst::CHUNKSIZE)
	,	isPaired(false)
	,	isSC(false)
    ,   methLevels(ref.cpgTable.size())
    ,   methLevelsStart(ref.cpgStartTable.size())
	,	bothStrandsFlag(bsFlag)
	,	r1FwdMatches(0)
	,	r1RevMatches(0)
	,	matchR1Fwd(true)
    //TODO
    ,   of("errOut.txt")
{
    isPaired = false;
    if (isGZ)
    {

        igz.open(filePath);

    } else {

        file.open(filePath);
    }

    // fill counting structure for parallelization
    for (unsigned int i = 0; i < CORENUM; ++i)
    {

        fwdMetaIDs[i] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        revMetaIDs[i] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        fwdMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 10);
        revMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 10);
        fwdMetaIDs[i].set_empty_key(ref.metaWindows.size() + 11);
        revMetaIDs[i].set_empty_key(ref.metaWindows.size() + 11);
    }
    // fill array mapping - locale specific filling
    lmap['A'%16] = 0;
    lmap['C'%16] = 1;
    lmap['G'%16] = 2;
    lmap['T'%16] = 3;

}
ReadQueue::ReadQueue(const char* filePath, const char* filePath2, RefGenome& reference, const bool isGZ, const bool bsFlag) :
        ref(reference)
    ,   readBuffer(MyConst::CHUNKSIZE)
    ,   readBuffer2(MyConst::CHUNKSIZE)
	,	isPaired(true)
	,	isSC(false)
    ,   methLevels(ref.cpgTable.size())
    ,   methLevelsStart(ref.cpgStartTable.size())
	,	bothStrandsFlag(bsFlag)
	,	r1FwdMatches(0)
	,	r1RevMatches(0)
	,	matchR1Fwd(true)
	// TODO
    ,   of("errOut.txt")
{

    isPaired = true;
    if (isGZ)
    {

        igz.open(filePath);
        igz2.open(filePath2);

    } else {

        file.open(filePath);
        file2.open(filePath2);
    }

    // fill counting structure for parallelization
    for (unsigned int i = 0; i < CORENUM; ++i)
    {

        // paired_fwdMetaIDs[i] = google::dense_hash_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>();
        // paired_revMetaIDs[i] = google::dense_hash_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>();
        // paired_fwdMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 10);
        // paired_revMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 10);
        // paired_fwdMetaIDs[i].set_empty_key(ref.metaCpGs.size() + 11);
        // paired_revMetaIDs[i].set_empty_key(ref.metaCpGs.size() + 11);
		paired_fwdMetaIDs[i] = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>();
		paired_revMetaIDs[i] = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>();
    }
    // fill array mapping - locale specific filling
    lmap['A'%16] = 0;
    lmap['C'%16] = 1;
    lmap['G'%16] = 2;
    lmap['T'%16] = 3;
}
ReadQueue::ReadQueue(const char* scOutputPath, RefGenome& reference, const bool isGZ, const bool bsFlag, const bool isP) :
        ref(reference)
    ,   readBuffer(MyConst::CHUNKSIZE)
    ,   readBuffer2(MyConst::CHUNKSIZE)
	,	isPaired(isP)
	,	isSC(true)
    ,   methLevels(ref.cpgTable.size())
    ,   methLevelsStart(ref.cpgStartTable.size())
	,	methLevelsSc(ref.cpgTable.size())
	,	scOutput(scOutputPath)
	,	bothStrandsFlag(bsFlag)
	,	r1FwdMatches(0)
	,	r1RevMatches(0)
	,	matchR1Fwd(true)
	// TODO
    ,   of("errOut.txt")
{

    // fill array mapping - locale specific filling
    lmap['A'%16] = 0;
    lmap['C'%16] = 1;
    lmap['G'%16] = 2;
    lmap['T'%16] = 3;
	// initialize sc output file
	scOutput << "\nSC_ID\tCount_Type\t";
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {
        scOutput << ref.chrMap[ref.cpgTable[cpgID].chrom] << "\t";
	}
	scOutput << "\nSC_ID\tCount_Type\t";
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {
        scOutput << ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 << "\t";
	}
}

bool ReadQueue::parseChunk(unsigned int& procReads)
{

    // std::cout << "Start reading chunk of reads\n";
    //
    std::string id;

    // counter on how many reads have been read so far
    unsigned int readCounter = 0;

    // read first line of read (aka @'SEQID')
    while (std::getline(file, id))
    {

        // read the next line (aka raw sequence)
        std::string seq;
        std::getline(file, seq);
        // construct read and push it to buffer
        readBuffer[readCounter] = Read(seq, id);
        // read the rest of read (aka +'SEQID' and quality score sequence)
        std::getline(file,id);
        std::getline(file,seq);

        ++readCounter;

        // if buffer is read completely, return
        if (readCounter >= MyConst::CHUNKSIZE)
        {
            procReads = MyConst::CHUNKSIZE;
            break;

        }
    }
    // if needed, read paired reads
    if (isPaired)
    {

        unsigned int readCounter2 = 0;
        // read first line of read (aka @'SEQID')
        while (std::getline(file2, id))
        {

            // read the next line (aka raw sequence)
            std::string seq;
            std::getline(file2, seq);
            // construct read and push it to buffer
            readBuffer2[readCounter2] = Read(seq, id);
            // read the rest of read (aka +'SEQID' and quality score sequence)
            std::getline(file2,id);
            std::getline(file2,seq);

            ++readCounter2;

            // if buffer is read completely, return
            if (readCounter2 >= MyConst::CHUNKSIZE)
            {
                procReads = MyConst::CHUNKSIZE;
                return true;

            }
        }
        // check if same number of reads is processed so far
        if (readCounter != readCounter2)
        {
            std::cerr << "Not the same number of reads available in the paired read files! \
                            Make sure that you paired all reads. \
                            Single reads have to be processed separately.\n\n";
            exit(1);
        }

    } else {

        if (readCounter >= MyConst::CHUNKSIZE)
            return true;
    }

    procReads = readCounter;
    return false;
}

bool ReadQueue::parseChunkGZ(unsigned int& procReads)
{

    std::string id;

    // counter on how many reads have been read so far
    unsigned int readCounter = 0;

    // read first line of read (aka @'SEQID')
    while (std::getline(igz, id))
    {

        // read the next line (aka raw sequence)
        std::string seq;
        std::getline(igz, seq);
        // construct read and push it to buffer
        readBuffer[readCounter] = Read(seq, id);
        // read the rest of read (aka +'SEQID' and quality score sequence)
        std::getline(igz,id);
        std::getline(igz,seq);

        ++readCounter;

        // if buffer is read completely, return
        if (readCounter >= MyConst::CHUNKSIZE)
        {
            procReads = MyConst::CHUNKSIZE;
            break;

        }
    }
    // if needed, read paired reads
    if (isPaired)
    {

        unsigned int readCounter2 = 0;
        // read first line of read (aka @'SEQID')
        while (std::getline(igz2, id))
        {

            // read the next line (aka raw sequence)
            std::string seq;
            std::getline(igz2, seq);
            // construct read and push it to buffer
            readBuffer2[readCounter2] = Read(seq, id);
            // read the rest of read (aka +'SEQID' and quality score sequence)
            std::getline(igz2,id);
            std::getline(igz2,seq);

            ++readCounter2;

            // if buffer is read completely, return
            if (readCounter2 >= MyConst::CHUNKSIZE)
            {
                procReads = MyConst::CHUNKSIZE;
                return true;

            }
        }
        // check if same number of reads is processed so far
        if (readCounter != readCounter2)
        {
            std::cerr << "Not the same number of reads available in the paired read files! \
                            Make sure that you paired all reads. \
                            Single reads have to be processed separately.\n\n";
            exit(1);
        }
    } else {

        if (readCounter >= MyConst::CHUNKSIZE)
            return true;
    }

    procReads = readCounter;
    return false;
}


void ReadQueue::decideStrand()
{
	std::cout << "\nFwd matches: " << r1FwdMatches;
	std::cout << "\nRev matches: " << r1RevMatches;
	const double odds = (double)(r1FwdMatches + 1)/(double)(r1RevMatches + 1);
	std::cout << "\nOdds of matching to fwd/rev strand: " << odds << "\n\n";
	if (odds > 0.005 && odds < 200)
	{
		std::cout << "Warning! Many of the reads are mapped in different orientation\
			Stranding might harm the prediction performance.\n\
			If you built a library where read 1 can occur as C->T or G->A converted, consider running the tool with\n\
			\"--non_stranded\"\n\
			flag.\n\n";
	}

	std::cout << "Deciding conversion status of read 1.\n";
	if (r1FwdMatches > r1RevMatches)
	{
		std::cout << "\tMatching read 1 as C->T converted.\n\n";
		matchR1Fwd = true;
	} else {
		std::cout << "\tMatching read 1 as G->A converted.\n\n";
		matchR1Fwd = false;
	}
}



bool ReadQueue::matchReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch, const bool getStranded)
{

    // reset all counters
    for (unsigned int i = 0; i < CORENUM; ++i)
    {
        matchStats[i] = 0;
        nonUniqueStats[i] = 0;
        noMatchStats[i] = 0;
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(CORENUM) schedule(static)
#endif
    for (unsigned int i = 0; i < procReads; ++i)
    {

        int threadnum = omp_get_thread_num();

        uint64_t& succMatchT = matchStats[threadnum];
        uint64_t& nonUniqueMatchT = nonUniqueStats[threadnum];
        uint64_t& unSuccMatchT = noMatchStats[threadnum];
        Read& r = readBuffer[i];

        const size_t readSize = r.seq.size();

        if (readSize < MyConst::READLEN - 10)
        {

            r.isInvalid = true;
            continue;
        }

        // flag stating if read contains N
        // reads with N are ignored
        bool nflag = false;

        // get correct offset for reverse strand (strand orientation must be correct)
        size_t revPos = readSize - 1;

        // std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
        // string containing reverse complement (under FULL alphabet)
        std::string revSeq;
        revSeq.resize(readSize);

        // construct reduced alphabet sequence for forward and reverse strand
        for (size_t pos = 0; pos < readSize; ++pos, --revPos)
        {

            switch (r.seq[pos])
            {
                case 'A':

                    revSeq[revPos] = 'T';
                    break;

                case 'C':

                    revSeq[revPos] = 'G';
                    break;

                case 'G':

                    revSeq[revPos] = 'C';
                    break;

                case 'T':

                    revSeq[revPos] = 'A';
                    break;

                case 'N':

                    nflag = true;
                    break;

                default:

                    std::cerr << "Unknown character '" << r.seq[pos] << "' in read with sequence id " << r.id << std::endl;
            }
        }

        if (nflag)
        {
            r.isInvalid = true;
            continue;
        }

        // std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
        // auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        //
        // of << runtime << "\n";
		//

        // set qgram threshold
		uint16_t qThreshold = MyConst::QTHRESH;
        // TODO
        // startTime = std::chrono::high_resolution_clock::now();
        MATCH::match matchFwd = 0;
        getSeedRefs(r.seq, readSize, qThreshold);
        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        // of << runtime << "\n";
        // if (runtime > 3000)
        // {
        //     of << r.seq << "\n";
        //     of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
        //     // uint16_t qThreshold = readSize - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT);
        //     // check for overflow (i.e. read is to small for lemma)
        //     // if (qThreshold > readSize)
        //     //     qThreshold = 0;
        //     uint64_t qcount = 0;
        //     for (const auto& m : fwdMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= qThreshold)
        //             ++qcount;
        //     }
        //     for (const auto& m : revMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= qThreshold)
        //             ++qcount;
        //     }
        //     of << "Meta CpGs passing q-gram (q=" << qThreshold << ") filter: " << qcount << "\n";
        // }
        // startTime = std::chrono::high_resolution_clock::now();
		int succQueryFwd = 0;
		if (bothStrandsFlag || getStranded || matchR1Fwd)
		{
			ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saFwd(r.seq, lmap);
			// endTime = std::chrono::high_resolution_clock::now();
			// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
			// of << runtime << "\t";
			// startTime = std::chrono::high_resolution_clock::now();
			// of << "--------------------------------\n\n";
			// of << "Matching read " << r.id << "\n\n";
			succQueryFwd = saQuerySeedSetRef(saFwd, matchFwd, qThreshold);
			// succQueryFwd = matchSingle(r.seq, qThreshold, saFwd, matchFwd, threadnum);
		}


        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        // of << runtime << "\n";
        //
        // startTime = std::chrono::high_resolution_clock::now();
        MATCH::match matchRev = 0;
        getSeedRefs(revSeq, readSize, qThreshold);
        // if (!succFlag)
        // {
        //     ++readCount;
        //     // of << "\n";
        //     r.isInvalid = true;
        //     continue;
        // }
        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        // of << runtime << "\n";
        // if (runtime > 3000)
        // {
        //     of << r.seq << "\n";
        //     of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
        //     // uint16_t qThreshold = readSize - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT);
        //     // check for overflow (i.e. read is to small for lemma)
        //     // if (qThreshold > readSize)
        //     //     qThreshold = 0;
        //     uint64_t qcount = 0;
        //     for (const auto& m : fwdMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= qThreshold)
        //             ++qcount;
        //     }
        //     for (const auto& m : revMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= qThreshold)
        //             ++qcount;
        //     }
        //     of << "Meta CpGs passing q-gram (q=" << qThreshold << ") filter: " << qcount << "\n";
        // }
        // of << "\nMeta id count fwd: "  << fwdMetaIDs[omp_get_thread_num()][319746] << "\n";
        // of << "\nMeta id count rev: "  << revMetaIDs[omp_get_thread_num()][319746] << "\n";
        // if (runtime > 100)
        // {
        //     of << r.seq << "\n";
        //     of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
        // }
        // startTime = std::chrono::high_resolution_clock::now();
        // std::cout << revSeq << "\n";
		int succQueryRev = 0;
		if (bothStrandsFlag || getStranded || !matchR1Fwd)
		{
			ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saRev(revSeq, lmap);
			// endTime = std::chrono::high_resolution_clock::now();
			// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
			// of << runtime << "\n";
			// startTime = std::chrono::high_resolution_clock::now();
			succQueryRev = saQuerySeedSetRef(saRev, matchRev, qThreshold);
			// succQueryRev = matchSingle(revSeq, qThreshold, saRev, matchRev, threadnum);
		}


        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        // of << runtime << "\n";
        // if (runtime > 5000)
        // {
            // of << r.seq << "\n";
            // of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
            // // uint16_t qThreshold = readSize - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT);
            // // check for overflow (i.e. read is to small for lemma)
            // // if (qThreshold > readSize)
            // //     qThreshold = 0;
            // qcount = 0;
            // for (const auto& m : fwdMetaIDs[omp_get_thread_num()])
            // {
            //     if (m.second >= qThreshold)
            //         ++qcount;
            // }
            // for (const auto& m : revMetaIDs[omp_get_thread_num()])
            // {
            //     if (m.second >= qThreshold)
            //         ++qcount;
            // }
            // of << "Meta CpGs passing q-gram (q=" << qThreshold << ") filter: " << qcount << "\n";
        // }

        // found match for fwd and rev automaton
        if (succQueryFwd == 1 && succQueryRev == 1)
        {

            uint8_t fwdErr = MATCH::getErrNum(matchFwd);
            uint8_t revErr = MATCH::getErrNum(matchRev);
            // of << "Found match with read and reverse complement. Errors (fwd/rev): " << fwdErr << "/" << revErr;
            // of << "\nMatching strands are (fwd/rev): " << MATCH::isFwd(matchFwd) << "/" << MATCH::isFwd(matchRev) << "\n";

            // check which one has fewer errors
            if (fwdErr < revErr)
            {

				if (getStranded)
#pragma omp atomic
					++r1FwdMatches;
                ++succMatchT;
                r.mat = matchFwd;
                // startTime = std::chrono::high_resolution_clock::now();
                computeMethLvl(matchFwd, r.seq);
                // endTime = std::chrono::high_resolution_clock::now();
                // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                // of << runtime << "\n";

            } else {

                if (fwdErr > revErr)
                {

					if (getStranded)
#pragma omp atomic
						++r1RevMatches;
                    ++succMatchT;
                    r.mat = matchRev;
                    // startTime = std::chrono::high_resolution_clock::now();
                    computeMethLvl(matchRev, revSeq);
                    // endTime = std::chrono::high_resolution_clock::now();
                    // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    // of << runtime << "\n";

                // if same number of errors, then not unique
                } else {

                    const uint32_t metaFwd = MATCH::getMetaID(matchFwd);
                    const uint32_t metaRev = MATCH::getMetaID(matchRev);
                    const uint64_t offFwd = MATCH::getOffset(matchFwd);
                    const uint64_t offRev = MATCH::getOffset(matchRev);
                    const bool m1_isFwd = MATCH::isFwd(matchFwd);
                    const bool m2_isFwd = MATCH::isFwd(matchRev);
                    const bool m1_isStart = MATCH::isStart(matchFwd);
                    const bool m2_isStart = MATCH::isStart(matchRev);
                    uint32_t m1_pos;
                    uint32_t m2_pos;
                    if (m1_isStart && m2_isStart)
                    {

                        m1_pos = offFwd;
                        m2_pos = offRev;

                    } else {

                        m1_pos = ref.metaWindows[metaFwd].startPos + offFwd;
                        m2_pos = ref.metaWindows[metaRev].startPos + offRev;
                    }
                    // test if same match in same region
                    if ((m1_isStart == m2_isStart) && (m1_isFwd == m2_isFwd) && (m1_pos == m2_pos))
                    {
                        ++succMatchT;
						if (getStranded)
#pragma omp atomic
							++r1FwdMatches;
                        r.mat = matchFwd;
                        // startTime = std::chrono::high_resolution_clock::now();
                        computeMethLvl(matchFwd, r.seq);
                        // endTime = std::chrono::high_resolution_clock::now();
                        // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                        // of << runtime << "\n";

                    } else {

                        ++nonUniqueMatchT;

                        r.isInvalid = true;
                    }
                }
            }
        // unique match on forward strand
        } else if (succQueryFwd == 1) {

            if (succQueryRev == -1)
            {
                if (MATCH::getErrNum(matchFwd) < MATCH::getErrNum(matchRev))
                {
                    // of << "Match with FWD automaton. Strand is " << MATCH::isFwd(matchFwd) << "\n";
                    ++succMatchT;
					if (getStranded)
#pragma omp atomic
						++r1FwdMatches;
                    r.mat = matchFwd;
                    // startTime = std::chrono::high_resolution_clock::now();
                    computeMethLvl(matchFwd, r.seq);
                    // endTime = std::chrono::high_resolution_clock::now();
                    // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    // of << runtime << "\n";
                } else {

                    ++nonUniqueMatchT;
                    r.isInvalid = true;
                }
            } else {

                // of << "Match with FWD automaton. Strand is " << MATCH::isFwd(matchFwd) << "\n";
                ++succMatchT;
				if (getStranded)
#pragma omp atomic
					++r1FwdMatches;
                r.mat = matchFwd;
                // startTime = std::chrono::high_resolution_clock::now();
                computeMethLvl(matchFwd, r.seq);
                // endTime = std::chrono::high_resolution_clock::now();
                // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                // of << runtime << "\n";
            }

        // unique match on backward strand
        } else if (succQueryRev == 1) {

            if (succQueryFwd == -1)
            {
                if (MATCH::getErrNum(matchRev) < MATCH::getErrNum(matchFwd))
                {
                    // of << "Match with REV automaton. Strand is " << MATCH::isFwd(matchRev) << "\n";
                    ++succMatchT;
					if (getStranded)
#pragma omp atomic
						++r1RevMatches;
                    r.mat = matchRev;
                    // startTime = std::chrono::high_resolution_clock::now();
                    computeMethLvl(matchRev, revSeq);
                    // endTime = std::chrono::high_resolution_clock::now();
                    // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    // of << runtime << "\n";
                } else {

                    ++nonUniqueMatchT;
                    r.isInvalid = true;
                }
            } else {

                // of << "Match with REV automaton. Strand is " << MATCH::isFwd(matchRev) << "\n";
                ++succMatchT;
				if (getStranded)
#pragma omp atomic
					++r1RevMatches;
                r.mat = matchRev;
                // startTime = std::chrono::high_resolution_clock::now();
                computeMethLvl(matchRev, revSeq);
                // endTime = std::chrono::high_resolution_clock::now();
                // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                // of << runtime << "\n";
            }

        // no match found at all
        } else {

            r.isInvalid = true;
            if (succQueryFwd == -1 || succQueryRev == -1)
            {

                // of << "Nonunique match.\n";
                ++nonUniqueMatchT;

            } else {

                // of << "No match.\n";
                ++unSuccMatchT;

        //     }
        // }

// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//                 // construct hash and look up the hash table entries
//                 size_t lPos = r.id.find_last_of('_');
//                 std::string stringOffset (r.id.begin() + lPos + 1, r.id.end());
//                 size_t rPos = r.id.find_last_of('R');
//                 std::string stringChrom (r.id.begin() + 1 + rPos, r.id.begin() + lPos);
//                 uint8_t chrom = std::stoul(stringChrom);
//                 unsigned long offset = std::stoul(stringOffset);
//                 of << "\nreal seq/real revSeq/sequence in genome: " << r.id << "\n" << r.seq << "\n" << revSeq << "\n" << std::string(ref.fullSeq[chrom].begin() + offset, ref.fullSeq[chrom].begin() + offset + 100) << "\n\n\n";
//
//
//                 uint64_t hVal;
// 				uint64_t sVal = ntHash::NTPS64(r.seq.data(), MyConst::SEED, MyConst::KMERLEN, hVal) % MyConst::HTABSIZE;
//                 auto startIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal];
//                 auto endIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal + 1];
//                 auto tit = ref.strandTable.begin() + ref.tabIndex[sVal];
//                 for (auto it = startIt; it != endIt; ++it, ++tit)
//                 {
//                     KMER_S::kmer& k = *it;
//                     const uint32_t m = KMER_S::getMetaCpG(k);
//                     const bool isStart = KMER_S::isStartCpG(k);
//                     if (!isStart)
//                     {
//                         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
//                         if (*tit)
//                         {
//                             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
//                             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
//                             of << std::string(stIt, enIt) << "\n";
//                         }
//                     }
//                 }
//                 of << "Last Sequence part:\n";
//                 hVal = 0;
// 				sVal = ntHash::NTPS64(r.seq.data()+70, MyConst::SEED, MyConst::KMERLEN, hVal) % MyConst::HTABSIZE;
//                 startIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal];
//                 endIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal + 1];
//                 tit = ref.strandTable.begin() + ref.tabIndex[sVal];
//                 for (auto it = startIt; it != endIt; ++it, ++tit)
//                 {
//                     KMER_S::kmer& k = *it;
//                     const uint32_t m = KMER_S::getMetaCpG(k);
//                     const bool isStart = KMER_S::isStartCpG(k);
//                     if (!isStart)
//                     {
//                         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
//                         if (*tit)
//                         {
//                             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
//                             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
//                             of << std::string(stIt, enIt) << "\n";
//                         }
//                     }
//                 }
//
//                 of << "\n\nReverse seq matches\n";
//                 hVal = 0;
// 				sVal = ntHash::NTPS64(revSeq.data(), MyConst::SEED, MyConst::KMERLEN, hVal) % MyConst::HTABSIZE;
//                 startIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal];
//                 endIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal + 1];
//                 tit = ref.strandTable.begin() + ref.tabIndex[sVal];
//                 for (auto it = startIt; it != endIt; ++it, ++tit)
//                 {
//                     KMER_S::kmer& k = *it;
//                     const uint64_t m = KMER_S::getMetaCpG(k);
//                     const bool isStart = KMER_S::isStartCpG(k);
//                     if (!isStart)
//                     {
//                         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
//                         if (*tit)
//                         {
//                             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
//                             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
//                             of << std::string(stIt, enIt) << "\n";
//                         }
//                     }
//                 }
//                 of << "Last Sequence part:\n";
//                 hVal = 0;
// 				sVal = ntHash::NTPS64(revSeq.data()+70, MyConst::SEED, MyConst::KMERLEN, hVal) % MyConst::HTABSIZE;
//                 startIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal];
//                 endIt = ref.kmerTableSmall.begin() + ref.tabIndex[sVal + 1];
//                 tit = ref.strandTable.begin() + ref.tabIndex[sVal];
//                 for (auto it = startIt; it != endIt; ++it, ++tit)
//                 {
//                     KMER_S::kmer& k = *it;
//                     const uint64_t m = KMER_S::getMetaCpG(k);
//                     const bool isStart = KMER_S::isStartCpG(k);
//                     if (!isStart)
//                     {
//                         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
//                         if (*tit)
//                         {
//                             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
//                             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
//                             of << std::string(stIt, enIt) << "\n";
//                         }
//                     }
//                 }
//                 of << "\n\n--------------------\n\n";
//
//
// // END PRAGMA OMP CRITICAL
// }
            }
            // if (unSuccMatch > 10)
            // {
            //     of.close();
            //     exit(1);
            // }
        }
    }

    // sum up counts
    for (unsigned int i = 0; i < CORENUM; ++i)
    {
        succMatch += matchStats[i];
        nonUniqueMatch += nonUniqueStats[i];
        unSuccMatch += noMatchStats[i];
    }
    return true;
}

bool ReadQueue::matchPairedReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch, uint64_t& succPairedMatch, uint64_t& tooShortCountMatch, const bool getStranded)
{


    // reset all counters
    for (unsigned int i = 0; i < CORENUM; ++i)
    {
        matchStats[i] = 0;
        nonUniqueStats[i] = 0;
        noMatchStats[i] = 0;
        matchPairedStats[i] = 0;
		tooShortCounts[i] = 0;
    }

	// TODO
	// std::ofstream of2 ("errOut2.txt");

#ifdef _OPENMP
#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,50)
#endif
    for (unsigned int i = 0; i < procReads; ++i)
    {

        int threadnum = omp_get_thread_num();

        uint64_t& succMatchT = matchStats[threadnum];
        uint64_t& nonUniqueMatchT = nonUniqueStats[threadnum];
        uint64_t& unSuccMatchT = noMatchStats[threadnum];
        uint64_t& succPairedMatchT = matchPairedStats[threadnum];
		uint64_t& tooShortCount = tooShortCounts[threadnum];
        Read& r1 = readBuffer[i];
        Read& r2 = readBuffer2[i];

        const size_t readSize1 = r1.seq.size();
        const size_t readSize2 = r2.seq.size();


        if (readSize1 < MyConst::READLEN - 20)
        {

            r1.isInvalid = true;
        }
        if (readSize2 < MyConst::READLEN - 20)
        {

            r2.isInvalid = true;
        }
		if (r1.isInvalid || r2.isInvalid)
		{
			++tooShortCount;
		}

        // get correct offset for reverse strand (strand orientation must be correct)
        size_t revPos = readSize1 - 1;

        // string containing reverse complement (under FULL alphabet)
        std::string revSeq1;
        revSeq1.resize(readSize1);

        // construct reduced alphabet sequence for forward and reverse strand
        for (size_t pos = 0; pos < readSize1; ++pos, --revPos)
        {

            switch (r1.seq[pos])
            {
                case 'A':

                    revSeq1[revPos] = 'T';
                    break;

                case 'C':

                    revSeq1[revPos] = 'G';
                    break;

                case 'G':

                    revSeq1[revPos] = 'C';
                    break;

                case 'T':

                    revSeq1[revPos] = 'A';
                    break;

                case 'N':

                    r1.isInvalid = true;
                    break;

                default:

                    std::cerr << "Unknown character '" << r1.seq[pos] << "' in read with sequence id " << r1.id << std::endl;
                    r1.isInvalid = true;
            }
        }

        // get correct offset for reverse strand (strand orientation must be correct)
        revPos = readSize2 - 1;

        // string containing reverse complement (under FULL alphabet)
        std::string revSeq2;
        revSeq2.resize(readSize2);

        // construct reduced alphabet sequence for forward and reverse strand
        for (size_t pos = 0; pos < readSize2; ++pos, --revPos)
        {

            switch (r2.seq[pos])
            {
                case 'A':

                    revSeq2[revPos] = 'T';
                    break;

                case 'C':

                    revSeq2[revPos] = 'G';
                    break;

                case 'G':

                    revSeq2[revPos] = 'C';
                    break;

                case 'T':

                    revSeq2[revPos] = 'A';
                    break;

                case 'N':

                    r2.isInvalid = true;
                    break;

                default:

                    std::cerr << "Unknown character '" << r2.seq[pos] << "' in read with sequence id " << r2.id << std::endl;
                    r2.isInvalid = true;
            }
        }

        if (r1.isInvalid || r2.isInvalid)
        {
            tooShortCount += 2;
            continue;
        }

        // set qgram threshold
        // uint16_t qThreshold = std::min(readSize1, readSize2) - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT) + 1;
        // if (qThreshold > readSize1)
        //     qThreshold = 0;
        //
		// TODO
		const uint16_t qThreshold = MyConst::QTHRESH;
		// of << "\nq-gram: " << qThreshold << "\n";

// #pragma omp critical
// {
//
//         // construct hash and look up the hash table entries
//         size_t lPos = r1.id.find_last_of('_');
//         std::string stringOffset (r1.id.begin() + lPos + 1, r1.id.end());
//         size_t rPos = r1.id.find_last_of('r');
//         std::string stringChrom (r1.id.begin() + 1 + rPos, r1.id.begin() + lPos);
//         uint8_t chrom = std::stoul(stringChrom);
//         unsigned long offset = std::stoul(stringOffset);
//         of << "\nRead1 real seq/real revSeq/sequence in genome: " << r1.id << "\n" << r1.seq << "\n" << revSeq1 << "\n" << std::string(ref.fullSeq[chrom].begin() + offset, ref.fullSeq[chrom].begin() + offset + 100) << "\n\n\n";
//         lPos = r2.id.find_last_of('_');
//         stringOffset = std::string(r2.id.begin() + lPos + 1, r2.id.end());
//         rPos = r2.id.find_last_of('r');
//         stringChrom = std::string(r2.id.begin() + 1 + rPos, r2.id.begin() + lPos);
//         chrom = std::stoul(stringChrom);
//         offset = std::stoul(stringOffset);
//         of << "\nRead2 real seq/real revSeq/sequence in genome: " << r2.id << "\n" << r2.seq << "\n" << revSeq2 << "\n" << std::string(ref.fullSeq[chrom].begin() + offset, ref.fullSeq[chrom].begin() + offset + 100) << "\n\n\n";
//
// }

    // POSSIBLE ORIENTATION 1
		std::vector<MATCH::match> matches1Fwd;
		std::vector<MATCH::match> matches2Rev;

		// std::chrono::high_resolution_clock::time_point startTime;
		// std::chrono::high_resolution_clock::time_point endTime;

		if (bothStrandsFlag || getStranded || matchR1Fwd)
		{

			// startTime = std::chrono::high_resolution_clock::now();

			const uint16_t qAdapt1 = getSeedRefsFirstRead(r1.seq, readSize1, qThreshold);
// #pragma omp critical
// 				{
// 					of << "Found metas after first seeding: " << paired_fwdMetaIDs[threadnum].size() << "\t" << paired_revMetaIDs[threadnum].size() << "\n";
// 					of << "Meta IDs fwd:\t";
// 					for (const auto& m : paired_fwdMetaIDs[threadnum])
// 					{
// 						of << m.first << "/" << ref.metaWindows[m.first].startPos << "\t";
// 						of << (uint32_t)std::get<0>(m.second) << "\t" << (uint32_t)std::get<1>(m.second) << "\t" << std::get<2>(m.second) << "\t" << std::get<3>(m.second) << "\n";
// 					}
// 					of << "\n";
// 					of << "Meta IDs rev:\t";
// 					for (const auto& m : paired_revMetaIDs[threadnum])
// 					{
// 						of << m.first  << "/" << ref.metaWindows[m.first].startPos << "\t";
// 						of << (uint32_t)std::get<0>(m.second) << "\t" << (uint32_t)std::get<1>(m.second) << "\t" << std::get<2>(m.second) << "\t" << std::get<3>(m.second) << "\n";
// 					}
// 					of << "\n";
// 				}
			// if (qAdapt1)
			// {

			// endTime = std::chrono::high_resolution_clock::now();
			// auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
			// of << "\nHashtable sizes after first read: " << paired_fwdMetaIDs[threadnum].size() << "/"<< paired_fwdMetaIDs[threadnum].size()<< "\t\tRuntime: "  << runtime << "\t\tCandidates1: " << candCount << "\n";

			if ((paired_fwdMetaIDs[threadnum].size() || paired_revMetaIDs[threadnum].size()) )
			{
				// startTime = std::chrono::high_resolution_clock::now();

				const bool hasCpG = getSeedRefsSecondRead(revSeq2, readSize2, qThreshold);
// #pragma omp critical
// 				{
// 					of << "Found metas after second seeding: " << paired_fwdMetaIDs[threadnum].size() << "\t" << paired_revMetaIDs[threadnum].size() << "\n";
// 					of << "Meta IDs fwd:\t";
// 					uint32_t wCount = 0;
// 					uint32_t wMax = 0;
// 					for (const auto& m : paired_fwdMetaIDs[threadnum])
// 					{
// 						if (std::get<0>(m.second) > wMax)
// 						{
// 							wMax = std::get<0>(m.second);
//
// 						}
// 						if (std::get<2>(m.second))
// 							++wCount;
//
// 						// of << m.first << "/" << ref.metaWindows[m.first].startPos << "\t";
// 						// of << (uint32_t)std::get<0>(m.second) << "\t" << (uint32_t)std::get<1>(m.second) << "\t" << std::get<2>(m.second) << "\t" << std::get<3>(m.second) << "\n";
// 					}
// 					// of << "\n";
// // 					of << "Meta IDs rev:\t";
// 					for (const auto& m : paired_revMetaIDs[threadnum])
// 					{
// 						if (std::get<0>(m.second) > wMax)
// 						{
// 							wMax = std::get<0>(m.second);
//
// 						}
// 						if (std::get<2>(m.second))
// 							++wCount;
//
// 						// of << m.first << "/" << ref.metaWindows[m.first].startPos << "\t";
// 						// of << (uint32_t)std::get<0>(m.second) << "\t" << (uint32_t)std::get<1>(m.second) << "\t" << std::get<2>(m.second) << "\t" << std::get<3>(m.second) << "\n";
// 					}
// 					for (const auto& m : paired_fwdMetaIDs[threadnum])
// 					{
// 						if (std::get<2>(m.second))
// 						{
// 							of2 << wMax - std::get<0>(m.second) << "\n";
// 						}
// 					}
// 					for (const auto& m : paired_revMetaIDs[threadnum])
// 					{
// 						if (std::get<2>(m.second))
// 						{
// 							of2 << wMax - std::get<0>(m.second) << "\n";
// 						}
// 					}
// 					// of << "\n";
// 					of << wCount << "\n";
// 				}
				// if (qAdapt2)
				// {
				if (hasCpG)
				{
					ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saFwd(r1.seq, lmap);
					ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saRev2(revSeq2, lmap);

					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << "\nHashtable sizes after second read: " << paired_fwdMetaIDs[threadnum].size() << "/"<< paired_fwdMetaIDs[threadnum].size()<< "\t\tRuntime: "  << runtime << "\t\tCandidates2: " << candCount << "\n";
					// startTime = std::chrono::high_resolution_clock::now();


					saQuerySeedSetRefFirst(saFwd, matches1Fwd, qThreshold);
					if (!matches1Fwd.empty())
						saQuerySeedSetRefSecond(saRev2, matches2Rev, qThreshold);
					// }

					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << runtime << "\n\n";
				}
				// }
			}
		}

    // POSSIBLE ORIENTATION 2
		std::vector<MATCH::match> matches1Rev;
		std::vector<MATCH::match> matches2Fwd;

		if (bothStrandsFlag || getStranded || !matchR1Fwd)
		{
			// startTime = std::chrono::high_resolution_clock::now();

			const uint16_t qAdapt1 = getSeedRefsFirstRead(revSeq1, readSize1, qThreshold);
			// if (qAdapt1)
			// {

			// endTime = std::chrono::high_resolution_clock::now();
			// auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
			// of << "\nHashtable sizes after first read: " << paired_fwdMetaIDs[threadnum].size() << "/"<< paired_fwdMetaIDs[threadnum].size()<< "\tRuntime: "  << runtime << "\n";

			if ((paired_fwdMetaIDs[threadnum].size() || paired_revMetaIDs[threadnum].size()))
			{
				// startTime = std::chrono::high_resolution_clock::now();


				const bool hasCpG = getSeedRefsSecondRead(r2.seq, readSize2, qThreshold);
				// if (qAdapt2)
				// {
				if (hasCpG)
				{
					ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saRev(revSeq1, lmap);
					ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saFwd2(r2.seq, lmap);

					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << "\nHashtable sizes after second read: " << paired_fwdMetaIDs[threadnum].size() << "/"<< paired_fwdMetaIDs[threadnum].size()<< "\tRuntime: "  << runtime << "\n";
					// startTime = std::chrono::high_resolution_clock::now();


					saQuerySeedSetRefFirst(saRev, matches1Rev, qThreshold);
					if (!matches1Rev.empty())
						saQuerySeedSetRefSecond(saFwd2, matches2Fwd, qThreshold);

					// }
					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << runtime << "\n";
				}
			}
			// }
		}

    // TEST IF MATCHING WAS SUCCESSFULL

        if ((matches1Fwd.size() == 0 && matches1Rev.size() == 0) || (matches2Fwd.size() == 0 && matches2Rev.size() == 0))
        {
// #pragma omp critical
// 			{
// 			of << "No Meta candidates found:\n\t" << r1.id << "\n\t" << r1.seq << "\n\t" << r2.id << "\n\t" << r2.seq << "\n\n";
// 			of << "-----------------------------\n\n\n\n";
// 			}
            r1.isInvalid = true;
            r2.isInvalid = true;
            unSuccMatchT += 2;
            continue;
        }

    // TRY TO PAIR MATCHES

        // of << "Matching sizes: " << matches1Fwd.size() << "/" << matches2Rev.size() << "\tTimings: ";
        // startTime = std::chrono::high_resolution_clock::now();
        // current best matching pair (sum of errors)
        int bestErrNum = 2*(MyConst::MISCOUNT + MyConst::ADDMIS) + 1;
        MATCH::match bestMatch1;
        MATCH::match bestMatch2;
        bool nonUniqueFlag = false;
        bool mat1OriginalStrand = true;
		// std::cout << matches1Fwd.size() << "\t" << matches1Rev.size() << "\t" << matches2Fwd.size() << "\t" << matches2Rev.size() << "\n";
		// printMatch(std::cout, matches1Fwd[0]);
		// printMatch(std::cout, matches2Rev[0]);

		if (bothStrandsFlag || getStranded || matchR1Fwd)
		{
			// uint32_t i = 0;
			for (auto m1It = matches1Fwd.begin(); m1It != matches1Fwd.end(); ++m1It)
			{
				// std::cout << i << "\n";
				MATCH::match& mat1 = *m1It;
				// ++i;
				for (auto m2It = matches2Rev.begin(); m2It != matches2Rev.end(); ++m2It)
				{
					// std::cout << i << "\n";
					MATCH::match& mat2Rev = *m2It;
					// ++i;
					int extractedMatchErrs = extractPairedMatch(mat1, mat2Rev);
					if (extractedMatchErrs >= 0)
					{
						if (extractedMatchErrs == bestErrNum)
						{
							nonUniqueFlag = true;

						} else if (extractedMatchErrs < bestErrNum) {

							bestErrNum = extractedMatchErrs;
							bestMatch1 = mat1;
							bestMatch2 = mat2Rev;
							nonUniqueFlag = false;

						}
					}
				}
			}
		}
        // endTime = std::chrono::high_resolution_clock::now();
        // auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        // of << runtime << "\n";
        // of << "Matching sizes: " << matches1Rev.size() << "/" << matches2Fwd.size() << "\tTimings: ";
        // startTime = std::chrono::high_resolution_clock::now();
		if (bothStrandsFlag || getStranded || !matchR1Fwd)
		{
			for (auto m1It = matches1Rev.begin(); m1It != matches1Rev.end(); ++m1It)
			{
				MATCH::match& mat1 = *m1It;
				for (auto m2It = matches2Fwd.begin(); m2It != matches2Fwd.end(); ++m2It)
				{
					MATCH::match& mat2Fwd = *m2It;
					int extractedMatchErrs = extractPairedMatch(mat1, mat2Fwd);
					if (extractedMatchErrs >= 0)
					{
						if (extractedMatchErrs == bestErrNum)
						{
							nonUniqueFlag = true;

						} else if (extractedMatchErrs < bestErrNum) {

							bestErrNum = extractedMatchErrs;
							bestMatch1 = mat1;
							bestMatch2 = mat2Fwd;
							mat1OriginalStrand = false;
							nonUniqueFlag = false;

						}
					}
				}
			}
		}
        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        // of << runtime << "\n\n";

// #pragma omp critical
// {
//
//         // construct hash and look up the hash table entries
//         size_t lPos = r1.id.find_last_of('_');
//         std::string stringOffset (r1.id.begin() + lPos + 1, r1.id.end());
//         size_t rPos = r1.id.find_last_of('r');
//         std::string stringChrom (r1.id.begin() + 1 + rPos, r1.id.begin() + lPos);
//         uint8_t chrom = std::stoul(stringChrom);
//         unsigned long offset = std::stoul(stringOffset);
//         of << "\nreal seq/real revSeq/sequence in genome: " << r1.id << "\n" << r1.seq << "\n" << revSeq1 << "\n" << std::string(ref.fullSeq[chrom].begin() + offset, ref.fullSeq[chrom].begin() + offset + 100) << "\n\n\n";
//
//
// 		uint64_t fhVal;
//         uint64_t hVal = ntHash::NTPS64(r1.seq.data(), MyConst::SEED, MyConst::KMERLEN, fhVal) % MyConst::HTABSIZE;
//         auto startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//         auto endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//         auto tit = ref.strandTable.begin() + ref.tabIndex[hVal];
//         for (auto it = startIt; it != endIt; ++it, ++tit)
//         {
//             KMER_S::kmer& k = *it;
//             const uint32_t m = KMER_S::getMetaCpG(k);
//             // const bool isStart = KMER_S::isStartCpG(k);
//             // if (!isStart)
//             // {
//                 const auto& w = ref.metaWindows[m];
//                 if (*tit)
//                 {
//                     auto stIt = ref.fullSeq[w.chrom].begin() + w.startPos;
//                     auto enIt = ref.fullSeq[w.chrom].begin() + MyConst::WINLEN + w.startPos;
//                     of << std::string(stIt, enIt) << "\n";
//                 }
//             // }
//         }
//         of << "Last Sequence part:\n";
//         hVal = ntHash::NTPS64(r1.seq.data()+65, MyConst::SEED, MyConst::KMERLEN, fhVal) % MyConst::HTABSIZE;
//         startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//         endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//         tit = ref.strandTable.begin() + ref.tabIndex[hVal];
//         for (auto it = startIt; it != endIt; ++it, ++tit)
//         {
//             KMER_S::kmer& k = *it;
//             const uint32_t m = KMER_S::getMetaCpG(k);
//             // const bool isStart = KMER_S::isStartCpG(k);
//             // if (!isStart)
//             // {
//                 const auto& w = ref.metaWindows[m];
//                 if (*tit)
//                 {
//                     auto stIt = ref.fullSeq[w.chrom].begin() + w.startPos;
//                     auto enIt = ref.fullSeq[w.chrom].begin() + MyConst::WINLEN + w.startPos;
//                     of << std::string(stIt, enIt) << "\n";
//                 }
//             // }
//         }
//
//         of << "\n\nReverse seq matches\n";
//         hVal = ntHash::NTPS64(revSeq1.data(), MyConst::SEED, MyConst::KMERLEN, fhVal) % MyConst::HTABSIZE;
//         startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//         endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//         tit = ref.strandTable.begin() + ref.tabIndex[hVal];
//         for (auto it = startIt; it != endIt; ++it, ++tit)
//         {
//             KMER_S::kmer& k = *it;
//             const uint64_t m = KMER_S::getMetaCpG(k);
//             // const bool isStart = KMER_S::isStartCpG(k);
//             // if (!isStart)
//             // {
//                 const auto& w = ref.metaWindows[m];
//                 if (*tit)
//                 {
//                     auto stIt = ref.fullSeq[w.chrom].begin() + w.startPos;
//                     auto enIt = ref.fullSeq[w.chrom].begin() + MyConst::WINLEN + w.startPos;
//                     of << std::string(stIt, enIt) << "\n";
//                 }
//             // }
//         }
//         of << "Last Sequence part:\n";
//         hVal = ntHash::NTPS64(revSeq1.data()+65, MyConst::SEED, MyConst::KMERLEN, fhVal) % MyConst::HTABSIZE;
//         startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//         endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//         tit = ref.strandTable.begin() + ref.tabIndex[hVal];
//         for (auto it = startIt; it != endIt; ++it, ++tit)
//         {
//             KMER_S::kmer& k = *it;
//             const uint64_t m = KMER_S::getMetaCpG(k);
//             // const bool isStart = KMER_S::isStartCpG(k);
//             // if (!isStart)
//             // {
//                 const auto& w = ref.metaWindows[m];
//                 if (*tit)
//                 {
//                     auto stIt = ref.fullSeq[w.chrom].begin() + w.startPos;
//                     auto enIt = ref.fullSeq[w.chrom].begin() + MyConst::WINLEN + w.startPos;
//                     of << std::string(stIt, enIt) << "\n";
//                 }
//             // }
//         }
// }


        // Check if no pairing possible
        if (bestErrNum == 2*(MyConst::MISCOUNT + MyConst::ADDMIS) + 1)
        {
// #pragma omp critical
// {
//             of << "\n\n\nNo pairing possible\n\n";
//             of << "Matches of read 1, ID " << r1.id << " \nfwd: " << r1.seq << "\n";
//             for (auto mat : matches1Fwd)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "Matches of read 1,  ID " << r1.id << " \nrev: " << revSeq1 << "\n";
//             for (auto mat : matches1Rev)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "\nMatches of read 2,  ID " << r2.id << " \nfwd: " << r2.seq << "\n";
//             for (auto mat : matches2Fwd)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "Matches of read 2,  ID " << r2.id << " \nrev: " << revSeq2 << "\n";
//             for (auto mat : matches2Rev)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "\n";
// // end pragma omp critical
// }
            //
            // if (extractSingleMatch(matches1Fwd, matches1Rev, r1, revSeq1))
            // {
            //     ++succMatchT;
// #pragma omp critical
// {
//                 of << "\tSuccessfull r1\n";
// }
            //
            // } else {
            //
            //     matches1Fwd.size() + matches1Rev.size() > 0 ? ++nonUniqueMatchT : ++unSuccMatchT;
// #pragma omp critical
// {
//                 of << "\tUnsuccessfull r1\n";
// }
            // }
            // if (extractSingleMatch(matches2Fwd, matches2Rev, r2, revSeq2))
            // {
            //
            //     ++succMatchT;
// #pragma omp critical
// {
//                 of << "\tSuccessfull r2\n";
// }
            // } else {
            //
            //     matches2Fwd.size() + matches2Rev.size() > 0 ? ++nonUniqueMatchT : ++unSuccMatchT;
// #pragma omp critical
// {
//                 of << "\tUnsuccessfull r1\n";
// }
            // }
			// of << "\nNo Match\n";
            r1.isInvalid = true;
            r2.isInvalid = true;
            unSuccMatchT += 2;

        } else if (nonUniqueFlag)
        {
// #pragma omp critical
// {
//             of << "\n\n\nNonunique pair\n\n";
//             of << "Matches of read 1,  ID " << r1.id << " \nfwd: " << r1.seq << "\n";
//             for (auto mat : matches1Fwd)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "Matches of read 1,  ID " << r1.id << " \nrev: " << revSeq1 << "\n";
//             for (auto mat : matches1Rev)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "\nMatches of read 2,  ID " << r2.id << " \nfwd: " << r2.seq << "\n";
//             for (auto mat : matches2Fwd)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "Matches of read 2,  ID " << r2.id << " \nrev: " << revSeq2 << "\n";
//             for (auto mat : matches2Rev)
//             {
//                 printMatch(of, mat);
//                 of << "\n";
//             }
//             of << "\n";
// // end pragma omp critical
// }
			// of << "\nWas nonunique\n";
            nonUniqueMatchT += 2;
            r1.isInvalid = true;
            r2.isInvalid = true;

        } else {

// #pragma omp critical
// {
// 			of << "\nWas match\n";
// }
            r1.mat = bestMatch1;
            r2.mat = bestMatch2;
			// TODO
// #pragma omp critical
// 			{
// 				of << "Paired Match\n";
// 				of << r1.id << "\n\t" << static_cast<uint64_t>(ref.metaWindows[MATCH::getMetaID(r1.mat)].chrom) << "\t" << MATCH::getOffset(r1.mat) + ref.metaWindows[MATCH::getMetaID(r1.mat)].startPos << "\t" << (MATCH::isFwd(r1.mat) ? "fwd" : "rev") << "\t" << (ref.metaWindows[MATCH::getMetaID(r1.mat)].startInd == MyConst::CPGDUMMY ? "F" : "T");
// 				of << r2.id << "\n\t" << static_cast<uint64_t>(ref.metaWindows[MATCH::getMetaID(r2.mat)].chrom) << "\t" << MATCH::getOffset(r2.mat) + ref.metaWindows[MATCH::getMetaID(r2.mat)].startPos << "\t" << (MATCH::isFwd(r2.mat) ? "fwd" : "rev") << "\t" << (ref.metaWindows[MATCH::getMetaID(r2.mat)].startInd == MyConst::CPGDUMMY ? "F" : "T");
// 			}

            if (mat1OriginalStrand)
            {
				if (getStranded)
				{
#pragma omp atomic
					++r1FwdMatches;
				}
                computeMethLvl(r1.mat, r1.seq);
                computeMethLvl(r2.mat, revSeq2);

            } else {

				if (getStranded)
				{
#pragma omp atomic
					++r1RevMatches;
				}
                computeMethLvl(r1.mat, revSeq1);
                computeMethLvl(r2.mat, r2.seq);

            }
			if (ref.metaWindows[MATCH::getMetaID(r1.mat)].startInd == MyConst::CPGDUMMY && ref.metaWindows[MATCH::getMetaID(r2.mat)].startInd == MyConst::CPGDUMMY)
			{
				tooShortCount += 2;
			} else {
				if (ref.metaWindows[MATCH::getMetaID(r1.mat)].startInd == MyConst::CPGDUMMY || ref.metaWindows[MATCH::getMetaID(r2.mat)].startInd == MyConst::CPGDUMMY)
				{
					++tooShortCount;
					++succMatchT;
				} else {

					++succPairedMatchT;
					succMatchT += 2;
				}
			}
        }
// #pragma omp critical
// {
// 		of << "---------------------------\n\n\n\n";
// }
    }

    // sum up counts
    for (unsigned int i = 0; i < CORENUM; ++i)
    {
        succMatch += matchStats[i];
        nonUniqueMatch += nonUniqueStats[i];
        unSuccMatch += noMatchStats[i];
        succPairedMatch += matchPairedStats[i];
		tooShortCountMatch += tooShortCounts[i];
    }
    return true;
}




bool ReadQueue::matchSCBatch(const char* scFile, const std::string scId, const bool isGZ)
{
    unsigned int readCounter = 0;
    unsigned int i = 0;
    // counter
    uint64_t succPairedMatch = 0;
    uint64_t succMatch = 0;
    uint64_t nonUniqueMatch = 0;
    uint64_t unSuccMatch = 0;
	uint64_t tooShortCount = 0;

    if (isGZ)
    {

        igz.open(scFile);

    } else {

        file.open(scFile);
    }
	r1FwdMatches = 0;
	r1RevMatches = 0;


    // fill counting structure for parallelization
    for (unsigned int j = 0; j < CORENUM; ++j)
    {

        fwdMetaIDs[j] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        revMetaIDs[j] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        fwdMetaIDs[j].set_deleted_key(ref.metaWindows.size() + 10);
        revMetaIDs[j].set_deleted_key(ref.metaWindows.size() + 10);
        fwdMetaIDs[j].set_empty_key(ref.metaWindows.size() + 11);
        revMetaIDs[j].set_empty_key(ref.metaWindows.size() + 11);

    }
	if (!bothStrandsFlag)
	{
		++i;
		isGZ ? parseChunkGZ(readCounter) : parseChunk(readCounter);
		matchReads(readCounter, succMatch, nonUniqueMatch, unSuccMatch, true);
		decideStrand();
        std::cout << "Processed " << MyConst::CHUNKSIZE * (i) << " paired reads\n";
	}

    while(isGZ ? parseChunkGZ(readCounter) : parseChunk(readCounter))
    {
        ++i;
		// if (i>2)
		// 	break;
        matchReads(readCounter, succMatch, nonUniqueMatch, unSuccMatch, false);
        std::cout << "Processed " << MyConst::CHUNKSIZE * (i) << " paired reads\n";
    }
    // match remaining reads
    matchReads(readCounter, succMatch, nonUniqueMatch, unSuccMatch, false);
	std::cout << "Processed " << MyConst::CHUNKSIZE * (i+1) << " paired reads\n";

	std::cout << "Finished " << std::string(scFile) << "\n\n";
	std::cout << "\nOverall number of reads for cell " << scId << ": (2*)" << MyConst::CHUNKSIZE * i + readCounter;
    std::cout << "\tOverall successfully matched: " << succMatch << "\n\tUnsuccessfully matched: " << unSuccMatch << "\n\tNonunique matches: " << nonUniqueMatch << "\n\nReads discarded as too short: " << tooShortCount << "\n\nFully matched pairs: " << succPairedMatch << "\n";

    if (isGZ)
    {

        igz.close();

    } else {

        file.close();
    }
	printSCMethylationLevels(scId);
	return true;
}
bool ReadQueue::matchSCBatchPaired(const char* scFile1, const char* scFile2, const std::string scId, const bool isGZ)
{
    unsigned int readCounter = 0;
    unsigned int i = 0;
    // counter
    uint64_t succPairedMatch = 0;
    uint64_t succMatch = 0;
    uint64_t nonUniqueMatch = 0;
    uint64_t unSuccMatch = 0;
	uint64_t tooShortCount = 0;

    if (isGZ)
    {

        igz.open(scFile1);
        igz2.open(scFile2);

    } else {

        file.open(scFile1);
        file2.open(scFile2);
    }

    // fill counting structure for parallelization
    for (unsigned int j = 0; j < CORENUM; ++j)
    {

		paired_fwdMetaIDs[j] = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>();
		paired_revMetaIDs[j] = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>();
    }
	if (!bothStrandsFlag)
	{
		++i;
		isGZ ? parseChunkGZ(readCounter) : parseChunk(readCounter);
		matchPairedReads(readCounter, succMatch, nonUniqueMatch, unSuccMatch, succPairedMatch, tooShortCount, true);
		decideStrand();
        std::cout << "Processed " << MyConst::CHUNKSIZE * (i) << " paired reads\n";
	}

    while(isGZ ? parseChunkGZ(readCounter) : parseChunk(readCounter))
    {
        ++i;
		// if (i>2)
		// 	break;
        matchPairedReads(readCounter, succMatch, nonUniqueMatch, unSuccMatch, succPairedMatch, tooShortCount, false);
        std::cout << "Processed " << MyConst::CHUNKSIZE * (i) << " paired reads\n";
    }
    // match remaining reads
    matchPairedReads(readCounter, succMatch, nonUniqueMatch, unSuccMatch, succPairedMatch, tooShortCount, false);
	std::cout << "Processed " << MyConst::CHUNKSIZE * (i+1) << " paired reads\n";

	std::cout << "Finished " << std::string(scFile1) << " + " << std::string(scFile2) << "\n\n";
	std::cout << "\nOverall number of reads for cell " << scId << ": (2*)" << MyConst::CHUNKSIZE * i + readCounter;
    std::cout << "\tOverall successfully matched: " << succMatch << "\n\tUnsuccessfully matched: " << unSuccMatch << "\n\tNonunique matches: " << nonUniqueMatch << "\n\nReads discarded as too short: " << tooShortCount << "\n\nFully matched pairs: " << succPairedMatch << "\n";

    if (isGZ)
    {

        igz.close();
        igz2.close();

    } else {

        file.close();
        file2.close();
    }
	printSCMethylationLevels(scId);
	return true;
}





void ReadQueue::printMethylationLevels(std::string& filename)
{

    // std::cout << "\nStart writing Methylation levels to \"" << filename << "_cpg_start.tsv\" and \"" << filename << "_cpg.tsv\"\n\n";
    // std::ofstream cpgFile(filename + "_cpg_start.tsv");
    //
    // // go over CpGs close to start of a chromosome
    // for (size_t cpgID = 0; cpgID < ref.cpgStartTable.size(); ++cpgID)
    // {
    //
    //     // print the position of the (C of the) CpG
    //     cpgFile << ref.chrMap[ref.cpgStartTable[cpgID].chrom] << "\t" << ref.cpgStartTable[cpgID].pos << "\t";
    //
    //     // print the counts
    //     // fwd counts
    //     cpgFile << methLevelsStart[cpgID].methFwd << "\t" << methLevelsStart[cpgID].unmethFwd << "\t";
    //     // rev counts
    //     cpgFile << methLevelsStart[cpgID].methRev << "\t" << methLevelsStart[cpgID].unmethRev << "\n";
    // }
    // cpgFile.close();

    std::cout << "\nStart writing Methylation levels to \"" << filename << "_cpg.tsv\"\n\n";
	std::ofstream cpgFile(filename + "_cpg.tsv");

    // go over remaining CpGs
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {

		// if (methLevels[cpgID].methFwd + methLevels[cpgID].unmethFwd == 0)
		// {
		// 	of << "Unmapped CpG fwd: " << ref.chrMap[ref.cpgTable[cpgID].chrom] << "\t" << ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 << "\n";
		// 	of << std::string(ref.fullSeq[ref.cpgTable[cpgID].chrom].begin() + ref.cpgTable[cpgID].pos, std::min(ref.fullSeq[ref.cpgTable[cpgID].chrom].begin() + ref.cpgTable[cpgID].pos + 2*MyConst::READLEN, ref.fullSeq[ref.cpgTable[cpgID].chrom].end())) << "\n\n";
		// }
		// if (methLevels[cpgID].methRev + methLevels[cpgID].unmethRev == 0)
		// {
		// 	of << "Unmapped CpG rev: " << ref.chrMap[ref.cpgTable[cpgID].chrom] << "\t" << ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 << "\n";
		// 	of << std::string(ref.fullSeq[ref.cpgTable[cpgID].chrom].begin() + ref.cpgTable[cpgID].pos, std::min(ref.fullSeq[ref.cpgTable[cpgID].chrom].begin() + ref.cpgTable[cpgID].pos + 2*MyConst::READLEN, ref.fullSeq[ref.cpgTable[cpgID].chrom].end())) << "\n\n";
		// }
        // print the position of the (C of the) CpG
        cpgFile << ref.chrMap[ref.cpgTable[cpgID].chrom] << "\t" << ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 << "\t";

        // print the counts
        // fwd counts
        cpgFile << methLevels[cpgID].methFwd << "\t" << methLevels[cpgID].unmethFwd << "\t";
        // rev counts
        cpgFile << methLevels[cpgID].methRev << "\t" << methLevels[cpgID].unmethRev << "\n";
    }
    cpgFile.close();
    std::cout << "Finished writing methylation levels to file\n\n";
}
void ReadQueue::printSCMethylationLevels(const std::string scID)
{

	scOutput << "\n" << scID << "\tmethFwd\t";
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {
        scOutput << methLevels[cpgID].methFwd << "\t";
    }
	scOutput << "\n" << scID << "\tunmethFwd\t";
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {
        scOutput << methLevels[cpgID].unmethFwd << "\t";
    }
	scOutput << "\n" << scID << "\tmethRev\t";
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {
        scOutput << methLevels[cpgID].methRev << "\t";
    }
	scOutput << "\n" << scID << "\tunmethFwd\t";
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {
        scOutput << methLevels[cpgID].unmethRev << "\t";
    }
}






inline int ReadQueue::saQuerySeedSetRef(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, MATCH::match& mat, uint16_t& qThreshold)
{

	// use counters to flag what has been processed so far
	auto& fwdMetaIDs_t = fwdMetaIDs[omp_get_thread_num()];
	auto& revMetaIDs_t = revMetaIDs[omp_get_thread_num()];

	// counter for how often we had a match
	std::array<uint8_t, MyConst::ADDMIS + MyConst::MISCOUNT + 1> multiMatch;
	multiMatch.fill(0);

	// will contain matches iff match is found for number of errors specified by index
	std::array<MATCH::match, MyConst::ADDMIS + MyConst::MISCOUNT + 1> uniqueMatches;
	// store the last match found in current MetaCpG
	uint8_t prevChr = 0;
	uint64_t prevOff = 0xffffffffffffffffULL;


	// check all fwd meta CpGs
	for (const auto& m : fwdMetaIDs_t)
	{
		// apply qgram lemma
		if (m.second < qThreshold)
			continue;

		auto startIt = ref.fullSeq[ref.metaWindows[m.first].chrom].begin() + ref.metaWindows[m.first].startPos;
		auto endIt = ref.fullSeq[ref.metaWindows[m.first].chrom].begin() + ref.metaWindows[m.first].startPos + MyConst::WINLEN + MyConst::MISCOUNT + MyConst::ADDMIS - 1;

		// check if CpG was too near to the end
		if (endIt > ref.fullSeq[ref.metaWindows[m.first].chrom].end())
		{
			endIt = ref.fullSeq[ref.metaWindows[m.first].chrom].end();
		}

		// use shift and to find all matchings
		std::vector<uint64_t> matchings;
		std::vector<uint8_t> errors;
		sa.querySeq(startIt, endIt, matchings, errors);

		size_t i = 0;
		// compare first found match with last found match of previous meta CpG
		if (matchings.size() > 0)
		{
			// compare chromosome and offset
			if (matchings[0] + ref.metaWindows[m.first].startPos == prevOff && ref.metaWindows[m.first].chrom == prevChr)
			{
				++i;
			}
		}
		// go through matching and see if we had such a match (with that many errors) before - if so,
		// return to caller reporting no match
		for (; i < matchings.size(); ++i)
		{

			// check if we had a match with that many errors before
			if (multiMatch[errors[i]])
			{

				MATCH::match& match_2 = uniqueMatches[errors[i]];
				const bool isStart = MATCH::isStart(match_2);
				const bool isFwd = MATCH::isFwd(match_2);
				// check if same k-mer (borders of meta CpGs)
				if (isFwd && !isStart && ref.metaWindows[MATCH::getMetaID(match_2)].startPos + MATCH::getOffset(match_2) == ref.metaWindows[m.first].startPos + matchings[i])
				{
					continue;

				} else {

					// check if this is a match without errors
					if (!errors[i])
					{

						// if so, return without a match
						return -1;

					}
					// set the number of matches with that many errors to 2
					// indicating that we do not have a unique match with that many errors
					multiMatch[errors[i]] = 2;
				}


			} else {

				// update qgram lemma
				// uint16_t newQ = sa.size() - MyConst::KMERLEN - (MyConst::KMERLEN * errors[i]);
				// check for overflow and if we improved old q
				// if (newQ < sa.size() && newQ > qThreshold)
				// 	qThreshold = newQ;


				// we don't have such a match yet,
				// so save this match at the correct position
				uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 1, 0, m.first);
				multiMatch[errors[i]] = 1;
			}
		}
		if (matchings.size() > 0)
		{

			prevChr = ref.metaWindows[m.first].chrom;
			prevOff = ref.metaWindows[m.first].startPos + matchings[matchings.size() - 1];

		} else {

			prevChr = 0;
			prevOff = 0xffffffffffffffffULL;
		}
	}
	prevChr = 0;
	prevOff = 0xffffffffffffffffULL;
	// go through reverse sequences
	for (const auto& m : revMetaIDs_t)
	{

		// apply qgram lemma
		if (m.second < qThreshold)
			continue;

		// retrieve sequence
		// retrieve sequence
		auto endIt = ref.fullSeq[ref.metaWindows[m.first].chrom].begin();

		if (ref.metaWindows[m.first].startPos > 0)
		{
			endIt += ref.metaWindows[m.first].startPos - 1;
		}

		auto startIt = ref.fullSeq[ref.metaWindows[m.first].chrom].begin() + ref.metaWindows[m.first].startPos + MyConst::WINLEN + MyConst::MISCOUNT + MyConst::ADDMIS - 1;
		if (startIt >= ref.fullSeq[ref.metaWindows[m.first].chrom].end())
		{
			startIt = ref.fullSeq[ref.metaWindows[m.first].chrom].end() - 1;
		}

		// use shift and to find all matchings
		std::vector<uint64_t> matchings;
		std::vector<uint8_t> errors;
		sa.queryRevSeq(startIt, endIt, matchings, errors);

		size_t i = 0;
		// compare first found match with last found match of previous meta CpG
		if (matchings.size() > 0)
		{
			// compare chromosome and offset
			if (matchings[0] + ref.metaWindows[m.first].startPos == prevOff && ref.metaWindows[m.first].chrom == prevChr)
			{
				++i;
			}
		}
		// go through matching and see if we had such a match (with that many errors) before - if so,
		// return to caller reporting no match
		for (; i < matchings.size(); ++i)
		{

			// check if we had a match with that many errors before
			if (multiMatch[errors[i]])
			{

				MATCH::match& match_2 = uniqueMatches[errors[i]];
				const bool isStart = MATCH::isStart(match_2);
				const bool isFwd = MATCH::isFwd(match_2);
				// check if same k-mer (borders of meta CpGs)
				if (!isFwd && !isStart && ref.metaWindows[MATCH::getMetaID(match_2)].startPos + MATCH::getOffset(match_2) == ref.metaWindows[m.first].startPos + matchings[i])
				{
					continue;

				} else {

					// check if this is a match without errors
					if (!errors[i])
					{

						// if so, return without a match
						return -1;

					}
					// set the number of matches with that many errors to 2
					// indicating that we do not have a unique match with that many errors
					multiMatch[errors[i]] = 2;
				}


			} else {

				// update qgram lemma
				// uint16_t newQ = sa.size() - MyConst::KMERLEN - (MyConst::KMERLEN * errors[i]);
				// check for overflow and if we improved old q
				// if (newQ < sa.size() && newQ > qThreshold)
				// 	qThreshold = newQ;

				// we don't have such a match yet,
				// so save this match at the correct position
				uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 0, 0, m.first);
				multiMatch[errors[i]] = 1;
			}
		}
		if (matchings.size() > 0)
		{

			prevChr = ref.metaWindows[m.first].chrom;
			prevOff = ref.metaWindows[m.first].startPos + matchings[matchings.size() - 1];

		} else {

			prevChr = 0;
			prevOff = 0xffffffffffffffffULL;
		}
	}

	// go through found matches for each [0,maxErrorNumber] and see if it is unique
	for (size_t i = 0; i < multiMatch.size(); ++i)
	{
		// there is no match with that few errors, search the one with more errors
		if (multiMatch[i] == 0)
		{
			continue;
		}
		mat = uniqueMatches[i];
		// if match is not unique, return unsuccessfull to caller
		if (multiMatch[i] > 1)
		{

			return -1;

		// exactly one with that many errors - return successfull
		} else {

			return 1;
		}

	}
	// we have not a single match at all, return unsuccessfull to caller
	return 0;
}






inline void ReadQueue::saQuerySeedSetRefFirst(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, std::vector<MATCH::match>& mats, const uint16_t& qThreshold)
{

	// use counters to flag what has been processed so far
	auto& fwdMetaIDs_t = paired_fwdMetaIDs[omp_get_thread_num()];
	auto& revMetaIDs_t = paired_revMetaIDs[omp_get_thread_num()];

	// store the last match found in current MetaCpG
	uint8_t prevChr = 0;
	uint64_t prevOff = 0xffffffffffffffffULL;

	std::vector<std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> > > fwdMetas(fwdMetaIDs_t.begin(), fwdMetaIDs_t.end());
	std::sort(fwdMetas.begin(), fwdMetas.end(), cmpMetaFirst);
	std::vector<std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> > > revMetas(revMetaIDs_t.begin(), revMetaIDs_t.end());
	std::sort(revMetas.begin(), revMetas.end(), cmpMetaFirst);
	// current best k-mer count of true match
	int32_t bmCount = 0;
	if (fwdMetas.size() > 0)
	{
		if (revMetas.size() > 0)
		{
			// match fwd first
			if (std::get<0>(fwdMetas[0].second) > std::get<0>(revMetas[0].second))
			{
				// check all fwd meta CpGs
				for (auto& meta : fwdMetas)
				{
					if (!matchFwdFirst(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}
				prevChr = 0;
				prevOff = 0xffffffffffffffffULL;
				// check all rev meta CpGs
				for (auto& meta : revMetas)
				{
					if (!matchRevFirst(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}

			// match rev first
			} else {
				// check all rev meta CpGs
				for (auto& meta : revMetas)
				{
					if (!matchRevFirst(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}
				prevChr = 0;
				prevOff = 0xffffffffffffffffULL;
				// check all fwd meta CpGs
				for (auto& meta : fwdMetas)
				{
					if (!matchFwdFirst(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}
			}
		// Just match fwd
		} else {
			for (auto& meta : fwdMetas)
			{
				if (!matchFwdFirst(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
					break;
			}
		}

	} else {
		// Just match rev
		if (revMetas.size() > 0)
		{
			for (auto& meta : revMetas)
			{
				if (!matchRevFirst(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
					break;
			}
		}
	}
}





inline void ReadQueue::saQuerySeedSetRefSecond(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, std::vector<MATCH::match>& mats, const uint16_t& qThreshold)
{
	auto& fwdMetaIDs_t = paired_fwdMetaIDs[omp_get_thread_num()];
	auto& revMetaIDs_t = paired_revMetaIDs[omp_get_thread_num()];

	// store the last match found in current MetaCpG
	uint8_t prevChr = 0;
	uint64_t prevOff = 0xffffffffffffffffULL;

	std::vector<std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> > > fwdMetas(fwdMetaIDs_t.begin(), fwdMetaIDs_t.end());
	std::sort(fwdMetas.begin(), fwdMetas.end(), cmpMetaSecond);
	std::vector<std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> > > revMetas(revMetaIDs_t.begin(), revMetaIDs_t.end());
	std::sort(revMetas.begin(), revMetas.end(), cmpMetaSecond);
	// current best k-mer count of true match
	int32_t bmCount = 0;
	if (fwdMetas.size() > 0)
	{
		if (revMetas.size() > 0)
		{
			// match fwd first
			if (std::get<0>(fwdMetas[0].second) > std::get<0>(revMetas[0].second))
			{
				// check all fwd meta CpGs
				for (auto& meta : fwdMetas)
				{
					if (!matchFwdSecond(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}
				prevChr = 0;
				prevOff = 0xffffffffffffffffULL;
				// check all rev meta CpGs
				for (auto& meta : revMetas)
				{
					if (!matchRevSecond(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}

			// match rev first
			} else {
				// check all rev meta CpGs
				for (auto& meta : revMetas)
				{
					if (!matchRevSecond(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}
				prevChr = 0;
				prevOff = 0xffffffffffffffffULL;
				// check all fwd meta CpGs
				for (auto& meta : fwdMetas)
				{
					if (!matchFwdSecond(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
						break;
				}
			}
		// Just match fwd
		} else {
			for (auto& meta : fwdMetas)
			{
				if (!matchFwdSecond(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
					break;
			}
		}

	} else {
		// Just match rev
		if (revMetas.size() > 0)
		{
			for (auto& meta : revMetas)
			{
				if (!matchRevSecond(meta, prevChr, prevOff, mats, bmCount, qThreshold, sa))
					break;
			}
		}
	}
}





inline void ReadQueue::getSeedRefs(const std::string& seq, const size_t& readSize, const uint16_t qThreshold)
{

	auto& fwdMetaIDs_t = fwdMetaIDs[omp_get_thread_num()];
	auto& revMetaIDs_t = revMetaIDs[omp_get_thread_num()];
	fwdMetaIDs_t.clear();
	revMetaIDs_t.clear();

	// retrieve kmers for first hash
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seq.data(), MyConst::SEED, MyConst::KMERLEN, fhVal);

	uint64_t key = sfVal % MyConst::HTABSIZE;

	uint64_t lastId = 0xffffffffffffffffULL;
	bool wasFwd = false;

	// maximum position until we can insert completely new meta cpgs
	uint32_t maxQPos = seq.size() - MyConst::KMERLEN + 1 - qThreshold;

	uint64_t endIdx = ref.tabIndex[key+1];
	for (uint64_t i = ref.tabIndex[key]; i < endIdx; ++i)
	{

		const uint32_t metaId = KMER_S::getMetaCpG(ref.kmerTableSmall[i]);
		const bool isFwd = ref.strandTable[i];
		// check if we visited meta CpG before
		if (metaId == lastId && isFwd == wasFwd)
		{
			continue;
		}

		// update vars for last checked metaCpG
		lastId = metaId;
		wasFwd = isFwd;

		if (isFwd)
		{
			++fwdMetaIDs_t[metaId];

		} else {

			++revMetaIDs_t[metaId];

		}
	}

	for (unsigned int cIdx = 0; cIdx < (seq.size() - MyConst::KMERLEN); ++cIdx)
	{

		// use rolling hash
		sfVal = ntHash::NTPS64(seq.data()+cIdx+1, MyConst::SEED, seq[cIdx], seq[cIdx + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);

		key = sfVal % MyConst::HTABSIZE;

		lastId = 0xffffffffffffffffULL;
		wasFwd = false;

		endIdx = ref.tabIndex[key+1];
		for (uint64_t i = ref.tabIndex[key]; i < endIdx; ++i)
		{

			const uint32_t metaId = KMER_S::getMetaCpG(ref.kmerTableSmall[i]);
			const bool isFwd = ref.strandTable[i];
			// check if we visited meta CpG before
			if (metaId == lastId && isFwd == wasFwd)
			{
				continue;
			}

			// update vars for last checked metaCpG
			lastId = metaId;
			wasFwd = isFwd;

			if (isFwd)
			{
				// check if it is at all possible to have newly inserted element passing q
				if (cIdx < maxQPos)
				{
					++fwdMetaIDs_t[metaId];

				} else {

					auto it = fwdMetaIDs_t.find(metaId);
					if (it != fwdMetaIDs_t.end())
					{
						++it->second;
					}
				}

			} else {

				if (cIdx < maxQPos)
				{
					++revMetaIDs_t[metaId];

				} else {

					auto it = revMetaIDs_t.find(metaId);
					if (it != revMetaIDs_t.end())
					{
						++it->second;
					}
				}

			}
		}
	}
}






inline uint16_t ReadQueue::getSeedRefsFirstRead(const std::string& seq, const size_t& readSize, const uint16_t qThreshold)
{

	// std::vector<uint16_t>& threadCountFwdStart = countsFwdStart[omp_get_thread_num()];
	// std::vector<uint16_t>& threadCountRevStart = countsRevStart[omp_get_thread_num()];
	// fill with zeroes
	// threadCountFwdStart.assign(ref.metaStartCpGs.size(), 0);
	// threadCountRevStart.assign(ref.metaStartCpGs.size(), 0);

	auto& fwdMetaIDs_t = paired_fwdMetaIDs[omp_get_thread_num()];
	auto& revMetaIDs_t = paired_revMetaIDs[omp_get_thread_num()];
	// fwdMetaIDs_t.clear();
	// revMetaIDs_t.clear();
	// fwdMetaIDs_t.resize(800);
	// revMetaIDs_t.resize(800);
	uint32_t kmerNum = 1;
	uint32_t bucketCount = 0;
	{
		// retrieve kmers for first hash
		uint64_t fhVal;
		uint64_t sfVal = ntHash::NTPS64(seq.data(), MyConst::SEED, MyConst::KMERLEN, fhVal);

		uint64_t key = sfVal % MyConst::HTABSIZE;
		bucketCount += ref.tabIndex[key+1] - ref.tabIndex[key];

		for (unsigned int cIdx = 0; cIdx < (seq.size() - MyConst::KMERLEN); ++cIdx)
		{

			// use rolling hash
			sfVal = ntHash::NTPS64(seq.data()+cIdx+1, MyConst::SEED, seq[cIdx], seq[cIdx + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);

			key = sfVal % MyConst::HTABSIZE;
			bucketCount += ref.tabIndex[key+1] - ref.tabIndex[key];
			++kmerNum;
		}
	}
// #pragma omp critical
// 	{
// 	of << bucketCount << "\n";
// 	}
	bucketCount = bucketCount / 5;
	auto newmap = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>(bucketCount);
	fwdMetaIDs_t.swap(newmap);
	newmap = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>(bucketCount);
	revMetaIDs_t.swap(newmap);

	// retrieve kmers for first hash
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seq.data(), MyConst::SEED, MyConst::KMERLEN, fhVal);

	uint64_t key = sfVal % MyConst::HTABSIZE;

	uint64_t lastId = 0xffffffffffffffffULL;
	bool wasFwd = false;
	bool wasStart = false;

	// compute bitmask for all positions where C occurrs
	uint32_t cMask = 0;
	// also compute perfect hash
	uint32_t kSeq = 0;
	for (unsigned int i = 0; i < MyConst::KMERLEN; ++i)
	{
		cMask = cMask << 1;
		if (seq[i] == 'C')
		{
			cMask |= 1;
		}
		kSeq = kSeq << 2;
		switch (seq[i])
		{

			case 'C':
			case 'T':

				kSeq += 3;
				break;

			case 'G':

				kSeq += 2;
				break;
		// end switch
		}
	}
	// uint16_t qAdapt = 0;
	// maximum position until we can insert completely new meta cpgs
	uint32_t maxQPos = seq.size() - MyConst::KMERLEN + 1 - qThreshold;
	// if (ref.filteredKmers.find(kSeq) != ref.filteredKmers.end())
	// {
	// 	++maxQPos;
	// 	++qAdapt;
	// }


	uint64_t endIdx = ref.tabIndex[key+1];
	for (uint64_t i = ref.tabIndex[key]; i < endIdx; ++i)
	{

		const KMER_S::kmer currentKmer = ref.kmerTableSmall[i];
		const uint32_t metaId = KMER_S::getMetaCpG(currentKmer);
		// test for asymmetric mapping
		if (cMask & MyConst::SEEDBITS & currentKmer.tmask)
		{
			continue;
		}
		const bool isFwd = ref.strandTable[i];
		const bool isStart = KMER_S::isStartCpG(currentKmer);
		// check if we visited meta CpG before
		if (metaId == lastId && isFwd == wasFwd && isStart == wasStart)
		{
			continue;
		}

		// update vars for last checked metaCpG
		lastId = metaId;
		wasFwd = isFwd;
		wasStart = isStart;
		// if (isStart)
		// {
		// 	if (isFwd)
		// 	{
		// 		++threadCountFwdStart[metaId];
        //
		// 	} else {
        //
		// 		++threadCountRevStart[metaId];
        //
		// 	}
        //
		// } else {

			if (isFwd)
			{
				++std::get<0>(fwdMetaIDs_t[metaId]);

			} else {

				++std::get<0>(revMetaIDs_t[metaId]);

			}
		// }
	}

	for (unsigned int cIdx = 0; cIdx < (seq.size() - MyConst::KMERLEN); ++cIdx)
	{

		kSeq = kSeq << 2;
		switch (seq[MyConst::KMERLEN + cIdx])
		{

			case 'C':
			case 'T':

				kSeq += 3;
				break;

			case 'G':

				kSeq += 2;
				break;
		// end switch
		}
		// TODO: check if correct for k-mer lengths < 32
		cMask = (cMask << 1) & MyConst::KMERMASK;
		if (seq[MyConst::KMERLEN + cIdx] == 'C')
		{
			cMask |= 1;
		}

		// use rolling hash
		sfVal = ntHash::NTPS64(seq.data()+cIdx+1, MyConst::SEED, seq[cIdx], seq[cIdx + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);

		key = sfVal % MyConst::HTABSIZE;

		// test if kmer is blacklisted
		if (ref.filteredKmers.find(kSeq) != ref.filteredKmers.end())
		{
			// ++maxQPos;
			// ++qAdapt;
			// if (qAdapt >= qThreshold)
			// {
			// 	return 0;
			// }
			continue;
		}
		lastId = 0xffffffffffffffffULL;
		wasFwd = false;
		wasStart = false;

		endIdx = ref.tabIndex[key+1];
		for (uint64_t i = ref.tabIndex[key]; i < endIdx; ++i)
		{

			const KMER_S::kmer currentKmer = ref.kmerTableSmall[i];
			const uint32_t metaId = KMER_S::getMetaCpG(currentKmer);
			if (cMask & MyConst::SEEDBITS & currentKmer.tmask)
			{
				continue;
			}
			const bool isFwd = ref.strandTable[i];
			const bool isStart = KMER_S::isStartCpG(currentKmer);
			// check if we visited meta CpG before
			if (metaId == lastId && isFwd == wasFwd && isStart == wasStart)
			{
				continue;
			}

			// update vars for last checked metaCpG
			lastId = metaId;
			wasFwd = isFwd;
			wasStart = isStart;
			// if (isStart)
			// {
            //
			// 	if (isFwd)
			// 	{
			// 		++threadCountFwdStart[metaId];
            //
			// 	} else {
            //
			// 		++threadCountRevStart[metaId];
            //
			// 	}
            //
			// } else {

				if (isFwd)
				{
					// check if it is at all possible to have newly inserted element passing q
					if (cIdx < maxQPos)
					{
						++std::get<0>(fwdMetaIDs_t[metaId]);

					} else {

						auto it = fwdMetaIDs_t.find(metaId);
						if (it != fwdMetaIDs_t.end())
						{
							++std::get<0>(it.value());
						}
					}

				} else {

					if (cIdx < maxQPos)
					{
						++std::get<0>(revMetaIDs_t[metaId]);

					} else {

						auto it = revMetaIDs_t.find(metaId);
						if (it != revMetaIDs_t.end())
						{
							++std::get<0>(it.value());
						}
					}

				}
			// }
		}
	}
	// if (qAdapt < qThreshold)
	// {
	// 	qAdapt = std::max(5, qThreshold - qAdapt);
	// } else {
	// 	qAdapt = 5;
	// 	// return 0;
	// }
	// TODO
	// std::cout << "\n\n\nFirst read candidates:\n\n";
	// for (auto hIt = fwdMetaIDs_t.begin(); hIt != fwdMetaIDs_t.end(); ++hIt)
	// {
	// 	std::cout << ref.chrMap[ref.cpgTable[ref.metaCpGs[hIt->first].start].chrom] << "\t";
	// 	std::cout << ref.cpgTable[ref.metaCpGs[hIt->first].start].pos << "\t";
	// 	std::cout << hIt->first << "\t";
	// 	std::cout << (uint64_t)std::get<0>(hIt->second) << "\t";
	// 	std::cout << (uint64_t)std::get<1>(hIt->second) << "\t";
	// 	std::cout << std::get<2>(hIt->second) << "\n";
	// }
	// for (auto hIt = revMetaIDs_t.begin(); hIt != revMetaIDs_t.end(); ++hIt)
	// {
	// 	std::cout << ref.chrMap[ref.cpgTable[ref.metaCpGs[hIt->first].start].chrom] << "\t";
	// 	std::cout << ref.cpgTable[ref.metaCpGs[hIt->first].start].pos << "\t";
	// 	std::cout << hIt->first << "\t";
	// 	std::cout << (uint64_t)std::get<0>(hIt->second) << "\t";
	// 	std::cout << (uint64_t)std::get<1>(hIt->second) << "\t";
	// 	std::cout << std::get<2>(hIt->second) << "\n";
	// }
	// of << "Estimate: " << bucketCount << "\n";
	// of << "Real size: " << fwdMetaIDs_t.size() << "/" << revMetaIDs_t.size() << "\n\n\n"; 
	auto newFwdMetas = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>(20);
	auto newRevMetas = tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>(20);
	// unsigned int c = 0;
	for (auto hIt = fwdMetaIDs_t.begin(); hIt != fwdMetaIDs_t.end(); ++hIt)
	{
		// of << hIt->first << "/" << ref.chrMap[ref.metaWindows[hIt->first].chrom] << "/" << ref.metaWindows[hIt->first].startPos << "\t" << (uint32_t)std::get<0>(hIt->second) << "\n";
		if (std::get<0>(hIt->second) >= qThreshold)
		{
			newFwdMetas[hIt->first] = {std::get<0>(hIt->second), 0, false, false};
		}
		// if (std::get<0>(hIt->second) < qThreshold)
		// 	fwdMetaIDs_t.erase(hIt);
		// else
		// 	++c;
	}
	for (auto hIt = revMetaIDs_t.begin(); hIt != revMetaIDs_t.end(); ++hIt)
	{
		// of << hIt->first << "/" << ref.chrMap[ref.metaWindows[hIt->first].chrom] << "/" << ref.metaWindows[hIt->first].startPos << "\t" << (uint32_t)std::get<0>(hIt->second) << "\n";
		if (std::get<0>(hIt->second) >= qThreshold)
		{
			newRevMetas[hIt->first] = {std::get<0>(hIt->second), 0, false, false};
		}
		// if (std::get<0>(hIt->second) < qThreshold)
		// 	revMetaIDs_t.erase(hIt);
		// else
		// 	++c;
	}
	// if (c == 0)
	// 	return;
	// fwdMetaIDs_t.rehash(1);
	// revMetaIDs_t.rehash(1);
	fwdMetaIDs_t.swap(newFwdMetas);
	revMetaIDs_t.swap(newRevMetas);
	// return qAdapt;
	//
	return 0;
}






inline bool ReadQueue::getSeedRefsSecondRead(const std::string& seq, const size_t& readSize, const uint16_t qThreshold)
{

	auto& fwdMetaIDs_t = paired_fwdMetaIDs[omp_get_thread_num()];
	auto& revMetaIDs_t = paired_revMetaIDs[omp_get_thread_num()];

	// retrieve kmers for first hash
	uint64_t fhVal;
	uint64_t sfVal = ntHash::NTPS64(seq.data(), MyConst::SEED, MyConst::KMERLEN, fhVal);

	uint64_t key = sfVal % MyConst::HTABSIZE;

	uint64_t lastId = 0xffffffffffffffffULL;
	bool wasFwd = false;

	// compute bitmask for all positions where C occurrs
	uint32_t cMask = 0;
	// also compute perfect hash
	uint32_t kSeq = 0;
	for (unsigned int i = 0; i < MyConst::KMERLEN; ++i)
	{
		cMask = cMask << 1;
		if (seq[i] == 'C')
		{
			cMask |= 1;
		}
		kSeq = kSeq << 2;
		switch (seq[i])
		{

			case 'C':
			case 'T':

				kSeq += 3;
				break;

			case 'G':

				kSeq += 2;
				break;
		// end switch
		}
	}
	// maximum position until we can insert completely new meta cpgs
	uint32_t maxQPos = seq.size() - MyConst::KMERLEN + 1 - qThreshold;
	// uint16_t qAdapt = 0;
	// if (ref.filteredKmers.find(kSeq) != ref.filteredKmers.end())
	// {
	// 	++qAdapt;
	// 	++maxQPos;
	// }

	// how many windows we need to look for around current meta CpG
	constexpr int contextWLen = (int)(((double)MyConst::MAXPDIST / MyConst::WINLEN)) + 1;
	// TODO
	// count how often we have a meta CpG candidate for matching
	// unsigned int candCount = 0;
	uint64_t endIdx = ref.tabIndex[key+1];
	for (uint64_t i = ref.tabIndex[key]; i < endIdx; ++i)
	{

		const KMER_S::kmer currentKmer = ref.kmerTableSmall[i];
		const uint32_t metaId = KMER_S::getMetaCpG(currentKmer);
		// test for asymmetric mapping
		if (cMask & MyConst::SEEDBITS & currentKmer.tmask)
		{
			continue;
		}
		const bool isFwd = ref.strandTable[i];
		// check if we visited meta CpG before
		if (metaId == lastId && isFwd == wasFwd)
		{
			continue;
		}

		// update vars for last checked metaCpG
		lastId = metaId;
		wasFwd = isFwd;

		// test if the current or its adjacent Meta CpGs fulfill qgram lemma for the first read
		auto foundMeta = fwdMetaIDs_t.end();
		// update counts for second read
		if (isFwd)
		{
			for (int cOff = -contextWLen; cOff <= contextWLen; ++cOff)
			{
				if ((foundMeta = fwdMetaIDs_t.find(metaId + cOff)) != fwdMetaIDs_t.end())
				{
					// if (std::get<0>(foundMeta->second) >= qThreshold)
					// {
						++std::get<1>(fwdMetaIDs_t[metaId]);
						break;
						// propagate information to adjacent meta CpGs if enough kmers are matched
						// if (std::get<1>(fwdMetaIDs_t[metaId]) == qThreshold)
						// {
						// 	// ++candCount;
						// 	std::get<2>(fwdMetaIDs_t[metaId]) = true;
						// 	for (int cOffUp = 1; cOffUp <= contextWLen; ++cOffUp)
						// 	{
						// 		auto it = fwdMetaIDs_t.find(metaId + cOffUp);
						// 		if (it != fwdMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 		it = fwdMetaIDs_t.find(metaId - cOffUp);
						// 		if (it != fwdMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 	}
						// }
						// break;
					// }
				}
			}
		} else {

			for (int cOff = -contextWLen; cOff <= contextWLen; ++cOff)
			{
				if ((foundMeta = revMetaIDs_t.find(metaId + cOff)) != revMetaIDs_t.end())
				{
					// if (std::get<0>(foundMeta->second) >= qThreshold)
					// {
						++std::get<1>(revMetaIDs_t[metaId]);
						break;
						// propagate information to adjacent meta CpGs if enough kmers are matched
						// if (std::get<1>(revMetaIDs_t[metaId]) == qThreshold)
						// {
						// 	// ++candCount;
						// 	std::get<2>(revMetaIDs_t[metaId]) = true;
						// 	for (int cOffUp = 1; cOffUp <= contextWLen; ++cOffUp)
						// 	{
						// 		auto it = revMetaIDs_t.find(metaId + cOffUp);
						// 		if (it != revMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 		it = revMetaIDs_t.find(metaId - cOffUp);
						// 		if (it != revMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 	}
						// }
						// break;
					// }
				}
			}
		}
	}

	for (unsigned int cIdx = 0; cIdx < (seq.size() - MyConst::KMERLEN); ++cIdx)
	{

		kSeq = kSeq << 2;
		switch (seq[MyConst::KMERLEN + cIdx])
		{

			case 'C':
			case 'T':

				kSeq += 3;
				break;

			case 'G':

				kSeq += 2;
				break;
		// end switch
		}
		cMask = (cMask << 1) & MyConst::KMERMASK;
		if (seq[MyConst::KMERLEN + cIdx] == 'C')
		{
			cMask |= 1;
		}
		// use rolling hash
		sfVal = ntHash::NTPS64(seq.data()+cIdx+1, MyConst::SEED, seq[cIdx], seq[cIdx + MyConst::KMERLEN], MyConst::KMERLEN, fhVal);

		key = sfVal % MyConst::HTABSIZE;

		// test if kmer is blacklisted
		if (ref.filteredKmers.find(kSeq) != ref.filteredKmers.end())
		{
			// ++maxQPos;
			// ++qAdapt;
			// if (qAdapt >= qThreshold)
			// {
			// 	return 0;
			// }
			continue;
		}
		lastId = 0xffffffffffffffffULL;
		wasFwd = false;

		endIdx = ref.tabIndex[key+1];
		for (uint64_t i = ref.tabIndex[key]; i < endIdx; ++i)
		{

			const KMER_S::kmer currentKmer = ref.kmerTableSmall[i];
			const uint32_t metaId = KMER_S::getMetaCpG(currentKmer);
			if (cMask & MyConst::SEEDBITS & currentKmer.tmask)
			{
				continue;
			}
			const bool isFwd = ref.strandTable[i];
			// check if we visited meta CpG before
			if (metaId == lastId && isFwd == wasFwd)
			{
				continue;
			}

			// update vars for last checked metaCpG
			lastId = metaId;
			wasFwd = isFwd;

			// check if it is at all possible to have newly inserted element passing q
			if (cIdx < maxQPos)
			{
				// test if the current or its adjacent Meta CpGs fulfill qgram lemma for the first read
				auto foundMeta = fwdMetaIDs_t.end();
				// update counts for second read
				if (isFwd)
				{
					for (int cOff = -contextWLen; cOff <= contextWLen; ++cOff)
					{
						if ((foundMeta = fwdMetaIDs_t.find(metaId + cOff)) != fwdMetaIDs_t.end())
						{
							// if (std::get<0>(foundMeta->second) >= qThreshold)
							// {
								++std::get<1>(fwdMetaIDs_t[metaId]);
								break;
							// }
						}
					}
				} else {

					for (int cOff = -contextWLen; cOff <= contextWLen; ++cOff)
					{
						if ((foundMeta = revMetaIDs_t.find(metaId + cOff)) != revMetaIDs_t.end())
						{
							// if (std::get<0>(foundMeta->second) >= qThreshold)
							// {
								++std::get<1>(revMetaIDs_t[metaId]);
								break;
							// }
						}
					}
				}

			} else {

				if (isFwd)
				{
					auto it = fwdMetaIDs_t.find(metaId);
					if (it != fwdMetaIDs_t.end())
					{
						++std::get<1>(it.value());
						// propagate information to adjacent meta CpGs if enough kmers are matched
						// if (std::get<1>(it->second) == qThreshold)
						// {
						// 	// ++candCount;
						// 	std::get<2>(it.value()) = true;
						// 	for (int cOff = 1; cOff <= contextWLen; ++cOff)
						// 	{
						// 		it = fwdMetaIDs_t.find(metaId + cOff);
						// 		if (it != fwdMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 		it = fwdMetaIDs_t.find(metaId - cOff);
						// 		if (it != fwdMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 	}
						// }
					}
				} else {

					auto it = revMetaIDs_t.find(metaId);
					if (it != revMetaIDs_t.end())
					{
						++std::get<1>(it.value());
						// propagate information to adjacent meta CpGs if enough kmers are matched
						// if (std::get<1>(it->second) == qThreshold)
						// {
						// 	// ++candCount;
						// 	std::get<2>(it.value()) = true;
						// 	for (int cOff = 1; cOff <= contextWLen; ++cOff)
						// 	{
						// 		it = revMetaIDs_t.find(metaId + cOff);
						// 		if (it != revMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 		it = revMetaIDs_t.find(metaId - cOff);
						// 		if (it != revMetaIDs_t.end())
						// 			std::get<2>(it.value()) = true;
						// 	}
						// }
					}
				}
			}
		}
	}
	// if (qAdapt < qThreshold)
	// {
	// 	qAdapt = std::max(5, qThreshold - qAdapt);
	// } else {
    //
	// 	qAdapt = 5;
	// }
	for (auto mInfo = fwdMetaIDs_t.begin(); mInfo != fwdMetaIDs_t.end(); ++mInfo)
	{
		// propagate information to adjacent meta CpGs if enough kmers are matched
		if (std::get<1>(mInfo->second) >= qThreshold)
		{
			std::get<2>(mInfo.value()) = true;
			for (int cOffUp = 1; cOffUp <= contextWLen; ++cOffUp)
			{
				auto it = fwdMetaIDs_t.find(mInfo->first + cOffUp);
				if (it != fwdMetaIDs_t.end())
					std::get<2>(it.value()) = true;
				it = fwdMetaIDs_t.find(mInfo->first - cOffUp);
				if (it != fwdMetaIDs_t.end())
					std::get<2>(it.value()) = true;
			}
		}
	}
	for (auto mInfo = revMetaIDs_t.begin(); mInfo != revMetaIDs_t.end(); ++mInfo)
	{
		// propagate information to adjacent meta CpGs if enough kmers are matched
		if (std::get<1>(mInfo->second) >= qThreshold)
		{
			std::get<2>(mInfo.value()) = true;
			for (int cOffUp = 1; cOffUp <= contextWLen; ++cOffUp)
			{
				auto it = revMetaIDs_t.find(mInfo->first + cOffUp);
				if (it != revMetaIDs_t.end())
					std::get<2>(it.value()) = true;
				it = revMetaIDs_t.find(mInfo->first - cOffUp);
				if (it != revMetaIDs_t.end())
					std::get<2>(it.value()) = true;
			}
		}
	}
	for (auto& m : fwdMetaIDs_t)
	{
		if (std::get<2>(m.second) && ref.metaWindows[m.first].startInd != MyConst::CPGDUMMY)
			return true;
	}
	for (auto& m : revMetaIDs_t)
	{
		if (std::get<2>(m.second) && ref.metaWindows[m.first].startInd != MyConst::CPGDUMMY)
			return true;
	}
		// std::cout << "\n\n\nSecond read candidates:\n\n";
	// for (auto hIt = fwdMetaIDs_t.begin(); hIt != fwdMetaIDs_t.end(); ++hIt)
	// {
	// 	std::cout << ref.chrMap[ref.cpgTable[ref.metaCpGs[hIt->first].start].chrom] << "\t";
	// 	std::cout << ref.cpgTable[ref.metaCpGs[hIt->first].start].pos << "\t";
	// 	std::cout << hIt->first << "\t";
	// 	std::cout << (uint64_t)std::get<0>(hIt->second) << "\t";
	// 	std::cout << (uint64_t)std::get<1>(hIt->second) << "\t";
	// 	std::cout << std::get<2>(hIt->second) << "\n";
	// }
	// for (auto hIt = revMetaIDs_t.begin(); hIt != revMetaIDs_t.end(); ++hIt)
	// {
	// 	std::cout << ref.chrMap[ref.cpgTable[ref.metaCpGs[hIt->first].start].chrom] << "\t";
	// 	std::cout << ref.cpgTable[ref.metaCpGs[hIt->first].start].pos << "\t";
	// 	std::cout << hIt->first << "\t";
	// 	std::cout << (uint64_t)std::get<0>(hIt->second) << "\t";
	// 	std::cout << (uint64_t)std::get<1>(hIt->second) << "\t";
	// 	std::cout << std::get<2>(hIt->second) << "\n";
	// }
	// return candCount;
	return false;
}






inline bool ReadQueue::extractSingleMatch(std::vector<MATCH::match>& fwdMatches, std::vector<MATCH::match>& revMatches, Read& r, std::string& revSeq)
{
	// Construct artificial best match
	MATCH::match bestMat = MATCH::constructMatch(0, MyConst::MISCOUNT + MyConst::ADDMIS + 1,0,0,0);
	bool isUnique = true;
	// indicates which of the reads had the best match
	// 1 = original read
	// 2 = reverse comp of read
	// 0 = none
	unsigned int matchedReadID = 0;
	// Extract best match from list of found matches of fwd reads
	for (MATCH::match mat : fwdMatches)
	{
		if (MATCH::getErrNum(mat) == MATCH::getErrNum(bestMat))
		{
			// Check if same match
			uint32_t matPos = ref.metaWindows[MATCH::getMetaID(mat)].startPos + MATCH::getOffset(mat);
			uint32_t bestMatPos = ref.metaWindows[MATCH::getMetaID(bestMat)].startPos + MATCH::getOffset(bestMat);

			// Dealing with large offsets, we need unsigned. Hence check both directions.
			// check within offset of one for security reasons (insertions deletions etc)
			if (matPos - bestMatPos > 1 && bestMatPos - matPos > 1)
			{
				isUnique = false;
			}
		} else if (MATCH::getErrNum(mat) < MATCH::getErrNum(bestMat))
		{
			// update bestMatch
			matchedReadID = 1;
			bestMat = mat;
			isUnique = true;
		}

	}
	for (MATCH::match mat : revMatches)
	{
		if (MATCH::getErrNum(mat) == MATCH::getErrNum(bestMat))
		{
			// Check if same match
			uint32_t matPos = ref.metaWindows[MATCH::getMetaID(mat)].startPos + MATCH::getOffset(mat);
			uint32_t bestMatPos = ref.metaWindows[MATCH::getMetaID(bestMat)].startPos + MATCH::getOffset(bestMat);

			// Dealing with large offsets, we need unsigned. Hence check both directions.
			// check within offset of one for security reasons (insertions deletions etc)
			if (matPos - bestMatPos > 1 && bestMatPos - matPos > 1)
			{
				isUnique = false;
			}
		} else if (MATCH::getErrNum(mat) < MATCH::getErrNum(bestMat))
		{
			// update bestMatch
			matchedReadID = 2;
			bestMat = mat;
			isUnique = true;
		}
	}
	// Note that either we MUST have found either a match or a nonunqiue matching
	// because only if the match lists contain at least 1 element we call this function
	// check if found match is unique
	if (isUnique)
	{
		r.mat = bestMat;
		if (matchedReadID == 1)
		{
			computeMethLvl(bestMat, r.seq);

		} else if (matchedReadID == 2) {

			computeMethLvl(bestMat, revSeq);
		} else {
			std::cerr << "You should not reach this code.\n\n";
		}

	} else {

		r.isInvalid = true;
	}
	return isUnique;
}





inline int ReadQueue::extractPairedMatch(MATCH::match& mat1, MATCH::match& mat2)
{

	// check if on same chromosome
	// if (ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(mat1)].start].chrom != ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(mat2)].start].chrom)
	if (ref.metaWindows[MATCH::getMetaID(mat1)].chrom != ref.metaWindows[MATCH::getMetaID(mat2)].chrom)
		return -1;
	// check if in range
	uint32_t mat1Pos = ref.metaWindows[MATCH::getMetaID(mat1)].startPos + MATCH::getOffset(mat1);
	uint32_t mat2Pos = ref.metaWindows[MATCH::getMetaID(mat2)].startPos + MATCH::getOffset(mat2);
	uint32_t matDist;
	if (mat2Pos > mat1Pos)
		matDist = mat2Pos - mat1Pos;
	else
		matDist = mat1Pos - mat2Pos;
	if (matDist <= MyConst::MAXPDIST + MyConst::READLEN)
	{
		return (MATCH::getErrNum(mat1) + MATCH::getErrNum(mat2));
	}
	return -1;
}







inline void ReadQueue::computeMethLvl(MATCH::match& mat, std::string& seq)
{

	// retrieve matched metaId
	bool isFwd = MATCH::isFwd(mat);
	uint32_t metaID = MATCH::getMetaID(mat);
	uint16_t offset = MATCH::getOffset(mat);
	uint8_t errNum = MATCH::getErrNum(mat);
	// if no cpg in window, do not carry out alignment
	if (ref.metaWindows[metaID].startInd == MyConst::CPGDUMMY)
	{
		return;
	}

	// if no errors -> simple lookup of sequences
	if (errNum == 0)
	{

		// retrieve chromosome and position of match
// 		if (isStart)
// 		{
// 			struct metaCpG& m = ref.metaStartCpGs[metaID];
// 			// uint8_t chrom = ref.cpgStartTable[m.start].chrom;
// 			for (uint32_t cpgId = m.start; cpgId <= m.end; ++cpgId)
// 			{
// 				// check if CpG is too far downstream of read match
// 				// i.e. no overlap
// 				if (ref.cpgStartTable[cpgId].pos < offset - seq.size())
// 					continue;
// 				// check if too far upstream
// 				if (ref.cpgStartTable[cpgId].pos + 1 > offset)
// 					break;
//
//
// 				// position of CpG in read
// 				// last term represents position of start of read in reference sequence
// 				uint32_t readCpGPos = ref.cpgStartTable[cpgId].pos - (offset - seq.size() + 1);
// 				if (isFwd)
// 				{
// 					if (seq[readCpGPos] == 'T')
// 					{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 						++methLevelsStart[cpgId].unmethFwd;
// 					} else {
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 						++methLevelsStart[cpgId].methFwd;
// 					}
// 				} else {
//
// 					if (seq[seq.size() - readCpGPos - 1] == 'T')
// 					{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 						++methLevelsStart[cpgId].unmethRev;
// 					} else {
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 						++methLevelsStart[cpgId].methRev;
// 					}
// 				}
// 			}
//
//
// 		// normal match
		// } else {

			// struct metaCpG& m = ref.metaCpGs[metaID];
			metaWindow& m = ref.metaWindows[metaID];
			// uint8_t chrom = ref.cpgTable[m.start].chrom;
			// uint32_t metaPos = ref.cpgTable[m.startInd].pos;
			uint32_t metaPos = m.startPos;
			const uint32_t minPos = metaPos + offset - (seq.size() - 1);
			const uint32_t maxPos = metaPos + offset;
			for (uint32_t cpgId = m.startInd; cpgId <= m.endInd; ++cpgId)
			{
				// check if CpG is too far downstream of read match
				// i.e. no overlap
				if (ref.cpgTable[cpgId].pos + MyConst::READLEN - 2 < minPos)
					continue;
				// check if too far upstream
				if (isFwd)
				{
					if (ref.cpgTable[cpgId].pos + MyConst::READLEN - 2 > maxPos)
						break;
				} else {
					if (ref.cpgTable[cpgId].pos + MyConst::READLEN - 1 > maxPos)
						break;
				}



				// position of CpG in read
				uint32_t readCpGPos = ref.cpgTable[cpgId].pos + MyConst::READLEN - 2 - (metaPos + offset - (seq.size() - 1));
				if (isFwd)
				{
					if (seq[readCpGPos] == 'T')
					{
#ifdef _OPENMP
#pragma omp atomic
#endif
						++methLevels[cpgId].unmethFwd;
// 						if (isSC)
// 						{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 							++methLevelsSc[cpgId].unmethFwd;
// 						}
					}
					else if (seq[readCpGPos] == 'C')
					{
#ifdef _OPENMP
#pragma omp atomic
#endif
						++methLevels[cpgId].methFwd;
// 						if (isSC)
// 						{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 							++methLevelsSc[cpgId].methFwd;
// 						}
					// TODO
					// } else {

// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 						{
// 						std::cerr << "WHAT fwd!?\n";
// 						}
					}
// 					// else std::cout << "This should not happen 1!\n";

				} else {

					if (seq[seq.size() - readCpGPos - 2] == 'T')
					{
#ifdef _OPENMP
#pragma omp atomic
#endif
						++methLevels[cpgId].unmethRev;
// 						if (isSC)
// 						{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 							++methLevelsSc[cpgId].unmethRev;
// 						}
					}
					else  if (seq[seq.size() - readCpGPos - 2] == 'C')
					{
#ifdef _OPENMP
#pragma omp atomic
#endif
						++methLevels[cpgId].methRev;
// 						if (isSC)
// 						{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 							++methLevelsSc[cpgId].methRev;
// 						}
					// TODO
// 					} else {
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 						{
// 						std::cerr << "WHAT rev!?\n";
// 						}
					}
				}
			}
		// }

	// if match was errornous
	} else {

// 		if (isStart)
// 		{
// 			struct metaCpG& m = ref.metaStartCpGs[metaID];
// 			uint8_t chrom = ref.cpgStartTable[m.start].chrom;
//
// 			const char* refSeq = ref.fullSeq[chrom].data() + offset;
//
// 			// init levenshtein DP algo
// 			LevenshtDP<uint16_t, MyConst::MISCOUNT + MyConst::ADDMIS> lev(seq, refSeq);
// 			std::vector<ERROR_T> alignment;
//
// 			// minimum position for overlap
// 			uint32_t minPos = offset - (seq.size() - 1);
// 			uint32_t maxPos = offset - 1;
// 			// positions of first/ last overlapping CpG
// 			uint32_t minIndex = m.start;
// 			uint32_t maxIndex = m.end;
// 			for (uint32_t cpgID = m.start; cpgID <= m.end; ++cpgID)
// 			{
// 				if (ref.cpgStartTable[cpgID].pos < minPos)
// 				{
// 					++minIndex;
// 				} else if (ref.cpgStartTable[cpgID].pos > maxPos)
// 				{
// 					maxIndex = cpgID - 1;
// 					break;
// 				}
// 			}
//
// 			if (isFwd)
// 			{
//
// 				// compute alignment
// 				lev.runDPFill<CompiFwd>(cmpFwd);
// 				lev.backtrackDP<CompiFwd>(cmpFwd, alignment);
// 				uint32_t refSeqPos = offset;
// 				uint32_t readSeqPos = seq.size() - 1;
// 				uint32_t alignPos = alignment.size() - 1;
//
// 				// go through all overlapping CpGs from back,
// 				// move through read and reference according to alignment
// 				// if position of CpG is hit, compare and count
// 				for (uint32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
// 				{
// 					// align until this CpG
// 					while (ref.cpgStartTable[cpgID].pos < refSeqPos)
// 					{
// 						switch (alignment[alignPos--])
// 						{
// 							case (MATCHING):
// 							case (MISMATCH):
// 								--readSeqPos;
// 								--refSeqPos;
// 								break;
// 							case (DELETION):
// 								--refSeqPos;
// 								break;
// 							case(INSERTION):
// 								--readSeqPos;
// 								break;
// 						}
// 						// check if we have a CpG aligned to the reference CpG
// 						if (seq[readSeqPos + 1] == 'G')
// 						{
// 							// check for methylated C
// 							if (seq[readSeqPos] == 'C')
// 							{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 								++methLevelsStart[cpgID].methFwd;
// 							}
// 							else if (seq[readSeqPos] == 'T')
// 							{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 								++methLevelsStart[cpgID].unmethFwd;
// 							}
//
// 						}
// 					}
// 				}
//
// 			} else {
//
// 				lev.runDPFillRev<CompiRev>(cmpRev);
// 				lev.backtrackDPRev<CompiRev>(cmpRev, alignment);
// 				uint32_t refSeqPos = offset;
// 				uint32_t readSeqPos = seq.size() - 1;
// 				uint32_t alignPos = alignment.size() - 1;
//
// 				// go through all overlapping CpGs from back,
// 				// move through read and reference according to alignment
// 				// if position of CpG is hit, compare and count
// 				for (uint32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
// 				{
// 					// align until this CpG
// 					while (ref.cpgStartTable[cpgID].pos < refSeqPos)
// 					{
// 						switch (alignment[alignPos--])
// 						{
// 							case (MATCHING):
// 							case (MISMATCH):
// 								--readSeqPos;
// 								--refSeqPos;
// 								break;
// 							case (DELETION):
// 								--refSeqPos;
// 								break;
// 							case(INSERTION):
// 								--readSeqPos;
// 								break;
// 						}
// 						// check if we have a CpG aligned to the reference CpG
// 						if (seq[readSeqPos] == 'G')
// 						{
// 							// check for methylated C (on reverse!)
// 							if (seq[readSeqPos + 1] == 'C')
// 							{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 								++methLevelsStart[cpgID].methRev;
// 							}
// 							else if (seq[readSeqPos] == 'T')
// 							{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 								++methLevelsStart[cpgID].unmethRev;
// 							}
//
// 						}
// 					}
// 				}
// 			}
//
//
// 		// not start CpG (normal)
// 		} else {

			// struct metaCpG& m = ref.metaCpGs[metaID];
			metaWindow& m = ref.metaWindows[metaID];
			uint8_t chrom = ref.cpgTable[m.startInd].chrom;
			// uint32_t metaPos = ref.cpgTable[m.startInd].pos;
			uint32_t metaPos = m.startPos;

			const char* refSeq = ref.fullSeq[chrom].data() + metaPos + offset;

			// init levenshtein DP algo
			LevenshtDP<uint16_t, MyConst::MISCOUNT + MyConst::ADDMIS> lev(seq, refSeq);
			std::vector<ERROR_T> alignment;

			// minimum position for overlap
			const uint32_t minPos = metaPos + offset - (seq.size() - 1);
			const uint32_t maxPos = metaPos + offset;
			// positions of first/ last overlapping CpG
			int32_t minIndex = m.startInd;
			int32_t maxIndex = m.endInd;
			for (int32_t cpgID = m.startInd; cpgID <= m.endInd; ++cpgID)
			{
				if (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 < minPos)
				{
					++minIndex;

				} else if (isFwd)
				{
					if (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 > maxPos)
					{
						maxIndex = cpgID - 1;
						break;
					}
				} else {
					if (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 > maxPos)
					{
						maxIndex = cpgID - 1;
						break;
					}
				}
			}

			if (isFwd)
			{

				// compute alignment
				lev.runDPFill<CompiFwd>(cmpFwd);
				lev.backtrackDP<CompiFwd>(cmpFwd, alignment);
				// current position in read and reference
				uint32_t refSeqPos = metaPos + offset;
				int32_t readSeqPos = seq.size() - 1;
				// current position in alignment (note that we align from right to left)
				int32_t alignPos = alignment.size() - 1;

				// go through all overlapping CpGs from back,
				// move through read and reference according to alignment
				// if position of CpG is hit, compare and count
				// TODO
				// std::string al (alignment.size(),'%');
				// bool hasShift = false;
				for (int32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
				{
					// align until this CpG
					while (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 < refSeqPos && alignPos >= 0)
					{
						switch (alignment[alignPos])
						{
							case (MATCHING):
								// TODO delete this
								// al[alignPos] = '=';
								// --readSeqPos;
								// --refSeqPos;
								// break;
							case (MISMATCH):
								// TODO
								// al[alignPos] = '!';
								--readSeqPos;
								--refSeqPos;
								break;
							case (DELETION):
								// al[alignPos] = '-';
								--refSeqPos;
								break;
							case(INSERTION):
								// al[alignPos] = '+';
								--readSeqPos;
								break;
						}
						if (readSeqPos < 0)
							break;
						if (readSeqPos == seq.size() - 1)
							continue;
						--alignPos;
					}
					if (readSeqPos < 0)
					{
						break;
					}
					// check if we have a CpG aligned to the reference CpG
					// if (seq[readSeqPos + 1] == 'G')
					// {
						// check for unmethylated C
						if (seq[readSeqPos] == 'C')
						{
#ifdef _OPENMP
#pragma omp atomic
#endif
							++methLevels[cpgID].methFwd;
// 							if (isSC)
// 							{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 								++methLevelsSc[cpgID].methFwd;
// 							}
						}
						else if (seq[readSeqPos] == 'T')
						{
#ifdef _OPENMP
#pragma omp atomic
#endif
							++methLevels[cpgID].unmethFwd;
// 							if (isSC)
// 							{
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
// 								++methLevelsSc[cpgID].unmethFwd;
// 							}
						// TODO
// 						} else {
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 							{
// 							// hasShift = true;
// 							std::cerr << "WHAT on fwd err!?\t" << seq[readSeqPos] << "\n";
// 							std::cerr << "\tPos: " << readSeqPos << "\n";
// 							std::cerr << "Sequence original/read:\n\t" << std::string(ref.fullSeq[chrom].data() + metaPos + offset - 100, 100) << "\n\t" << seq << "\n";
// 							}
						}

					// }
				}
				// TODO
				// if (hasShift)
				// 	std::cerr << "\t" << al << "\n\n";

			} else {

				lev.runDPFillRev<CompiRev>(cmpRev);
				lev.backtrackDPRev<CompiRev>(cmpRev, alignment);
				uint32_t refSeqPos = metaPos + offset;
				int32_t readSeqPos = seq.size() - 1;
				int32_t alignPos = alignment.size() - 1;
				std::reverse(seq.begin(),seq.end());

				// go through all overlapping CpGs from back,
				// move through read and reference according to alignment
				// if position of CpG is hit, compare and count
				// TODO
				// std::string al (alignment.size(),'%');
				// bool hasShift = false;
				for (int32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
				{
					// align until this CpG
					while (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 < refSeqPos && alignPos >= 0)
					{
						switch (alignment[alignPos])
						{
							case (MATCHING):
								// TODO delete this
								// al[alignPos] = '=';
								// --readSeqPos;
								// --refSeqPos;
								// break;
							case (MISMATCH):
								// TODO
								// al[alignPos] = '!';
								--readSeqPos;
								--refSeqPos;
								break;
							case (DELETION):
								// al[alignPos] = '-';
								--refSeqPos;
								break;
							case(INSERTION):
								// al[alignPos] = '+';
								--readSeqPos;
								break;
						}
						if (readSeqPos < 0)
							break;
						if (readSeqPos == seq.size() - 1)
							continue;
						--alignPos;
					}
					// check if we have a CpG aligned to the reference CpG
					// if (seq[readSeqPos] == 'G')
					// {
						// check for unmethylated C
						if (seq[readSeqPos + 1] == 'C')
						{
#ifdef _OPENMP
#pragma omp atomic
#endif
							++methLevels[cpgID].methRev;
							if (isSC)
							{
#ifdef _OPENMP
#pragma omp atomic
#endif
								++methLevelsSc[cpgID].methRev;
							}
						}
						else if (seq[readSeqPos + 1] == 'T')
						{
#ifdef _OPENMP
#pragma omp atomic
#endif
							++methLevels[cpgID].unmethRev;
							if (isSC)
							{
#ifdef _OPENMP
#pragma omp atomic
#endif
								++methLevelsSc[cpgID].unmethRev;
							}
						// TODO
// 						} else {
//
// 							break;
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// 							{
// 							// hasShift = true;
// 							std::cerr << "WHAT on rev err!?\t" << seq[readSeqPos + 1] << "\n";
// 							std::cerr << "\tPos: " << readSeqPos + 1 << "\n";
// 							std::cerr << "Sequence original/read:\n\t" << std::string(ref.fullSeq[chrom].da
//                                                     std::cerr << "Sequence original/read:\n\t" << std::string(ref.fullSeq[chrom].data() + metaPos + offset - 100, 100) << "\n\t" << seq << "\n\n";
                               // TODO
                               // if (hasShift)
                               //      std::cerr << "\t" << al << "\n\n";
						}
					// }
				}
			}
		// }
	}
}


inline bool ReadQueue::matchFwdFirst(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa)
{
	constexpr int contextWLen = (int)(((double)MyConst::MAXPDIST / MyConst::WINLEN)) + 1;
	auto& fwdMetaIDs_t = paired_fwdMetaIDs[omp_get_thread_num()];
	if (std::get<0>(meta.second) < std::max(bmCount - (int32_t)MyConst::KMERDIST, (int32_t)qThreshold))
		return false;
	if (!std::get<2>(meta.second))
		return true;
	auto& mIt = fwdMetaIDs_t[meta.first];

	auto startIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin() + ref.metaWindows[meta.first].startPos;
	auto endIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin() + ref.metaWindows[meta.first].startPos + MyConst::WINLEN + MyConst::MISCOUNT + MyConst::ADDMIS - 1;
	if (endIt > ref.fullSeq[ref.metaWindows[meta.first].chrom].end())
	{
		endIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].end();
	}

	// use shift and to find all matchings
	std::vector<uint64_t> matchings;
	std::vector<uint8_t> errors;

	sa.querySeq(startIt, endIt, matchings, errors);

	size_t i = 0;
	// compare first found match with last found match of previous meta CpG
	if (matchings.size() > 0)
	{
		// compare chromosome and offset
		if (matchings[0] + ref.metaWindows[meta.first].startPos == prevOff && ref.metaWindows[meta.first].chrom == prevChr)
		{
			++i;
		}
	}
	// translate found matchings
	for (; i < matchings.size(); ++i)
	{
		mats.push_back(MATCH::constructMatch(matchings[i], errors[i], 1, 0, meta.first));
	}
	if (matchings.size() > 0)
	{
		bmCount = std::max(bmCount, (int32_t)std::get<0>(meta.second));

		// update boolean for this and adjacent windows
		std::get<3>(mIt) = true;
		for (int cOff = 1; cOff <= contextWLen; ++cOff)
		{
			auto it = fwdMetaIDs_t.find(meta.first + cOff);
			if (it != fwdMetaIDs_t.end())
				std::get<3>(it.value()) = true;
			it = fwdMetaIDs_t.find(meta.first - cOff);
			if (it != fwdMetaIDs_t.end())
				std::get<3>(it.value()) = true;
		}
		prevChr = ref.metaWindows[meta.first].chrom;
		prevOff = ref.metaWindows[meta.first].startPos + matchings[matchings.size() - 1];

	} else {

		prevChr = 0;
		prevOff = 0xffffffffffffffffULL;
	}
	return true;
}




inline bool ReadQueue::matchRevFirst(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa)
{
	constexpr int contextWLen = (int)(((double)MyConst::MAXPDIST / MyConst::WINLEN)) + 1;
	auto& revMetaIDs_t = paired_revMetaIDs[omp_get_thread_num()];

	if (std::get<0>(meta.second) < std::max(bmCount - (int32_t)MyConst::KMERDIST, (int32_t)qThreshold))
		return false;
	if (!std::get<2>(meta.second))
		return true;

	auto& mIt = revMetaIDs_t[meta.first];

	// retrieve sequence
	auto endIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin();

	if (ref.metaWindows[meta.first].startPos > 0)
	{
		endIt += ref.metaWindows[meta.first].startPos - 1;
	}

	auto startIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin() + ref.metaWindows[meta.first].startPos + MyConst::WINLEN + MyConst::MISCOUNT + MyConst::ADDMIS - 1;
	if (startIt >= ref.fullSeq[ref.metaWindows[meta.first].chrom].end())
	{
		startIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].end() - 1;
	}

	// use shift and to find all matchings
	std::vector<uint64_t> matchings;
	std::vector<uint8_t> errors;

	sa.queryRevSeq(startIt, endIt, matchings, errors);

	size_t i = 0;
	// compare first found match with last found match of previous meta CpG
	if (matchings.size() > 0)
	{
		// compare chromosome and offset
		if (matchings[0] + ref.metaWindows[meta.first].startPos == prevOff && ref.metaWindows[meta.first].chrom == prevChr)
		{
			++i;
		}
	}
	// translate found matchings
	for (; i < matchings.size(); ++i)
	{
		mats.push_back(MATCH::constructMatch(matchings[i], errors[i], 0, 0, meta.first));
	}
	if (matchings.size() > 0)
	{

		bmCount = std::max(bmCount, (int32_t)std::get<0>(meta.second));
		std::get<3>(mIt) = true;
		for (int cOff = 1; cOff <= contextWLen; ++cOff)
		{
			auto it = revMetaIDs_t.find(meta.first + cOff);
			if (it != revMetaIDs_t.end())
				std::get<3>(it.value()) = true;
			it = revMetaIDs_t.find(meta.first - cOff);
			if (it != revMetaIDs_t.end())
				std::get<3>(it.value()) = true;
		}
		prevChr = ref.metaWindows[meta.first].chrom;
		prevOff = ref.metaWindows[meta.first].startPos + matchings[matchings.size() - 1];

	} else {

		prevChr = 0;
		prevOff = 0xffffffffffffffffULL;
	}
	return true;
}


inline bool ReadQueue::matchFwdSecond(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa)
{
	if (std::get<1>(meta.second) < std::max(bmCount - (int32_t)MyConst::KMERDIST, (int32_t)qThreshold))
		return false;
	if (!std::get<3>(meta.second))
		return true;

	auto startIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin() + ref.metaWindows[meta.first].startPos;
	auto endIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin() + ref.metaWindows[meta.first].startPos + MyConst::WINLEN + MyConst::MISCOUNT + MyConst::ADDMIS - 1;
	if (endIt > ref.fullSeq[ref.metaWindows[meta.first].chrom].end())
	{
		endIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].end();
	}

	// use shift and to find all matchings
	std::vector<uint64_t> matchings;
	std::vector<uint8_t> errors;
	sa.querySeq(startIt, endIt, matchings, errors);

	size_t i = 0;
	// compare first found match with last found match of previous meta CpG
	if (matchings.size() > 0)
	{
		// compare chromosome and offset
		if (matchings[0] + ref.metaWindows[meta.first].startPos == prevOff && ref.metaWindows[meta.first].chrom == prevChr)
		{
			++i;
		}
	}
	// translate found matchings
	for (; i < matchings.size(); ++i)
	{
		mats.push_back(MATCH::constructMatch(matchings[i], errors[i], 1, 0, meta.first));
	}
	if (matchings.size() > 0)
	{

		bmCount = std::max(bmCount, (int32_t)std::get<1>(meta.second));
		prevChr = ref.metaWindows[meta.first].chrom;
		prevOff = ref.metaWindows[meta.first].startPos + matchings[matchings.size() - 1];

	} else {

		prevChr = 0;
		prevOff = 0xffffffffffffffffULL;
	}
	return true;
}


inline bool ReadQueue::matchRevSecond(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa)
{
	if (std::get<1>(meta.second) < std::max(bmCount - (int32_t)MyConst::KMERDIST, (int32_t)qThreshold))
		return false;
	if (!std::get<3>(meta.second))
		return true;

	// retrieve sequence
	auto endIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin() + ref.metaWindows[meta.first].startPos - 1;
	auto startIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].begin() + ref.metaWindows[meta.first].startPos + MyConst::WINLEN + MyConst::MISCOUNT + MyConst::ADDMIS - 1;
	if (startIt >= ref.fullSeq[ref.metaWindows[meta.first].chrom].end())
	{
		startIt = ref.fullSeq[ref.metaWindows[meta.first].chrom].end() - 1;
	}

	// use shift and to find all matchings
	std::vector<uint64_t> matchings;
	std::vector<uint8_t> errors;
	sa.queryRevSeq(startIt, endIt, matchings, errors);

	size_t i = 0;
	// compare first found match with last found match of previous meta CpG
	if (matchings.size() > 0)
	{
		// compare chromosome and offset
		if (matchings[0] + ref.metaWindows[meta.first].startPos == prevOff && ref.metaWindows[meta.first].chrom == prevChr)
		{
			++i;
		}
	}
	// translate found matchings
	for (; i < matchings.size(); ++i)
	{
		mats.push_back(MATCH::constructMatch(matchings[i], errors[i], 0, 0, meta.first));
	}
	if (matchings.size() > 0)
	{

		bmCount = std::max(bmCount, (int32_t)std::get<1>(meta.second));
		prevChr = ref.metaWindows[meta.first].chrom;
		prevOff = ref.metaWindows[meta.first].startPos + matchings[matchings.size() - 1];

	} else {

		prevChr = 0;
		prevOff = 0xffffffffffffffffULL;
	}
	return true;
}
