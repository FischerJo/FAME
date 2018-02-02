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

#include "ReadQueue.h"

ReadQueue::ReadQueue(const char* filePath, RefGenome& reference, bool isGZ) :
        ref(reference)
    ,   readBuffer(MyConst::CHUNKSIZE)
    ,   methLevels(ref.cpgTable.size())
    ,   methLevelsStart(ref.cpgStartTable.size())
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

        // fwdMetaIDs[i] = std::unordered_map<uint32_t, uint16_t, MetaHash>();
        // revMetaIDs[i] = std::unordered_map<uint32_t, uint16_t, MetaHash>();
        // fwdMetaIDs[i] = spp::sparse_hash_map<uint32_t, uint16_t, MetaHash>();
        // revMetaIDs[i] = spp::sparse_hash_map<uint32_t, uint16_t, MetaHash>();
        fwdMetaIDs[i] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        revMetaIDs[i] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        fwdMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 1);
        revMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 1);
        fwdMetaIDs[i].set_empty_key(ref.metaCpGs.size() + 2);
        revMetaIDs[i].set_empty_key(ref.metaCpGs.size() + 2);
        countsFwd[i] = std::vector<uint16_t>();
        countsRev[i] = std::vector<uint16_t>();
        countsFwdStart[i] = std::vector<uint16_t>();
        countsRevStart[i] = std::vector<uint16_t>();
    }
    // fill array mapping - locale specific filling
    lmap['A'] = 0;
    lmap['C'] = 1;
    lmap['G'] = 2;
    lmap['T'] = 3;

    // TODO
    // write to statFile
    // for (size_t i = 0; i < MyConst::HTABSIZE; ++i)
    // {
        // if (ref.tabIndex[i+1] - ref.tabIndex[i] > 2000000)
        // {
        //
        //     auto startIt = ref.kmerTable.begin() + ref.tabIndex[i];
        //     auto endIt = ref.kmerTable.begin() + ref.tabIndex[i + 1];
        //     auto tit = ref.strandTable.begin() + ref.tabIndex[i];
        //     for (auto it = startIt; it != endIt; ++it, ++tit)
        //     {
        //         KMER::kmer& k = *it;
        //         const uint64_t m = KMER::getMetaCpG(k);
        //         const bool isStart = KMER::isStartCpG(k);
        //         if (!isStart)
        //         {
        //             const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
        //             uint64_t offset = KMER::getOffset(k);
        //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + offset + startCpg.pos;
        //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + offset + MyConst::KMERLEN + startCpg.pos;
        //             if (*tit)
        //             {
        //                 statFile << std::string(stIt, enIt) << "\n";
        //             } else {
        //
        //                 statFile << "REV---" << std::string(stIt, enIt) << "\n";
        //             }
        //         }
        //     }
        // }
    //     statFile << ref.tabIndex[i+1] - ref.tabIndex[i] << "\n";
    // }
    // std::vector<uint32_t> counts (10000000, 0);
    // for (size_t i = 0; i < MyConst::HTABSIZE; ++i)
    //     ++counts[ref.tabIndex[i+1] - ref.tabIndex[i]];
    // for (size_t j = 0; j < counts.size(); ++j)
    //     statFile << j << "\t" << counts[j] << "\n";
    // statFile.close();
}
ReadQueue::ReadQueue(const char* filePath, const char* filePath2, RefGenome& reference, bool isGZ) :
        ref(reference)
    ,   readBuffer(MyConst::CHUNKSIZE)
    ,   readBuffer2(MyConst::CHUNKSIZE)
    ,   methLevels(ref.cpgTable.size())
    ,   methLevelsStart(ref.cpgStartTable.size())
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

        // fwdMetaIDs[i] = std::unordered_map<uint32_t, uint16_t, MetaHash>();
        // revMetaIDs[i] = std::unordered_map<uint32_t, uint16_t, MetaHash>();
        // fwdMetaIDs[i] = spp::sparse_hash_map<uint32_t, uint16_t, MetaHash>();
        // revMetaIDs[i] = spp::sparse_hash_map<uint32_t, uint16_t, MetaHash>();
        fwdMetaIDs[i] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        revMetaIDs[i] = google::dense_hash_map<uint32_t, uint16_t, MetaHash>();
        fwdMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 1);
        revMetaIDs[i].set_deleted_key(ref.metaCpGs.size() + 1);
        fwdMetaIDs[i].set_empty_key(ref.metaCpGs.size() + 2);
        revMetaIDs[i].set_empty_key(ref.metaCpGs.size() + 2);
        countsFwd[i] = std::vector<uint16_t>();
        countsRev[i] = std::vector<uint16_t>();
        countsFwdStart[i] = std::vector<uint16_t>();
        countsRevStart[i] = std::vector<uint16_t>();
    }
    // fill array mapping - locale specific filling
    lmap['A'] = 0;
    lmap['C'] = 1;
    lmap['G'] = 2;
    lmap['T'] = 3;

    // TODO
    // write to statFile
    // for (size_t i = 0; i < MyConst::HTABSIZE; ++i)
    // {
        // if (ref.tabIndex[i+1] - ref.tabIndex[i] > 2000000)
        // {
        //
        //     auto startIt = ref.kmerTable.begin() + ref.tabIndex[i];
        //     auto endIt = ref.kmerTable.begin() + ref.tabIndex[i + 1];
        //     auto tit = ref.strandTable.begin() + ref.tabIndex[i];
        //     for (auto it = startIt; it != endIt; ++it, ++tit)
        //     {
        //         KMER::kmer& k = *it;
        //         const uint64_t m = KMER::getMetaCpG(k);
        //         const bool isStart = KMER::isStartCpG(k);
        //         if (!isStart)
        //         {
        //             const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
        //             uint64_t offset = KMER::getOffset(k);
        //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + offset + startCpg.pos;
        //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + offset + MyConst::KMERLEN + startCpg.pos;
        //             if (*tit)
        //             {
        //                 statFile << std::string(stIt, enIt) << "\n";
        //             } else {
        //
        //                 statFile << "REV---" << std::string(stIt, enIt) << "\n";
        //             }
        //         }
        //     }
        // }
    //     statFile << ref.tabIndex[i+1] - ref.tabIndex[i] << "\n";
    // }
    // std::vector<uint32_t> counts (10000000, 0);
    // for (size_t i = 0; i < MyConst::HTABSIZE; ++i)
    //     ++counts[ref.tabIndex[i+1] - ref.tabIndex[i]];
    // for (size_t j = 0; j < counts.size(); ++j)
    //     statFile << j << "\t" << counts[j] << "\n";
    // statFile.close();
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
        readBuffer[readCounter] = std::move(Read(seq, id));
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
            readBuffer2[readCounter2] = std::move(Read(seq, id));
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
        readBuffer[readCounter] = std::move(Read(seq, id));
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
            readBuffer2[readCounter2] = std::move(Read(seq, id));
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

bool ReadQueue::matchReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch)
{

    // TODO
    unsigned int pCount = 0;

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


        if (readSize < 90)
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

        // set qgram threshold
        uint16_t qThreshold = readSize - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT) + 1;
        if (qThreshold > readSize)
            qThreshold = 0;
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
        ShiftAnd<MyConst::MISCOUNT> saFwd(r.seq, lmap);
        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // of << runtime << "\t";
        // startTime = std::chrono::high_resolution_clock::now();
        int succQueryFwd = saQuerySeedSetRef(saFwd, matchFwd, qThreshold);
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
        ShiftAnd<MyConst::MISCOUNT> saRev(revSeq, lmap);
        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        // of << runtime << "\n";
        // startTime = std::chrono::high_resolution_clock::now();
        int succQueryRev = saQuerySeedSetRef(saRev, matchRev, qThreshold);
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

        // found match for fwd and rev strand
        if (succQueryFwd == 1 && succQueryRev == 1)
        {

            uint8_t fwdErr = MATCH::getErrNum(matchFwd);
            uint8_t revErr = MATCH::getErrNum(matchRev);
            // of << "Found match with read and reverse complement. Errors (fwd/rev): " << fwdErr << "/" << revErr;
            // of << "\nMatching strands are (fwd/rev): " << MATCH::isFwd(matchFwd) << "/" << MATCH::isFwd(matchRev) << "\n";

            // check which one has fewer errors
            if (fwdErr < revErr)
            {

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

                        m1_pos = ref.cpgTable[ref.metaCpGs[metaFwd].start].pos + offFwd;
                        m2_pos = ref.cpgTable[ref.metaCpGs[metaRev].start].pos + offRev;
                    }
                    // test if same match in same region
                    if ((m1_isStart == m2_isStart) && (m1_isFwd == m2_isFwd) && (m1_pos == m2_pos))
                    {
                        ++succMatchT;
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
//                 uint64_t hVal = ntHash::NTP64(r.seq.data()) % MyConst::HTABSIZE;
//                 auto startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//                 auto endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//                 auto tit = ref.strandTable.begin() + ref.tabIndex[hVal];
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
//                 hVal = ntHash::NTP64(r.seq.data()+70) % MyConst::HTABSIZE;
//                 startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//                 endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//                 tit = ref.strandTable.begin() + ref.tabIndex[hVal];
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
//                 hVal = ntHash::NTP64(revSeq.data()) % MyConst::HTABSIZE;
//                 startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//                 endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//                 tit = ref.strandTable.begin() + ref.tabIndex[hVal];
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
//                 hVal = ntHash::NTP64(revSeq.data()+70) % MyConst::HTABSIZE;
//                 startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
//                 endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
//                 tit = ref.strandTable.begin() + ref.tabIndex[hVal];
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
// // END PRAGMA OMP CRITICAL
// }
            }
            // if (unSuccMatch > 10)
            // {
            //     of.close();
            //     exit(1);
            // }
        }
// TODO
//         if (!r.isInvalid)
//         {
//             if (ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(r.mat)].start].chrom != 21)
//             {
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
//                 ++wrongChrCount;
//             }
//         }
//         // validate
//         if (!r.isInvalid)
//         {
//             size_t lPos = r.id.find_last_of('_');
//             std::string stringOffset (r.id.begin() + lPos + 1, r.id.end());
//             size_t rPos = r.id.find_last_of('R');
//             std::string stringChrom (r.id.begin() + 1 + rPos, r.id.begin() + lPos);
//             uint8_t chrom = std::stoul(stringChrom);
//             unsigned long offset = std::stoul(stringOffset);
//             struct CpG& cpg = ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(r.mat)].start];
//             uint64_t matchedOff = MATCH::getOffset(r.mat) + cpg.pos;
//             // if (!MATCH::isFwd(r.mat))
//             // {
//             //     if (MATCH::getErrNum(r.mat) > 0)
//             //     {
//             //         of << "This should be a non-unique match. Bad.\n";
//             //         of << "Number of errors: " << MATCH::getErrNum(r.mat) << "\n";
//             //         of << "Offset should be " << offset << " but is " << matchedOff << "\n";
//             //         continue;
//             //     }
//             // }
//             if (matchedOff > 100 && (matchedOff > offset + 102 || matchedOff < offset + 98))
//             {
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
//                 ++pCount;
//                 // std::cout << "Wrong match in read " << std::string(r.id.begin(), r.id.begin() + r.id.find_first_of('_')) << ". Offset should be " << offset << " in Chromosome " << static_cast<uint16_t>(chrom) << " but is " << matchedOff << " in chromosome " << static_cast<uint16_t>(cpg.chrom) << "\n";
//                 // std::cout << r.id << "\n";
//                 // std::cout << "Number of errors: " << MATCH::getErrNum(r.mat) << "\n";
//                 // if (succQueryRev == 1)
//                 //     std::cout << "It's a reverse match...\n";
//                 // std::cout << "Is forward: " << MATCH::isFwd(r.mat) << "\n";
//                 // std::cout << "Read sequence:\n" << r.seq << "\n\n";
//                 // std::cout << "Reference sequence (match):\n" << std::string(ref.fullSeq[cpg.chrom].begin() + matchedOff - 99, ref.fullSeq[cpg.chrom].begin() + matchedOff + 1) << "\n\n";
//                 // std::cout << "Reference sequence (source):\n" << std::string(ref.fullSeq[chrom].begin() + offset, ref.fullSeq[chrom].begin() + offset + 100) << "\n\n\n";
//             }
//             // if (pCount > 5)
//             // {
//             //
//             //     of.close();
//             //     exit(1);
//             // }
//         }
//
//         // BITMASK COMPARISON ON FULL ALPHABET FOR REMAINING SEEDS
//         // bitMatching(r, fwdSeedsK, fwdSeedsS);
//         // bitMatchingRev(r, revSeedsK, revSeedsS);
//
    }
//     std::cout << pCount << "\n";

    // sum up counts
    for (unsigned int i = 0; i < CORENUM; ++i)
    {
        succMatch += matchStats[i];
        nonUniqueMatch += nonUniqueStats[i];
        unSuccMatch += noMatchStats[i];
    }
    // of.close();
    return true;
}

bool ReadQueue::matchPairedReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch, uint64_t& succPairedMatch)
{


    // reset all counters
    for (unsigned int i = 0; i < CORENUM; ++i)
    {
        matchStats[i] = 0;
        nonUniqueStats[i] = 0;
        noMatchStats[i] = 0;
        matchPairedStats[i] = 0;
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
        uint64_t& succPairedMatchT = matchPairedStats[threadnum];
        Read& r1 = readBuffer[i];
        Read& r2 = readBuffer2[i];

        const size_t readSize1 = r1.seq.size();
        const size_t readSize2 = r2.seq.size();


        if (readSize1 < 90)
        {

            r1.isInvalid = true;
        }
        if (readSize2 < 90)
        {

            r2.isInvalid = true;
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
            unSuccMatchT += 2;
            continue;
        }

        // Stores the found matches for first and second read, resp.
        std::vector<MATCH::match> matches1Fwd;
        std::vector<MATCH::match> matches1Rev;
        matches1Fwd.reserve(20);
        matches1Rev.reserve(20);
        std::vector<MATCH::match> matches2Fwd;
        std::vector<MATCH::match> matches2Rev;
        matches2Fwd.reserve(20);
        matches2Rev.reserve(20);
        // set qgram threshold
        uint16_t qThreshold = readSize1 - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT) + 1;
        if (qThreshold > readSize1)
            qThreshold = 0;

    // MATCH FIRST READ

        getSeedRefs(r1.seq, readSize1, qThreshold);
        ShiftAnd<MyConst::MISCOUNT> saFwd(r1.seq, lmap);
        saQuerySeedSetRefPaired(saFwd, matches1Fwd, qThreshold);

        getSeedRefs(revSeq1, readSize1, qThreshold);
        ShiftAnd<MyConst::MISCOUNT> saRev(revSeq1, lmap);
        saQuerySeedSetRefPaired(saRev, matches1Rev, qThreshold);

        if (matches1Fwd.size() == 0 && matches1Rev.size() == 0)
        {
            r1.isInvalid = true;
            unSuccMatchT += 2;
            continue;
        }


    // MATCH SECOND READ

        getSeedRefs(r2.seq, readSize1, qThreshold);
        ShiftAnd<MyConst::MISCOUNT> saFwd(r2.seq, lmap);
        saQuerySeedSetRefPaired(saFwd, matches2Fwd, qThreshold);

        getSeedRefs(revSeq2, readSize1, qThreshold);
        ShiftAnd<MyConst::MISCOUNT> saRev(revSeq2, lmap);
        saQuerySeedSetRefPaired(saRev, matches2Rev, qThreshold);

        if (matches2Fwd.size() == 0 && matches2Rev.size() == 0)
        {
            r2.isInvalid = true;
            unSuccMatchT += 2;
            continue;
        }


        // TODO:
        // make timer for this

        // current best matching pair (sum of errors)
        int bestErrNum = 2*MyConst::MISCOUNT + 1;
        MATCH::match bestMatch1;
        MATCH::match bestMatch2;
        bool nonUniqueFlag = false;

        for (MATCH::match& mat1 : matches1Fwd)
        {
            for (MATCH::match& mat2Fwd : matches2Fwd)
            {
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
                        nonUniqueFlag = false;

                    }
                }
            }
            for (MATCH::match& mat2Rev : matches2Rev)
            {
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
        for (MATCH::match& mat1 : matches1Rev)
        {
            for (MATCH::match& mat2Fwd : matches2Fwd)
            {
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
                        nonUniqueFlag = false;

                    }
                }
            }
            for (MATCH::match& mat2Rev : matches2Rev)
            {
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


        // construct hash and look up the hash table entries
        // size_t lPos = r1.id.find_last_of('_');
        // std::string stringOffset (r1.id.begin() + lPos + 1, r1.id.end());
        // size_t rPos = r1.id.find_last_of('R');
        // std::string stringChrom (r1.id.begin() + 1 + rPos, r1.id.begin() + lPos);
        // uint8_t chrom = std::stoul(stringChrom);
        // unsigned long offset = std::stoul(stringOffset);
        // of << "\nreal seq/real revSeq/sequence in genome: " << r1.id << "\n" << r1.seq << "\n" << revSeq1 << "\n" << std::string(ref.fullSeq[chrom].begin() + offset, ref.fullSeq[chrom].begin() + offset + 100) << "\n\n\n";
        // of << "\nMeta CpG 14031:" << "\n" << std::string(ref.fullSeq[0].begin() + ref.cpgTable[ref.metaCpGs[14031].start].pos, ref.fullSeq[0].begin() +  ref.cpgTable[ref.metaCpGs[14031].start].pos + 150) << "\n\n\n";
        // of << "\nMeta CpG 14031 start at " << ref.cpgTable[ref.metaCpGs[14031].start].pos << "\n\n\n";
        //
        //
        // uint64_t hVal = ntHash::NTP64(r1.seq.data()) % MyConst::HTABSIZE;
        // auto startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
        // auto endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
        // auto tit = ref.strandTable.begin() + ref.tabIndex[hVal];
        // for (auto it = startIt; it != endIt; ++it, ++tit)
        // {
        //     KMER_S::kmer& k = *it;
        //     const uint32_t m = KMER_S::getMetaCpG(k);
        //     const bool isStart = KMER_S::isStartCpG(k);
        //     if (!isStart)
        //     {
        //         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
        //         if (*tit)
        //         {
        //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
        //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
        //             of << std::string(stIt, enIt) << "\n";
        //         }
        //     }
        // }
        // of << "Last Sequence part:\n";
        // hVal = ntHash::NTP64(r1.seq.data()+70) % MyConst::HTABSIZE;
        // startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
        // endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
        // tit = ref.strandTable.begin() + ref.tabIndex[hVal];
        // for (auto it = startIt; it != endIt; ++it, ++tit)
        // {
        //     KMER_S::kmer& k = *it;
        //     const uint32_t m = KMER_S::getMetaCpG(k);
        //     const bool isStart = KMER_S::isStartCpG(k);
        //     if (!isStart)
        //     {
        //         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
        //         if (*tit)
        //         {
        //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
        //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
        //             of << std::string(stIt, enIt) << "\n";
        //         }
        //     }
        // }
        //
        // of << "\n\nReverse seq matches\n";
        // hVal = ntHash::NTP64(revSeq1.data()) % MyConst::HTABSIZE;
        // startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
        // endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
        // tit = ref.strandTable.begin() + ref.tabIndex[hVal];
        // for (auto it = startIt; it != endIt; ++it, ++tit)
        // {
        //     KMER_S::kmer& k = *it;
        //     const uint64_t m = KMER_S::getMetaCpG(k);
        //     const bool isStart = KMER_S::isStartCpG(k);
        //     if (!isStart)
        //     {
        //         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
        //         if (*tit)
        //         {
        //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
        //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
        //             of << std::string(stIt, enIt) << "\n";
        //         }
        //     }
        // }
        // of << "Last Sequence part:\n";
        // hVal = ntHash::NTP64(revSeq1.data()+70) % MyConst::HTABSIZE;
        // startIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal];
        // endIt = ref.kmerTableSmall.begin() + ref.tabIndex[hVal + 1];
        // tit = ref.strandTable.begin() + ref.tabIndex[hVal];
        // for (auto it = startIt; it != endIt; ++it, ++tit)
        // {
        //     KMER_S::kmer& k = *it;
        //     const uint64_t m = KMER_S::getMetaCpG(k);
        //     const bool isStart = KMER_S::isStartCpG(k);
        //     if (!isStart)
        //     {
        //         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
        //         if (*tit)
        //         {
        //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
        //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + 2*MyConst::READLEN - 2 + startCpg.pos;
        //             of << std::string(stIt, enIt) << "\n";
        //         }
        //     }
        // }
        // of << "\n\n--------------------\n\n";


        // Check if no pairing possible
        if (bestErrNum == 2*MyConst::MISCOUNT + 1)
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

            if (extractSingleMatch(matches1Fwd, matches1Rev, r1, revSeq1))
            {
                ++succMatchT;
// #pragma omp critical
// {
//                 of << "\tSuccessfull r1\n";
// }

            } else {

                matches1Fwd.size() + matches1Rev.size() > 0 ? ++nonUniqueMatchT : ++unSuccMatchT;
// #pragma omp critical
// {
//                 of << "\tUnsuccessfull r1\n";
// }
            }
            if (extractSingleMatch(matches2Fwd, matches2Rev, r2, revSeq2))
            {

                ++succMatchT;
// #pragma omp critical
// {
//                 of << "\tSuccessfull r2\n";
// }
            } else {

                matches2Fwd.size() + matches2Rev.size() > 0 ? ++nonUniqueMatchT : ++unSuccMatchT;
// #pragma omp critical
// {
//                 of << "\tUnsuccessfull r1\n";
// }
            }

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
            nonUniqueMatchT += 2;
            r1.isInvalid = true;
            r2.isInvalid = true;

        } else {

            r1.mat = bestMatch1;
            r2.mat = bestMatch2;
            computeMethLvl(r1.mat, r1.seq);
            computeMethLvl(r2.mat, r2.seq);
            ++succPairedMatchT;
            succMatchT += 2;
        }
    }

    // sum up counts
    for (unsigned int i = 0; i < CORENUM; ++i)
    {
        succMatch += matchStats[i];
        nonUniqueMatch += nonUniqueStats[i];
        unSuccMatch += noMatchStats[i];
        succPairedMatch += matchPairedStats[i];
    }
    return true;
}



void ReadQueue::printMethylationLevels(std::string& filename)
{

    std::cout << "\nStart writing Methylation levels to \"" << filename << "_cpg_start.tsv\" and \"" << filename << "_cpg.tsv\"\n\n";
    std::ofstream cpgFile(filename + "_cpg_start.tsv");

    // go over CpGs close to start of a chromosome
    for (size_t cpgID = 0; cpgID < ref.cpgStartTable.size(); ++cpgID)
    {

        // print the position of the (C of the) CpG
        cpgFile << static_cast<uint16_t>(ref.cpgStartTable[cpgID].chrom) << "\t" << ref.cpgStartTable[cpgID].pos << "\t";

        // print the counts
        // fwd counts
        cpgFile << methLevelsStart[cpgID].methFwd << "\t" << methLevelsStart[cpgID].unmethFwd << "\t";
        // rev counts
        cpgFile << methLevelsStart[cpgID].methRev << "\t" << methLevelsStart[cpgID].unmethRev << "\n";
    }
    cpgFile.close();

    cpgFile.open(filename + "_cpg.tsv");

    // go over remaining CpGs
    for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
    {

        // print the position of the (C of the) CpG
        cpgFile << static_cast<uint16_t>(ref.cpgTable[cpgID].chrom) << "\t" << ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 << "\t";

        // print the counts
        // fwd counts
        cpgFile << methLevels[cpgID].methFwd << "\t" << methLevels[cpgID].unmethFwd << "\t";
        // rev counts
        cpgFile << methLevels[cpgID].methRev << "\t" << methLevels[cpgID].unmethRev << "\n";
    }
    cpgFile.close();
    std::cout << "Finished writing methylation levels to file\n\n";
}


// void ReadQueue::printStatistics(std::vector<std::vector<KMER::kmer> > seedsK)
// {
//
//     unsigned int count = 0;
//
//     // statistics on metaCpG repetitions
//     std::vector<unsigned int> metaRep (n, 0);
//
//     // count occurences of meta CpGs
//     for (unsigned int i = 0; i < seedsK.size(); ++i)
//     {
//
//         std::unordered_map<uint32_t, unsigned int> metaCounts;
//
//         // count up general seed counter
//         count += seedsK[i].size();
//         for (unsigned int j = 0; j < seedsK[i].size(); ++j)
//         {
//
//             uint64_t metaId = KMER::getMetaCpG(seedsK[i][j]);
//
//             auto insert = metaCounts.emplace(metaId, 1);
//
//             // if insertion fails because metaCpG is already inserted, count up everything
//             if (!insert.second)
//             {
//                 ++((insert.first)->second);
//             }
//         }
//         for (auto& c : metaCounts)
//             continue;
//         }
//         // i.e. check for all pairs the range criterion
//         // TODO:
//         // make timer for this
//         // for (MATCH::match mat1 : matches1)
//         // {
//         //     for (MATCH::match mat2 : matches2)
//         //     {
//         //         // TODO:
//         //         // Merge paired read matches
//         //     }
//         // }
//     }
//
//     // sum up counts
//     for (unsigned int i = 0; i < CORENUM; ++i)
//     {
//         succMatch += matchStats[i];
//         nonUniqueMatch += nonUniqueStats[i];
//         unSuccMatch += noMatchStats[i];
//     }
//     return true;
// }
//
// void ReadQueue::printMethylationLevels(std::string& filename)
// {
//
//     std::cout << "\nStart writing Methylation levels to \"" << filename << "_cpg_start.tsv\" and \"" << filename << "_cpg.tsv\"\n\n";
//     std::ofstream cpgFile(filename + "_cpg_start.tsv");
//
//     // go over CpGs close to start of a chromosome
//     for (size_t cpgID = 0; cpgID < ref.cpgStartTable.size(); ++cpgID)
//     {
//
//         // print the position of the (C of the) CpG
//         cpgFile << static_cast<uint16_t>(ref.cpgStartTable[cpgID].chrom) << "\t" << ref.cpgStartTable[cpgID].pos << "\t";
//
//         // print the counts
//         // fwd counts
//         cpgFile << methLevelsStart[cpgID].methFwd << "\t" << methLevelsStart[cpgID].unmethFwd << "\t";
//         // rev counts
//         cpgFile << methLevelsStart[cpgID].methRev << "\t" << methLevelsStart[cpgID].unmethRev << "\n";
//     }
//     cpgFile.close();
//
//     cpgFile.open(filename + "_cpg.tsv");
//
//     // go over remaining CpGs
//     for (size_t cpgID = 0; cpgID < ref.cpgTable.size(); ++cpgID)
//     {
//
//         // print the position of the (C of the) CpG
//         cpgFile << static_cast<uint16_t>(ref.cpgTable[cpgID].chrom) << "\t" << ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 << "\t";
//
//         // print the counts
//         // fwd counts
//         cpgFile << methLevels[cpgID].methFwd << "\t" << methLevels[cpgID].unmethFwd << "\t";
//         // rev counts
//         cpgFile << methLevels[cpgID].methRev << "\t" << methLevels[cpgID].unmethRev << "\n";
//     }
//     cpgFile.close();
//     std::cout << "Finished writing methylation levels to file\n\n";
// }
//
//
// void ReadQueue::printStatistics(std::vector<std::vector<KMER::kmer> > seedsK)
// {
//
//     unsigned int count = 0;
//
//     // statistics on metaCpG repetitions
//     std::vector<unsigned int> metaRep (n, 0);
//
//     // count occurences of meta CpGs
//     for (unsigned int i = 0; i < seedsK.size(); ++i)
//     {
//
//         std::unordered_map<uint32_t, unsigned int> metaCounts;
//
//         // count up general seed counter
//         count += seedsK[i].size();
//         for (unsigned int j = 0; j < seedsK[i].size(); ++j)
//         {
//
//             uint64_t metaId = KMER::getMetaCpG(seedsK[i][j]);
//
//             auto insert = metaCounts.emplace(metaId, 1);
//
//             // if insertion fails because metaCpG is already inserted, count up everything
//             if (!insert.second)
//             {
//                 ++((insert.first)->second);
//             }
//         }
//         for (auto& c : metaCounts)
//         {
//             if (c.second <= n)
//             {
//
//                 ++metaRep[c.second - 1];
//
//             }
//         }
//
//     }
//
//     // print everything
//     for (unsigned int i = 0; i < metaRep.size(); ++i)
//     {
//         statFile << metaRep[i] << "\t";
//     }
//     statFile << "\n";
//
//     countFile << count << "\t";
//
// }
