#include <iostream>
#include <chrono>

#include "ReadQueue.h"

ReadQueue::ReadQueue(const char* filePath, RefGenome& reference, bool isGZ) :
        statFile("seedStats.tsv")
    // ,   countFile("seedCounts.tsv")
    ,   ref(reference)
    ,   readBuffer(MyConst::CHUNKSIZE)
        // TODO
    // ,   of  ("meta_info_human_k30_m1.txt")
{
    if (isGZ)
    {

        igz.open(filePath);

    } else {

        file.open(filePath);
    }

    // fill counting structure for parallelization
    for (unsigned int i = 0; i < CORENUM; ++i)
    {

        fwdMetaIDs[i] = std::unordered_map<uint32_t, uint16_t, MetaHash>();
        revMetaIDs[i] = std::unordered_map<uint32_t, uint16_t, MetaHash>();
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
    // statFile.close();
}

bool ReadQueue::parseChunk(unsigned int& procReads)
{

    std::cout << "Start reading chunk of reads\n";

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
            return true;

        }
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
            return true;

        }
    }

    procReads = readCounter;
    return false;
}

bool ReadQueue::matchReads(const unsigned int& procReads)
{

    // TODO
    uint64_t succMatchFwd = 0;
    uint64_t succMatchRev = 0;
    uint64_t unSuccMatch = 0;
    uint64_t nonUniqueMatch = 0;
    uint64_t readCount = 0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(CORENUM) schedule(dynamic,1000)
#endif
    for (unsigned int i = 0; i < procReads; ++i)
    {

        Read& r = readBuffer[i];

        const size_t readSize = r.seq.size();

        if (readSize < MyConst::KMERLEN)
        {

            r.isInvalid = true;
            continue;
        }
        // RETRIEVE SEEDS
        //
        // std::vector<char> redSeq (readSize);
        // std::vector<char> redRevSeq (readSize);

        // flag stating if read contains N
        // reads with N are ignored
        bool nflag = false;

        // get correct offset for reverse strand (strand orientation must be correct)
        size_t revPos = readSize - 1;

        std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
        // string containing reverse complement (under FULL alphabet)
        std::string revSeq;
        revSeq.resize(readSize);

        bool isCpG = false;
        bool hadCorT = false;
        // construct reduced alphabet sequence for forward and reverse strand
        for (size_t pos = 0; pos < readSize; ++pos, --revPos)
        {

            switch (r.seq[pos])
            {
                case 'A':

                    // redSeq[pos] = 'A';
                    // redRevSeq[revPos] = 'T';
                    revSeq[revPos] = 'T';
                    hadCorT = false;
                    break;

                case 'C':

                    // redSeq[pos] = 'T';
                    // redRevSeq[revPos] = 'G';
                    revSeq[revPos] = 'G';
                    hadCorT = true;
                    break;

                case 'G':

                    // redSeq[pos] = 'G';
                    // redRevSeq[revPos] = 'T';
                    revSeq[revPos] = 'C';
                    hadCorT ? isCpG = true : false;
                    hadCorT = false;
                    break;

                case 'T':

                    // redSeq[pos] = 'T';
                    // redRevSeq[revPos] = 'A';
                    revSeq[revPos] = 'A';
                    hadCorT = true;
                    break;

                case 'N':

                    nflag = true;
                    break;

                default:

                    std::cerr << "Unknown character '" << r.seq[pos] << "' in read with sequence id " << r.id << std::endl;
            }
        }

        // if N is present in read or read has no CpG, skip
        if (nflag || !isCpG)
        {
            r.isInvalid = true;
            continue;
        }

        std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        // of << runtime << "\t";
        // startTime = std::chrono::high_resolution_clock::now();

        // std::vector<std::vector<KMER::kmer> > fwdSeedsK;
        // std::vector<std::vector<bool> > fwdSeedsS;
        // ref.getSeeds(redSeq, fwdSeedsK, fwdSeedsS);
        // std::vector<std::vector<KMER::kmer> > revSeedsK;
        // std::vector<std::vector<bool> > revSeedsS;
        // ref.getSeeds(redRevSeq, revSeedsK, revSeedsS);

        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // of << runtime << "\t";

        // printStatistics(fwdSeedsK);

        // startTime = std::chrono::high_resolution_clock::now();
        // // FILTER SEEDS BY COUNTING LEMMA
        // filterHeuSeeds(fwdSeedsK, fwdSeedsS, readSize);
        // filterHeuSeeds(revSeedsK, revSeedsS, readSize);
        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // of << runtime << "\t";

        // startTime = std::chrono::high_resolution_clock::now();
        // produce shift-and automaton for forward and reverse sequence of this read
        // ShiftAnd<MyConst::MISCOUNT> saFwd(r.seq, lmap);
        // ShiftAnd<MyConst::MISCOUNT> saRev(revSeq, lmap);

        // match
        // note that this is only used when initialized by saQuerySeedSet, based on boolean return of this function
        // MATCH::match matchFwd = 0;
        // MATCH::match matchRev = 0;
        // query seeds to shift-and automaton

        // int succQueryFwd = saQuerySeedSet(saFwd, fwdSeedsK, fwdSeedsS, matchFwd);


        // TODO
        startTime = std::chrono::high_resolution_clock::now();
        MATCH::match matchFwd = 0;
        bool succFlag = getSeedRefs(r.seq, readSize);
        if (!succFlag)
        {
            ++readCount;
            // of << "\n";
            r.isInvalid = true;
            continue;
        }
        endTime = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // if (runtime > 100)
        // {
        //     of << r.seq << "\n";
        //     of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
        // }
        // of << runtime << "\t";
        startTime = std::chrono::high_resolution_clock::now();
        ShiftAnd<MyConst::MISCOUNT> saFwd(r.seq, lmap);
        endTime = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // of << runtime << "\t";
        startTime = std::chrono::high_resolution_clock::now();
        int succQueryFwd = saQuerySeedSetRef(saFwd, matchFwd, computeQgramThresh(r.seq));
        endTime = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // if (runtime > 100)
        // {
        //     of << r.seq << "\n";
        //     of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
        //     uint64_t qcount = 0;
        //     for (const auto& m : fwdMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= 40)
        //             ++qcount;
        //     }
        //     for (const auto& m : revMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= 40)
        //             ++qcount;
        //     }
        //     of << "Meta CpGs passing q-gram filter: " << qcount << "\n";
        // }
        // of << runtime << "\t";

        startTime = std::chrono::high_resolution_clock::now();
        MATCH::match matchRev = 0;
        succFlag = getSeedRefs(revSeq, readSize);
        if (!succFlag)
        {
            ++readCount;
            // of << "\n";
            r.isInvalid = true;
            continue;
        }
        endTime = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // if (runtime > 100)
        // {
        //     of << r.seq << "\n";
        //     of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
        // }
        // of << runtime << "\t";
        startTime = std::chrono::high_resolution_clock::now();
        ShiftAnd<MyConst::MISCOUNT> saRev(revSeq, lmap);
        endTime = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // of << runtime << "\t";
        startTime = std::chrono::high_resolution_clock::now();
        int succQueryRev = saQuerySeedSetRef(saRev, matchRev, computeQgramThresh(revSeq));
        endTime = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // if (runtime > 100)
        // {
        //     of << r.seq << "\n";
        //     of << "Overall meta CpGs: " << fwdMetaIDs[omp_get_thread_num()].size() + revMetaIDs[omp_get_thread_num()].size() << "\n";
        //     uint64_t qcount = 0;
        //     for (const auto& m : fwdMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= 40)
        //             ++qcount;
        //     }
        //     for (const auto& m : revMetaIDs[omp_get_thread_num()])
        //     {
        //         if (m.second >= 40)
        //             ++qcount;
        //     }
        //     of << "Meta CpGs passing q-gram filter: " << qcount << "\n";
        // }
        // of << runtime << "\n";


        // of << "\n\n--------------------\n\n";






        // startTime = std::chrono::high_resolution_clock::now();
        // int succQueryRev = saQuerySeedSet(saRev, revSeedsK, revSeedsS, matchRev);
        // endTime = std::chrono::high_resolution_clock::now();
        // runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        // of <<  runtime << "\n";
        // if (runtime > 100)
        // {
        //     of << "REV---:" << r.seq << std::endl;
        //     for (auto s : revSeedsK)
        //     {
        //
        //         of << s.size() << " ; ";
        //     }
        //
        //     of << "\n\n";
        // }

        // found match for fwd and rev strand
        if (succQueryFwd == 1 && succQueryRev == 1)
        {

            uint8_t fwdErr = MATCH::getErrNum(matchFwd);
            uint8_t revErr = MATCH::getErrNum(matchRev);

            // check which one has fewer errors
            if (fwdErr < revErr)
            {

                ++succMatchFwd;
                r.mat = matchFwd;

            } else {

                if (fwdErr > revErr)
                {

                    ++succMatchRev;
                    r.mat = matchRev;

                // if same number of errors, then not unique
                } else {

                    // of << MATCH::getOffset(matchFwd) << " match fwdAutomaton " << MATCH::getErrNum(matchFwd) << "errors --- strand " << MATCH::isFwd << "\n" << MATCH::getOffset(matchRev) << " match revAutomaton " << MATCH::getErrNum(matchRev) << "errors --- strand: " << MATCH::isFwd << "\n";

                    const uint64_t offFwd = MATCH::getOffset(matchFwd);
                    const uint64_t offRev = MATCH::getOffset(matchRev);
                    // est if same match in same region
                    if (offFwd == offRev)
                    {
                        ++succMatchFwd;
                        r.mat = matchFwd;

                    } else {

                        ++nonUniqueMatch;

                        // std::cout << "Nonunique match: fwd and rev with same errorcount in match\n";
                        r.isInvalid = true;
                    }
                }
            }
        // unique match on forward strand
        } else if (succQueryFwd == 1) {

            ++succMatchFwd;
            r.mat = matchFwd;

        // unique match on backward strand
        } else if (succQueryRev == 1) {

            ++succMatchRev;
            r.mat = matchRev;

        // no match found at all
        } else {

            r.isInvalid = true;
            // TODO
            if (succQueryFwd == -1 || succQueryRev == -1)
            {

                ++nonUniqueMatch;

            } else {

                ++unSuccMatch;

                // construct hash and look up the hash table entries
                // of << "\nreal seq: " << r.seq << "\n" << r.id << "\n";
                // uint64_t hVal = ntHash::NTP64(r.seq.data()) % MyConst::HTABSIZE;
                // auto startIt = ref.kmerTable.begin() + ref.tabIndex[hVal];
                // auto endIt = ref.kmerTable.begin() + ref.tabIndex[hVal + 1];
                // auto tit = ref.strandTable.begin() + ref.tabIndex[hVal];
                // for (auto it = startIt; it != endIt; ++it, ++tit)
                // {
                //     KMER::kmer& k = *it;
                //     const uint64_t m = KMER::getMetaCpG(k);
                //     const bool isStart = KMER::isStartCpG(k);
                //     if (!isStart)
                //     {
                //         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
                //         uint64_t offset = KMER::getOffset(k);
                //         if (*tit)
                //         {
                //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + offset + startCpg.pos;
                //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + offset + MyConst::KMERLEN + startCpg.pos;
                //             // auto enIt = ref.fullSeq[startCpg.chrom].begin() + offset + startCpg.pos + 100;
                //             of << std::string(stIt, enIt) << " offset in meta " << offset << " ||| offset in seq " << offset + startCpg.pos << " ||| metaCpG " << m << "\n";
                //         } else {
                //
                //             of << "rev off in meta " << offset << " ||| rev off in seq " << offset + startCpg.pos << "\n";
                //         }
                //     }
                // }
                // of << "Last Sequence part:\n";
                // hVal = ntHash::NTP64(r.seq.data()+69) % MyConst::HTABSIZE;
                // startIt = ref.kmerTable.begin() + ref.tabIndex[hVal];
                // endIt = ref.kmerTable.begin() + ref.tabIndex[hVal + 1];
                // tit = ref.strandTable.begin() + ref.tabIndex[hVal];
                // for (auto it = startIt; it != endIt; ++it, ++tit)
                // {
                //     KMER::kmer& k = *it;
                //     const uint64_t m = KMER::getMetaCpG(k);
                //     const bool isStart = KMER::isStartCpG(k);
                //     if (!isStart)
                //     {
                //         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
                //         uint64_t offset = KMER::getOffset(k);
                //         if (*tit)
                //         {
                //             auto stIt = ref.fullSeq[startCpg.chrom].begin() + offset + startCpg.pos;
                //             auto enIt = ref.fullSeq[startCpg.chrom].begin() + offset + MyConst::KMERLEN + startCpg.pos;
                //             of << std::string(stIt, enIt) << " offset in meta " << offset << " ||| offset in seq " << offset + startCpg.pos << " ||| metaCpG " << m << "\n";
                //         } else {
                //
                //             of << "rev off in meta " << offset << " ||| rev off in seq " << offset + startCpg.pos << "\n";
                //         }
                //     }
            //     }
            //
            //     of << "\n\nreal reverse seq: " << revSeq << "\n";
            //     hVal = ntHash::NTP64(revSeq.data()) % MyConst::HTABSIZE;
            //     startIt = ref.kmerTable.begin() + ref.tabIndex[hVal];
            //     endIt = ref.kmerTable.begin() + ref.tabIndex[hVal + 1];
            //     tit = ref.strandTable.begin() + ref.tabIndex[hVal];
            //     for (auto it = startIt; it != endIt; ++it, ++tit)
            //     {
            //         KMER::kmer& k = *it;
            //         const uint64_t m = KMER::getMetaCpG(k);
            //         const bool isStart = KMER::isStartCpG(k);
            //         if (!isStart)
            //         {
            //             const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
            //             uint64_t offset = KMER::getOffset(k);
            //             if (*tit)
            //             {
            //                 auto stIt = ref.fullSeq[startCpg.chrom].begin() + offset + startCpg.pos;
            //                 auto enIt = ref.fullSeq[startCpg.chrom].begin() + offset + MyConst::KMERLEN + startCpg.pos;
            //                 of << std::string(stIt, enIt) << " offset in meta " << offset << " ||| offset in seq " << offset + startCpg.pos << " ||| metaCpG " << m << "\n";
            //             } else {
            //
            //                 of << "rev off in meta " << offset << " ||| rev off in seq " << offset + startCpg.pos << "\n";
            //             }
            //         }
            //     }
            //     of << "Last Sequence part:\n";
            //     hVal = ntHash::NTP64(revSeq.data()+69) % MyConst::HTABSIZE;
            //     startIt = ref.kmerTable.begin() + ref.tabIndex[hVal];
            //     endIt = ref.kmerTable.begin() + ref.tabIndex[hVal + 1];
            //     tit = ref.strandTable.begin() + ref.tabIndex[hVal];
            //     for (auto it = startIt; it != endIt; ++it, ++tit)
            //     {
            //         KMER::kmer& k = *it;
            //         const uint64_t m = KMER::getMetaCpG(k);
            //         const bool isStart = KMER::isStartCpG(k);
            //         if (!isStart)
            //         {
            //             const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
            //             uint64_t offset = KMER::getOffset(k);
            //             if (*tit)
            //             {
            //                 auto stIt = ref.fullSeq[startCpg.chrom].begin() + offset + startCpg.pos;
            //                 auto enIt = ref.fullSeq[startCpg.chrom].begin() + offset + MyConst::KMERLEN + startCpg.pos;
            //                 of << std::string(stIt, enIt) << " offset in meta " << offset << " ||| offset in seq " << offset + startCpg.pos << " ||| metaCpG " << m << "\n";
            //             } else {
            //
            //                 of << "rev off in meta " << offset << " ||| rev off in seq " << offset + startCpg.pos << "\n";
            //             }
            //         }
            //     }
            //     of << "\n\n---------------------------------\n\n";
            }

            // if (unSuccMatch > 10)
            // {
            //     of.close();
            //     exit(1);
            // }
        }
        // validate
        // if (!r.isInvalid)
        // {
            // std::string stringOffset (r.id.begin() + r.id.find_first_of('_', 1) + 1, r.id.end());
            // unsigned long offset = std::stoul(stringOffset);
            // uint64_t matchedOff = MATCH::getOffset(r.mat);
            // // if (!MATCH::isFwd(r.mat))
            // // {
            // //     if (MATCH::getErrNum(r.mat) > 0)
            // //     {
            // //         of << "This should be a non-unique match. Bad.\n";
            // //         of << "Number of errors: " << MATCH::getErrNum(r.mat) << "\n";
            // //         of << "Offset should be " << offset << " but is " << matchedOff << "\n";
            // //         continue;
            // //     }
            // // }
            // if (matchedOff > 100 && (matchedOff > offset + 110 || matchedOff < offset))
            // {
            //     ++pCount;
            //     std::cout << "Wrong match. Offset should be " << offset << " but is " << matchedOff << "\n";
            //     std::cout << "Number of errors: " << MATCH::getErrNum(r.mat) << "\n";
            //     std::cout << "Read sequence:\n" << r.seq << "\n\n";
            //     std::cout << "Reference sequence (match):\n" << std::string(ref.fullSeq[0].begin() + matchedOff - 99, ref.fullSeq[0].begin() + matchedOff) << "\n\n";
            //     std::cout << "Reference sequence (source):\n: " << std::string(ref.fullSeq[0].begin() + offset, ref.fullSeq[0].begin() + offset + 100) << "\n\n";
            // }
            // if (pCount > 5)
            // {
            //
            //     of.close();
            //     exit(1);
            // }
        // }

        // BITMASK COMPARISON ON FULL ALPHABET FOR REMAINING SEEDS
        // bitMatching(r, fwdSeedsK, fwdSeedsS);
        // bitMatchingRev(r, revSeedsK, revSeedsS);

    }

    of.close();
    std::cout << "Successfully matched (Fwd|Rev): " << succMatchFwd << "|" << succMatchRev << " / Unsuccessfully matched: " << unSuccMatch << " / Nonunique matches: " << nonUniqueMatch << " / Discarded because of too many meta CpGs: " << readCount << "\n\n";
    // TODO: Go over read set once and register the CpG matchings
    return true;
}


void ReadQueue::printStatistics(std::vector<std::vector<KMER::kmer> > seedsK)
{

    unsigned int count = 0;

    // statistics on metaCpG repetitions
    std::vector<unsigned int> metaRep (n, 0);

    // count occurences of meta CpGs
    for (unsigned int i = 0; i < seedsK.size(); ++i)
    {

        std::unordered_map<uint32_t, unsigned int> metaCounts;

        // count up general seed counter
        count += seedsK[i].size();
        for (unsigned int j = 0; j < seedsK[i].size(); ++j)
        {

            uint64_t metaId = KMER::getMetaCpG(seedsK[i][j]);

            auto insert = metaCounts.emplace(metaId, 1);

            // if insertion fails because metaCpG is already inserted, count up everything
            if (!insert.second)
            {
                ++((insert.first)->second);
            }
        }
        for (auto& c : metaCounts)
        {
            if (c.second <= n)
            {

                ++metaRep[c.second - 1];

            }
        }

    }

    // print everything
    for (unsigned int i = 0; i < metaRep.size(); ++i)
    {
        statFile << metaRep[i] << "\t";
    }
    statFile << "\n";

    countFile << count << "\t";

}
