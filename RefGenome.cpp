
#include <chrono>

#include "RefGenome.h"


RefGenome::RefGenome(std::vector<struct CpG>&& cpgTab, std::vector<struct CpG>&& cpgStartTab, std::vector<std::vector<char> >& genomeSeq) :
        cpgTable(cpgTab)
    ,   cpgStartTable(cpgStartTab)
    ,   genomeBit()
    ,   fullSeq(std::move(genomeSeq))
    ,   tabIndex(MyConst::HTABSIZE + 1, 0)
    ,   kmerTable()
    ,   strandTable()
    ,   metaCpGs()
    ,   metaStartCpGs()
{
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
    // filter out highly repetitive sequences
    // filterHashTable();
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
    std::ofstream of(filepath, std::ofstream::binary);

    // save CONSTANTS
    auto htabs = MyConst::HTABSIZE;
    of.write(reinterpret_cast<char*>(&htabs), sizeof(htabs));
    auto readl = MyConst::READLEN;
    of.write(reinterpret_cast<char*>(&readl), sizeof(readl));
    auto winl = MyConst::WINLEN;
    of.write(reinterpret_cast<char*>(&winl), sizeof(winl));
    auto kmerl = MyConst::KMERLEN;
    of.write(reinterpret_cast<char*>(&kmerl), sizeof(kmerl));

    // store NON-start CpGs
    size_t cpgNum = cpgTable.size();
    of.write(reinterpret_cast<char*>(&cpgNum), sizeof(cpgNum));
    for (auto cpg : cpgTable)
    {
        uint8_t chr = cpg.chrom;
        of.write(reinterpret_cast<char*>(&chr), sizeof(chr));
        uint32_t pos = cpg.pos;
        of.write(reinterpret_cast<char*>(&pos), sizeof(pos));
    }
    // store start CpGs
    cpgNum = cpgStartTable.size();
    of.write(reinterpret_cast<char*>(&cpgNum), sizeof(cpgNum));
    for (auto cpg : cpgStartTable)
    {
        uint8_t chr = cpg.chrom;
        of.write(reinterpret_cast<char*>(&chr), sizeof(chr));
        uint32_t pos = cpg.pos;
        of.write(reinterpret_cast<char*>(&pos), sizeof(pos));
    }
    // write reference sequence
    size_t chromNum = fullSeq.size();
    of.write(reinterpret_cast<char*>(&chromNum), sizeof(chromNum));
    for (std::vector<char> chromSeq : fullSeq)
    {
        size_t chromLen = chromSeq.size();
        of.write(reinterpret_cast<char*>(&chromLen), sizeof(chromLen));
        of.write(chromSeq.data(), chromLen);
    }
    // tabIndex
    size_t tabLen = tabIndex.size();
    of.write(reinterpret_cast<char*>(&tabLen), sizeof(tabLen));
    of.write(reinterpret_cast<char*>(tabIndex.data()), sizeof(tabIndex[0])*tabLen);

    // store kmers
    size_t kmerNum = kmerTable.size();
    of.write(reinterpret_cast<char*>(&kmerNum), sizeof(kmerNum));
    of.write(reinterpret_cast<char*>(kmerTable.data()), sizeof(kmerTable[0])*kmerNum);

    // store strands
    std::copy(strandTable.begin(), strandTable.end(), std::ostreambuf_iterator<char>(of));

    // store meta CpGs
    size_t metaCpGNum = metaCpGs.size();
    of.write(reinterpret_cast<char*>(&metaCpGNum), sizeof(metaCpGNum));
    for (auto m : metaCpGs)
    {
        uint32_t start = m.start;
        of.write(reinterpret_cast<char*>(&start), sizeof(start));
        uint32_t end = m.end;
        of.write(reinterpret_cast<char*>(&end), sizeof(end));
    }

    // store start meta CpGs
    metaCpGNum = metaStartCpGs.size();
    of.write(reinterpret_cast<char*>(&metaCpGNum), sizeof(metaCpGNum));
    for (auto m : metaStartCpGs)
    {
        uint32_t start = m.start;
        of.write(reinterpret_cast<char*>(&start), sizeof(start));
        uint32_t end = m.end;
        of.write(reinterpret_cast<char*>(&end), sizeof(end));
    }
    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    std::cout << "Finished writing index to file in " << runtime << "s\n\n";
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

    // load NON-start CpGs
    size_t cpgNum;
    ifs.read(reinterpret_cast<char*>(&cpgNum), sizeof(cpgNum));
    cpgTable.reserve(cpgNum);
    for (size_t i = 0; i < cpgNum; ++i)
    {
        uint8_t chr;
        ifs.read(reinterpret_cast<char*>(&chr), sizeof(chr));
        uint32_t pos;
        ifs.read(reinterpret_cast<char*>(&pos), sizeof(pos));
        cpgTable.push_back({chr, pos});
    }
    // load start CpGs
    ifs.read(reinterpret_cast<char*>(&cpgNum), sizeof(cpgNum));
    cpgStartTable.reserve(cpgNum);
    for (size_t i = 0; i < cpgNum; ++i)
    {
        uint8_t chr;
        ifs.read(reinterpret_cast<char*>(&chr), sizeof(chr));
        uint32_t pos;
        ifs.read(reinterpret_cast<char*>(&pos), sizeof(pos));
        cpgTable.push_back({chr, pos});
    }
    // write reference sequence
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
    kmerTable.resize(kmerNum);
    ifs.read(reinterpret_cast<char*>(kmerTable.data()), sizeof(kmerTable[0])*kmerNum);

    // load strands
    strandTable.resize(kmerNum);
    strandTable.insert(strandTable.begin(), std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());

    // load meta CpGs
    size_t metaCpGNum;
    ifs.read(reinterpret_cast<char*>(&metaCpGNum), sizeof(metaCpGNum));
    metaCpGs.reserve(metaCpGNum);
    for (size_t i = 0; i < metaCpGNum; ++i)
    {
        uint32_t start;
        ifs.read(reinterpret_cast<char*>(&start), sizeof(start));
        uint32_t end;
        ifs.read(reinterpret_cast<char*>(&end), sizeof(end));
        metaCpGs.push_back({start, end});
    }

    // load start meta CpGs
    ifs.read(reinterpret_cast<char*>(&metaCpGNum), sizeof(metaCpGNum));
    metaStartCpGs.reserve(metaCpGNum);
    for (size_t i = 0; i < metaCpGNum; ++i)
    {
        uint32_t start;
        ifs.read(reinterpret_cast<char*>(&start), sizeof(start));
        uint32_t end;
        ifs.read(reinterpret_cast<char*>(&end), sizeof(end));
        metaStartCpGs.push_back({start, end});
    }
    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    std::cout << "Finished reading index to file in " << runtime << "s\n\n";
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

        while (cpgStartInd < cpgStartTable.size() && cpgStartTable[cpgStartInd].chrom == currChr)
        {

            ++cpgStartInd;
        }

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


void RefGenome::generateBitStrings(std::vector<std::vector<char> >& genomeSeq)
{

    unsigned int genSeqIndex = 0;
    // generate genome encodings with bitmasks
    for (std::vector<char>& seq : genomeSeq)
    {
        const unsigned int segLen = seq.size() / 32;
        genomeBit.emplace_back(seq.size());

        for (unsigned int i = 0; i < segLen; ++i)
        {

            genomeBit[genSeqIndex].setBitStrN( std::string(seq.data() + (i*32), 32), i);
        }
        if ((seq.size() % 32) != 0)
        {

            genomeBit[genSeqIndex].setBitStrLast( std::string(seq.data() + (segLen*32), seq.size() % 32));
        }
        ++genSeqIndex;
    }

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
            // if we can read the full sequence breadth after the CpG
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

void RefGenome::ntHashLast(const std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG, const uint32_t& metacpg, uint32_t&& metaOff)
{
    // kmers to be skipped
    const int skipKmer = max(lastPos - pos, 0);
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
    uint64_t rhVal = ntHash::NTP64(seqStartRev);
    // first kmer on reverse complement corresponds to last kmer in forward sequence
    uint64_t kPosRev = off - MyConst::KMERLEN;

    // update kmer table
    kmerTable[--tabIndex[rhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
    strandTable[tabIndex[rhVal % MyConst::HTABSIZE]] = false;


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        --kPosRev;
        ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
        // update kmer table
        kmerTable[--tabIndex[rhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPosRev + metaOff));
        strandTable[tabIndex[rhVal % MyConst::HTABSIZE]] = false;
    }

    // initial hash forward
    uint64_t fhVal = ntHash::NTP64(seqStart);
    uint64_t kPos = lasN;

    // update kmer table
    kmerTable[--tabIndex[fhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
    strandTable[tabIndex[fhVal % MyConst::HTABSIZE]] = true;

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ++kPos;
        ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        // update kmer table
        kmerTable[--tabIndex[fhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(0, metacpg, kPos + metaOff));
        strandTable[tabIndex[fhVal % MyConst::HTABSIZE]] = true;
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
    uint64_t rhVal = ntHash::NTP64(seqStartRev);
    // first kmer on reverse complement corresponds to last kmer in forward sequence
    uint64_t kPosRev = off - MyConst::KMERLEN;

    // update kmer table
    kmerTable[--tabIndex[rhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPosRev));
    strandTable[tabIndex[rhVal % MyConst::HTABSIZE]] = false;


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        --kPosRev;
        ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
        // update kmer table
        kmerTable[--tabIndex[rhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPosRev));
        strandTable[tabIndex[rhVal % MyConst::HTABSIZE]] = false;
    }

    // initial hash forward
    uint64_t fhVal = ntHash::NTP64(seqStart);
    uint64_t kPos = lasN;

    // update kmer table
    kmerTable[--tabIndex[fhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPos));
    strandTable[tabIndex[fhVal % MyConst::HTABSIZE]] = true;

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ++kPos;
        ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        // update kmer table
        kmerTable[--tabIndex[fhVal % MyConst::HTABSIZE]] = std::move(KMER::constructKmer(1, metacpg, kPos));
        strandTable[tabIndex[fhVal % MyConst::HTABSIZE]] = true;
    }
    lastPos = off - MyConst::KMERLEN + 1;

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
    uint64_t rhVal = ntHash::NTP64(seqStartRev);

    // update kmer table
    ++tabIndex[rhVal % MyConst::HTABSIZE];


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
        ++tabIndex[rhVal % MyConst::HTABSIZE];
    }

    // initial hash forward
    uint64_t fhVal = ntHash::NTP64(seqStart);

    ++tabIndex[fhVal % MyConst::HTABSIZE];

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        ++tabIndex[fhVal % MyConst::HTABSIZE];
    }
    lastPos = off - MyConst::KMERLEN + 1;

}
void RefGenome::ntCountLast(std::vector<char>& seq, uint32_t& lastPos, const unsigned int& pos, const unsigned int& bpsAfterCpG)
{
    // kmers to be skipped
    const int skipKmer = max(lastPos - pos, 0);
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
    uint64_t rhVal = ntHash::NTP64(seqStartRev);

    // update indices
    ++tabIndex[rhVal % MyConst::HTABSIZE];


    // hash kmers of backward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ntHash::NTP64(rhVal, seqStartRev[i], seqStartRev[MyConst::KMERLEN + i]);
        // update indices
        ++tabIndex[rhVal % MyConst::HTABSIZE];
    }

    // initial hash forward
    uint64_t fhVal = ntHash::NTP64(seqStart);

    // update indices
    ++tabIndex[fhVal % MyConst::HTABSIZE];

    // hash kmers of forward strand
    for (unsigned int i = 0; i < (contextLen - MyConst::KMERLEN); ++i)
    {
        ntHash::NTP64(fhVal, seqStart[i], seqStart[MyConst::KMERLEN + i]);
        // update indices
        ++tabIndex[fhVal % MyConst::HTABSIZE];
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
    strandTable.resize(filterItK - kmerTable.begin());

    std::chrono::high_resolution_clock::time_point filterEndTime = std::chrono::high_resolution_clock::now();

    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(filterEndTime - filterStartTime).count();


    std::cout << "Hash table size after filter (running " << runtime << "s): " << tabIndex[MyConst::HTABSIZE] << "\n\n";
}
