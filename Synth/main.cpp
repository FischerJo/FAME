
#include <fstream>
#include <iostream>

#include "SynthDS.h"



// --------------- MAIN -----------------
//
//
int main(int argc, char** argv)
{

    if (argc < 2)
    {
        std::cerr << "Output filename required! Terminating...\n";
        exit(1);
    }
    // SynthDS synthGen(refLen);
    SynthDS synthGen(argv[1], mthRate);

    std::vector<std::pair<size_t,size_t> > offsets;
    // error offsets
    std::vector<std::array<int, errNum> > errOffs;
    std::vector<std::string> fwdReads = synthGen.genReadsFwdRef(readLen, readNum, errNum, offsets, errOffs);



    std::ofstream ofsReads(std::string(argv[2]) + "_fwd.fastq");

    for (size_t i = 0; i < fwdReads.size(); ++i)
    {
        // generate fastq format of reads
        ofsReads << '@' << i;
        for (unsigned int e = 0; e < errNum; ++e)
            ofsReads << "_" << errOffs[i][e];
        ofsReads << "_CHR" << offsets[i].second << "_" << offsets[i].first << "\n";
        ofsReads << fwdReads[i] << "\n";
        ofsReads << '+' << i << "\n";
        // produce dummy quality scores
        for (size_t i = 0; i < readLen; ++i)
            ofsReads << "%";
        ofsReads << "\n";
    }
    ofsReads.close();
    offsets.clear();
    errOffs.clear();



    std::vector<std::string> revReads = synthGen.genReadsRevRef(readLen, readNum, errNum, offsets, errOffs);
    ofsReads.open(std::string(argv[2]) + "_rev.fastq");
    for (size_t i = 0; i < revReads.size(); ++i)
    {
        // generate fastq format of reads
        ofsReads << '@' << i;
        for (unsigned int e = 0; e < errNum; ++e)
            ofsReads << "_" << errOffs[i][e];
        ofsReads << "_CHR" << offsets[i].second << "_" << offsets[i].first << "\n";
        ofsReads << revReads[i] << "\n";
        ofsReads << '+' << i << "\n";
        // produce dummy quality scores
        for (size_t i = 0; i < readLen; ++i)
            ofsReads << "%";
        ofsReads << "\n";
    }
    ofsReads.close();
    ofsReads.open(std::string(argv[2]) + "_cpginfo_fwd.tsv");
    ofsReads << "Chromosome\tPosition\tUnmethylated\tMethylated\n";
    for (auto& cpg : synthGen.cpgMethRateFwd)
    {
        ofsReads << (cpg.first >> 32) << "\t" << (cpg.first & 0x00000000ffffffffULL) << "\t" << cpg.second.first << "\t" << cpg.second.second << "\n";
    }
    ofsReads.close();
    ofsReads.open(std::string(argv[2]) + "_cpginfo_rev.tsv");
    ofsReads << "Chromosome\tPosition\tUnmethylated\tMethylated\n";
    for (auto& cpg : synthGen.cpgMethRateRev)
    {
        ofsReads << (cpg.first >> 32) << "\t" << (cpg.first & 0x00000000ffffffffULL) << "\t" << cpg.second.first << "\t" << cpg.second.second << "\n";
    }
    return 0;

}
