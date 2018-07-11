
#include <fstream>
#include <iostream>

#include "SynthDS.h"



// --------------- MAIN -----------------
//
//
int main(int argc, char** argv)
{

    // std::string genomeFile = "";
    // std::string outputFile = "out";
    //
    // unsigned int errNum = 0;
    // double methRate = 0.5;
    // double convRate = 1.0;
    // unsigned long readLen = 100;
    // unsigned long refLen = 1000000000;
    // unsigned long readNum = 10000000;
    //
    // // flags for parameters
    // bool hasErrNum = false;
    // bool hasMethRate = false;
    // bool hasConvRate = false;
    // bool hasReadLen = false;
    // bool hasRefLen = false;
    // bool hasReadNum = false;
    //
    // if (argc == 1)
    // {
    //     printHelp();
    //     return 0;
    // }
    //
    // for (int i = 1; i < argc; ++i)
    // {
    //
    //     if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help")
    //     {
    //         printHelp();
    //         return 0;
    //     }
    //
    //     if (std::string(argv[i]) == "--genome" || std::string(argv[i]) == "-g")
    //     {
    //         if (i + 1 < argc)
    //         {
    //             genomeFile = std::string(argv[++i]);
    //
    //         } else {
    //
    //             std::cerr << "No filepath for option \"" << argv[i] << "\" provided! Terminating...\n\n";
    //             exit(1);
    //         }
    //         continue;
    //     }
    //
    //     if (std::string(argv[i]) == "-o" || std::string(argv[i]) == "--out_basename")
    //     {
    //         if (i + 1 < argc)
    //         {
    //
    //             outputFile = argv[++i];
    //
    //         } else {
    //
    //             std::cerr << "No filepath for option \"" << argv[i] << "\" provided! Terminating...\n\n";
    //             exit(1);
    //         }
    //         continue;
    //     }
    //
    //     if (std::string(argv[i]) == "-e" || std::string(argv[i]) == "--max_err")
    //     {
    //         if (i + 1 < argc)
    //         {
    //             errNum = std::stoui(argv[i+1]);
    //
    //         } else {
    //
    //             std::cerr << "No integer following option \"" << argv[i] << "\"! Terminating...\n\n";
    //             exit(1);
    //         }
    //     }
    //
    //     if (std::string(argv[i]) == "--meth_rate")
    //     {
    //         if (i + 1 < argc)
    //         {
    //             methRate = std::stod(argv[i+1]);
    //
    //         } else {
    //
    //             std::cerr << "No probability following option \"" << argv[i] << "\"! Terminating...\n\n";
    //             exit(1);
    //         }
    //     }
    //
    //     if (std::string(argv[i]) == "--conv_rate")
    //     {
    //         if (i + 1 < argc)
    //         {
    //             convRate = std::stod(argv[i+1]);
    //
    //         } else {
    //
    //             std::cerr << "No probability following option \"" << argv[i] << "\"! Terminating...\n\n";
    //             exit(1);
    //         }
    //     }
    //
    //     // no such option
    //     std::cerr << "Don't know the option \"" << std::string(argv[i]) << "\", maybe you forgot a flag?\n\n";
    //     exit(1);
    //
    // }
    //
    //
    // Start processing
    SynthDS synthGen(argv[1], mthRate);



    if (argv[3] == std::string("P"))
    {
        std::cout << "Generating paired end reads.\n\n";
        std::pair<std::vector<std::pair<size_t,size_t> >, std::vector<std::pair<size_t,size_t> > >  offsets;
        std::pair<std::vector<std::list<int> >, std::vector<std::list<int> > > errOffs;
        std::pair<std::vector<std::string>, std::vector<std::string> > pairedReads = synthGen.genReadsPairedRef(readLen, readNum, offsets, errOffs);
        std::ofstream ofsReads1(std::string(argv[2]) + "_p1.fastq");
        std::ofstream ofsReads2(std::string(argv[2]) + "_p2.fastq");
        for (size_t i = 0; i < pairedReads.first.size(); ++i)
        {
            // generate fastq format of reads
            ofsReads1 << '@' << i;
            for (int errPos : errOffs.first[i])
                ofsReads1 << "_" << errPos;
            ofsReads1 << "_CHR" << offsets.first[i].second << "_" << offsets.first[i].first << "\n";
            ofsReads1 << pairedReads.first[i] << "\n";
            ofsReads1 << '+' << i << "\n";
            // produce dummy quality scores
            for (size_t j = 0; j < readLen; ++j)
                ofsReads1 << "%";
            ofsReads1 << "\n";
            // generate fastq format of reads
            ofsReads2 << '@' << i;
            for (int errPos : errOffs.second[i])
                ofsReads2 << "_" << errPos;
            ofsReads2 << "_CHR" << offsets.second[i].second << "_" << offsets.second[i].first << "\n";
            ofsReads2 << pairedReads.second[i] << "\n";
            ofsReads2 << '+' << i << "\n";
            // produce dummy quality scores
            for (size_t j = 0; j < readLen; ++j)
                ofsReads2 << "%";
            ofsReads2 << "\n";
        }
        ofsReads1.close();
        ofsReads2.close();

    } else {

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
            for (size_t j = 0; j < readLen; ++j)
                ofsReads << "%";
            ofsReads << "\n";
        }
        ofsReads.close();
        offsets.clear();
        errOffs.clear();
    }



    std::ofstream ofsReads(std::string(argv[2]) + "_cpginfo_fwd.tsv");
    ofsReads << "Chromosome\tPosition\tUnmethylated\tMethylated\tSampleRate\n";
    for (auto& cpg : synthGen.cpgMethRateFwd)
    {
        const uint64_t offset = cpg.first & 0x00000000ffffffffULL;
        ofsReads << (cpg.first >> 32) << "\t" << offset << "\t" << cpg.second.unmethCount << "\t" << cpg.second.methCount << "\t" << cpg.second.sampleRate << "\n";
    }
    ofsReads.close();
    ofsReads.open(std::string(argv[2]) + "_cpginfo_rev.tsv");
    ofsReads << "Chromosome\tPosition\tUnmethylated\tMethylated\tSampleRate\n";
    for (auto& cpg : synthGen.cpgMethRateRev)
    {
        const uint64_t offset = cpg.first & 0x00000000ffffffffULL;
        ofsReads << (cpg.first >> 32) << "\t" << offset << "\t" << cpg.second.unmethCount << "\t" << cpg.second.methCount << "\t" << cpg.second.sampleRate << "\n";
    }
    ofsReads.close();
    return 0;

}
