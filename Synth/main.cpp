
#include <fstream>

#include "SynthDS.h"

// --------------- PARAMS ---------------

constexpr size_t refLen = 1000000000;
constexpr double mthRate = 0.5;
constexpr size_t readLen = 100;
constexpr size_t readNum = 10000000;
constexpr unsigned int errNum = 2;


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
    std::vector<std::string> fwdReads = synthGen.genReadsFwdRef(readLen, readNum, errNum, offsets);


    // std::ofstream ofsRef("genRef_1err.fasta");
    // ofsRef << ">CHRSYNTH1\n";
    // ofsRef << synthGen.getRef();
    // ofsRef.close();

    std::ofstream ofsReads(std::string(argv[1]) + "_fwd.fastq");

    for (size_t i = 0; i < fwdReads.size(); ++i)
    {
        // generate fastq format of reads
        ofsReads << '@' << i << "_CHR" << offsets[i].second << "_" << offsets[i].first << "\n";
        ofsReads << fwdReads[i] << "\n";
        ofsReads << '+' << i << "\n";
        // produce dummy quality scores
        for (size_t i = 0; i < readLen; ++i)
            ofsReads << "%";
        ofsReads << "\n";
    }
    ofsReads.close();
    offsets.clear();
    std::vector<std::string> revReads = synthGen.genReadsRevRef(readLen, readNum, errNum, offsets);
    ofsReads.open(std::string(argv[1] + "_rev.fastq");
    for (size_t i = 0; i < fwdReads.size(); ++i)
    {
        // generate fastq format of reads
        ofsReads << '@' << i << "_CHR" << offsets[i].second << "_" << offsets[i].first << "\n";
        ofsReads << fwdReads[i] << "\n";
        ofsReads << '+' << i << "\n";
        // produce dummy quality scores
        for (size_t i = 0; i < readLen; ++i)
            ofsReads << "%";
        ofsReads << "\n";
    }
    ofsReads.close();
    return 0;

}
