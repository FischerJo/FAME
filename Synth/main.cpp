
#include <fstream>

#include "SynthDS.h"


constexpr size_t refLen = 1000000000;


// --------------- MAIN -----------------
//
//
int main(int argc, char** argv)
{

    SynthDS synthGen(refLen);

    std::vector<size_t> offsets;
    std::vector<std::string> fwdReads = synthGen.genReadsFwdFixed(100, 100000, 1, offsets);

    std::vector<std::string> revReads = synthGen.genReadsRevFixed(100, 100000, 1);


    std::ofstream ofsRef("genRef_1err.fasta");
    ofsRef << ">CHRSYNTH1\n";
    ofsRef << synthGen.getRef();
    ofsRef.close();

    std::ofstream ofsReads("genReadsFwd_1err.fastq");

    for (size_t i = 0; i < fwdReads.size(); ++i)
    {
        // generate fastq format of reads
        ofsReads << '@' << i << "_" << offsets[i] << std::endl;
        ofsReads << fwdReads[i] << std::endl;
        ofsReads << '+' << i << std::endl;
        ofsReads << 'N' << std::endl;
    }
    ofsReads.close();
    ofsReads.open("genReadsRev_1err.fastq");
    for (size_t i = 0; i < revReads.size(); ++i)
    {
        const size_t id = i + fwdReads.size();
        // generate fastq format of reads
        ofsReads << '@' << id << std::endl;
        ofsReads << revReads[i] << std::endl;
        ofsReads << '+' << id << std::endl;
        ofsReads << 'N' << std::endl;
    }

    ofsReads.close();
    return 0;

}
