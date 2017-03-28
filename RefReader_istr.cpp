

#include "RefReader_istr.h"

using namespace std;


void readReference(const char* const filename, vector<struct CpG>& cpgTab, vector<const char*>& genSeq, vector<std::size_t> genSeqLen)
{

    string line;

    // stores the length of the sequence read so far
    // will be resetted whenever a new chromosome is read
    unsigned int seqLen;

    // stores the sequence of each chromosome
    string seq;
    seq.reserve(MyConst::CHROMMAX);

    // stores the chromosome index we are currently reading
    uint8_t chrIndex = 0;

    // flag stating if we currently read a real chromosome assembly
    bool contFlag = false;
    // flag stating if the last character read in the previous line was a c
    bool lastC = false;

    while (getline(ifs, line)) {

        // Test for id tag line
        if (line.begin() != line.end() && *(line.begin()) == '>' && (line.begin() + 1) != line.end())
        {

            // if we read primary assembly previously, write it to vectors
            if (contFlag)
            {
                // TODO
                genSeq.push_back(seq.shrink_to_fit().c_str());
                genSeqLen.push_back(seqLen);
                // reset buffers
                // TODO FIX THIS
                seq = "";
                seq.reserve(MyConst::CHROMMAX);
                seqLen = 0;

            }
            // throw out unlocalized contigs
            if (*(++(line.begin())) != 'C')
            {
                contFlag = false;
                continue;

            } else {

                ++chrIndex;
                contFlag = true;
                continue;
            }

        }

        // check if we are in real primary chromosome assembly
        if (!contFlag)
        {
            continue;
        }

        readLine(line, lastC,

    }


    ofs << "Number of CpGs: " << cpg_num << "\n\n";

    ofs << "Number of Cs\tFrequency\n";

    for (int i = 0; i <= cpg_adj_c.size(); ++i)
    {

        ofs << i << "\t" << cpg_adj_c[i] << "\n";
    }

    ofs << "\nNumber of CpGs\tFrequency\n";
    for (int i = 0; i <= cpg_adj_cpg.size(); ++i)
    {

        ofs << i << "\t" << cpg_adj_cpg[i] << "\n";
    }

}
