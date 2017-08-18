#ifndef SYNTHDS_H
#define SYNTHDS_H


#include <random> // mersenne twister
#include <string>
#include <vector>
#include <array>
#include <unordered_map>


// --------------- PARAMS ---------------

constexpr size_t refLen = 1000000000;
constexpr double mthRate = 0.6;
constexpr size_t readLen = 100;
constexpr size_t readNum = 20000000;
constexpr unsigned int errNum = 2;
#define CORENUM  1

// class representing a synthetic dataset
class SynthDS
{

    public:

        // Ctor -----

        SynthDS() = delete;

        // generates a random reference sequence of the specified length
        //
        // ARGUMENTS:
        //          refLen      length of the reference sequence to be produced
        SynthDS(const size_t refLen);

        // generates a random reference sequence of the specified length
        //
        // ARGUMENTS:
        //          refLen      length of the reference sequence to be produced
        //          seed        initial seed for the internal pseudo random number generator
        SynthDS(const size_t refLen, const unsigned int seed);

        // loads reference specified by file into DS
        //
        // ARGUMENTS:
        //          genFile     file containing reference in fasta format
        //          methRate    specified methylation rate [0, 1.0] that should be reflected in generated reads
        //          convRate    specified bisulfite conversion rate [0, 1.0]
        SynthDS(const char* genFile, const double methRate, const double convRate = 1.0);

        // ----------


        // generate set of reads of given length, drawn from the reference or reverse complement
        // and with errors introduced
        //
        // ARGUMENTS:
        //          readLen     length of the returned reads
        //          readNum     number of reads to be generated
        //          maxErrNum   maximum number of errors that can be introduced per read
        //          offsets     will hold the offsets where the read stems from
        //
        // RETURN:
        //          vector of length readNum holding the generated reads
        std::vector<std::string> genReadsFwdFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<size_t>& offsets);
        std::vector<std::string> genReadsRevFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum);

        // generate set of reads of given length drawn from a LOADED reference
        std::vector<std::string> genReadsFwdRef(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<std::pair<size_t, size_t> >& offsets, std::vector<std::array<int, errNum> >& errOffs);
        std::vector<std::string> genReadsRevRef(const size_t readLen, const size_t readNum, const unsigned int maxErrNum, std::vector<std::pair<size_t, size_t> >& offsets, std::vector<std::array<int, errNum> >& errOffs);

        inline std::string& getRef() { return refFwd;}

        // storing methylation rates the same way we do in Metal
        // key is 32bit most significant chromosome; 32bit least significant position in chromosome
        // value is first: unmethylated second: methylated count
        std::unordered_map<uint64_t, std::pair<uint16_t, uint16_t> > cpgMethRateFwd;
        std::unordered_map<uint64_t, std::pair<uint16_t, uint16_t> > cpgMethRateRev;


    // private:

        // generate the reference sequence using a rng initialized with the given seed
        void initReference(const size_t refLen, const unsigned int seed);

        // read the reference from file
        // discards all nonlocalized contig stuff, i.e. only read reference starting with >C
        void loadRefSeq(const char* genFile);

        // reference sequence forward and reverse strand
        std::string refFwd;
        std::string refRev;
        // set of chrosomal sequences for given reference sequence
        std::vector<std::string> refSeqFwd;
        std::vector<std::string> refSeqRev;

        // set of random number generators
        std::array<std::mt19937, CORENUM> randGen;

        // alphabet mapping from random numbers to letters
        std::uniform_int_distribution<int> toIndex;
        std::bernoulli_distribution methToss;
        std::bernoulli_distribution convToss;
        const std::array<char, 4> alphabet;

};

#endif /* SYNTHDS_H */
