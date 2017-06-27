#ifndef SYNTHDS_H
#define SYNTHDS_H


#include <random> // mersenne twister
#include <string>
#include <vector>
#include <array>


#define CORENUM  16

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

        // ----------


        // generate set of reads of given length, drawn from the reference or reverse complement
        // and with errors introduced
        //
        // ARGUMENTS:
        //          readLen     length of the returned reads
        //          readNum     number of reads to be generated
        //          maxErrNum   maximum number of errors that can be introduced per read
        //
        // RETURN:
        //          vector of length readNum holding the generated reads
        std::vector<std::string> genReadsFwdFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum);
        std::vector<std::string> genReadsRevFixed(const size_t readLen, const size_t readNum, const unsigned int maxErrNum);

        // generate set of reads with length in between specified boundaries, drawn from the reference
        // and with errors introduced
        //
        // ARGUMENTS:
        //          readLenMin  minimum length of the returned reads
        //          readLenMax  maximum length of the returned reads
        //          readNum     number of reads to be generated
        //          maxErrNum   maximum number of errors that can be introduced per read
        //
        // RETURN:
        //          vector of length readNum holding the generated reads
        // std::vector<std::string> genReadsRange(const size_t readLenMin, const size_t readLenMax, const size_t readNum, const unsigned int maxErrNum);

        inline std::string& getRef() { return refFwd;}

    private:

        // generate the reference sequence using a rng initialized with the given seed
        void initReference(const size_t refLen, const unsigned int seed);

        // reference sequence forward and reverse strand
        std::string refFwd;
        std::string refRev;

        // set of random number generators
        std::array<std::mt19937, CORENUM> randGen;

        // alphabet mapping from random numbers to letters
        std::uniform_int_distribution<int> toIndex;
        const std::array<char, 4> alphabet;

};

#endif /* SYNTHDS_H */
