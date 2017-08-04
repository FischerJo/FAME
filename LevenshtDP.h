#ifndef LEVENSHTDP_H
#define LEVENSHTDP_H

#include <string>
#include <limits>   // numeric limits (max)
#include <iostream> // cerr


#include "BandedMatrix.h"


// types of errors allowed
enum ERROR_T {MATCH, MISMATCH, INSERTION, DELETION};

// class for dynamic programming approach for computing Levenshtein
// distance and the corresponding alignment for two strings
//
// template parameter is the size type T of the underlying matrix
// i.e. should be the minimal size type that fits the max of the two string sizes
// and the bandwidth of the banded alignment
// i.e. the number of errors allowed (hence we only fill some diagonals in the matrix)
template <typename T, size_t band>
class LevenshtDP
{

    public:

        // ------ Ctors ------

        LevenshtDP() = delete;

        LevenshtDP(std::string& rowStr, const char* colStr);

        // -------------------


        // compute the levenshtein distance
        void runDPFill();

        // return edit distance
        // undefined behaviour if runDPFill was not run beforehand
        T getEditDist()
        {
            T minimum = std::numeric_limits<T>::max();
            for (long offset = 0; offset <= static_cast<long>(band); ++offset)
                minimum = std::min(minimum, dpMatrix(rowPat.size(), rowPat.size() + offset));
            return minimum;
        }

        // backtrack the result to obtain alignment
        //
        // RETURN:
        //      vector of size of rowPat containing a value of type ERROR_T
        //      the value represents the type of transition made from i-th to
        //      i+1 th character of rowPat to match with colPat
        std::vector<ERROR_T> backtrackDP();


    private:

        // recursive formulation of the levenshtein distance
        // DP drops in here
        //
        // Formula:
        //
        //          LevRec(i,j) = min {LevRec(i-1,j) + 1, LevRec(i,j-1) + 1,
        //                              LevRec(i-1,j-1) + isEqual(rowPat(i),colPat(j)) }
        //
        // ARGUMENTS:
        //          indices of positions in strings
        //          i       position in rowPat
        //          j       position in colPat
        //
        // RETURN:
        //          solution to frecurrence
        inline T LevRec(long i, long j)
        {
            // TODO test printout
            T matchFlag = rowPat[i-1] == colPat[j-1] ? static_cast<T>(0) : static_cast<T>(1);
            return std::min({dpMatrix(i-1,j) + 1, dpMatrix(i,j-1) + 1, dpMatrix(i-1,j-1) + matchFlag});
        }

        // the two strings that are compared
        // rowPat is represented by rows of the DP matrix
        // colPat is represented by columns in DP matrix
        std::string& rowPat;
        const char* colPat;

        // underlying dp matrix
        // note that bandwidth is extended by one for the edge cases
        BandedMatrix<T, band + 1> dpMatrix;

};


template <typename T, size_t band>
LevenshtDP<T, band>::LevenshtDP(std::string& rowStr, const char* colStr) :
        rowPat(rowStr)
    ,   colPat(colStr)
    ,   dpMatrix(rowStr.size() + 1, rowStr.size() + 1 + band, static_cast<T>(0))
{
    // check if template param is correct size type
    if (!std::is_integral<T>::value)
    {
        std::cerr << "\n\n(LevenshtDP) Template parameter is not a valid integral type! Terminating...\n\n";
        exit(1);
    }
}

template <typename T, size_t band>
void LevenshtDP<T, band>::runDPFill()
{

    // INIT BORDERS

    dpMatrix(0,0) = 0;

    // initialize first row
    for (long col = 1; col <= static_cast<long>(band); ++col)
    {
        dpMatrix(0, col) = col;
    }
    // init first column
    for (long row = 1; row <= static_cast<long>(band); ++row)
    {
        dpMatrix(row, 0) = row;
    }
    // init outermost lefthanded band
    for (long col = 0; col < static_cast<long>(rowPat.size() - band); ++col)
    {
        dpMatrix(col + (band + 1), col) = std::numeric_limits<T>::max() - col - band - 1;
    }
    // init outermost righthanded band
    for (long row = 0; row < static_cast<long>(rowPat.size()); ++row)
    {
        dpMatrix(row, row + (band + 1)) = std::numeric_limits<T>::max() - row - band - 1;
    }


    // FILL MATRIX
    for (long row = 1; row <= static_cast<long>(rowPat.size()); ++row)
    {

        for (long offset = -band; offset <= static_cast<long>(band); ++offset)
        {
            // skip borders
            if (row + offset <= 0)
                continue;

            dpMatrix(row, row + offset) = LevRec(row, row + offset);
        }
    }

    // TEST PRINTOUT
    // for (long row = 0; row <= static_cast<long>(rowPat.size()); ++row)
    // {
    //     for (long offset = -band - 1; offset <= static_cast<long>(band + 1); ++offset)
    //     {
    //         // skip borders
    //         if (row + offset < 0)
    //             continue;
    //         if (row  + offset > rowPat.size() + band)
    //             continue;
    //         std::cout << dpMatrix(row, row + offset) << "\t";
    //     }
    //     std::cout << "\n";
    // }
}

#endif /* LEVENSHTDP_H */
