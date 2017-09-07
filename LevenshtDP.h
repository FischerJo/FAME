//	Metal - A fast methylation alignment and calling tool for WGBS data.
//	Copyright (C) 2017  Jonas Fischer
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	Jonas Fischer	jonaspost@web.de

#ifndef LEVENSHTDP_H
#define LEVENSHTDP_H

#include <string>
#include <limits>   // numeric limits (max)
#include <iostream> // cerr


#include "BandedMatrix.h"


// types of errors allowed
enum ERROR_T {MATCHING, MISMATCH, INSERTION, DELETION};

// class for dynamic programming approach for computing Levenshtein
// distance and the corresponding alignment for two strings
// it is REQUIRED that the first given string (rowStr in Ctor) must be fully
// matched!
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

        LevenshtDP(const std::string& rowStr, const char* colStr);

        // -------------------


        // compute the levenshtein distance
        //
        // ARGUMENTS:
        //      comp        comparison function for the letters
        //
        template <typename C>
        void runDPFill(C& comp);
        // for querying from right to left
        template <typename C>
        void runDPFillRev(C& comp);

        // return edit distance
        // undefined behaviour if runDPFill was not run beforehand
        T getEditDist()
        {
            T minimum = std::numeric_limits<T>::max();
            for (long offset = -static_cast<long>(band); offset <= static_cast<long>(band); ++offset)
                minimum = std::min(minimum, dpMatrix(rowPat.size(), rowPat.size() + offset));
            return minimum;
        }

        // backtrack the result to obtain alignment
        //
        // ARGUMENTS:
        //      comp        comparison function for the letters
        //      align       gives the alignment using error types as flags SHOULD BE EMPTY ON CALL
        //
        // MODIFICATIONS:
        //      alignment will bevector of size of rowPat containing a value of type ERROR_T
        //      the value represents the type of transition made from i-th to
        //      i+1 th character of rowPat to match with colPat
        template <typename C>
        void backtrackDP(C& comp, std::vector<ERROR_T>& align);
        template <typename C>
        void backtrackDPRev(C& comp, std::vector<ERROR_T>& align);


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
        //          comp    comparison function for the letters
        //
        //
        // RETURN:
        //          solution to recurrence
        template <typename C>
        inline T LevRec(long i, long j, C& comp)
        {
            T mismatchTest = comp(rowPat[rowPat.size() - i], *(colPat-j+1));
            return std::min({dpMatrix(i-1,j) + 1, dpMatrix(i,j-1) + 1, dpMatrix(i-1,j-1) + mismatchTest});
        }
        // if we want to align reverse sequence
        template <typename C>
        inline T LevRecRev(long i, long j, C& comp)
        {
            T mismatchTest = comp(rowPat[i-1], *(colPat-j+1));
            return std::min({dpMatrix(i-1,j) + 1, dpMatrix(i,j-1) + 1, dpMatrix(i-1,j-1) + mismatchTest});
        }

        // the two strings that are compared
        // rowPat is represented by rows of the DP matrix
        // colPat is represented by columns in DP matrix
        const std::string& rowPat;
        const char* colPat;

        // underlying dp matrix
        // note that bandwidth is extended by one for the edge cases
        BandedMatrix<T, band + 1> dpMatrix;

};


template <typename T, size_t band>
LevenshtDP<T, band>::LevenshtDP(const std::string& rowStr, const char* colStr) :
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
template <typename C>
void LevenshtDP<T, band>::runDPFill(C& comp)
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

            dpMatrix(row, row + offset) = LevRec<C>(row, row + offset, comp);
        }
    }

    // // TEST PRINTOUT
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
template <typename T, size_t band>
template <typename C>
void LevenshtDP<T, band>::runDPFillRev(C& comp)
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

            dpMatrix(row, row + offset) = LevRecRev<C>(row, row + offset, comp);
        }
    }
    // // TEST PRINTOUT
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



template <typename T, size_t band>
template <typename C>
void LevenshtDP<T, band>::backtrackDP(C& comp, std::vector<ERROR_T>& alignment)
{

    // the trace of errors
    // note that insertion means we have an additional character in the PATTERN (rowPat)
    // and deletion analogously
    // we start indexing for the character in the back
    // that is alignment[0] is the error type of the last character
    alignment.reserve(rowPat.size() + band);
    // find minimum in the last part of table
    T minimum = std::numeric_limits<T>::max();
    // stores the index of the cell of the minimum
    size_t col = 0;
    for (long offset = -static_cast<long>(band); offset <= static_cast<long>(band); ++offset)
    {
        const T& valRef = dpMatrix(rowPat.size(), rowPat.size() + offset);
        if (valRef < minimum)
        {

            col = rowPat.size() + offset;
            minimum = valRef;
        }
    }

    // BACKTRACK
    // starting from the minimum cell in the last column
    for (long row = rowPat.size(); row > 0;)
    {

        // check if characters match
        T mismatchFlag = comp(rowPat[rowPat.size() - row], *(colPat-col+1));

        // see where does it came from
        if (dpMatrix(row-1,col) < dpMatrix(row,col-1))
        {
            // we have an insertion in the pattern
            if (dpMatrix(row-1,col) + 1 < dpMatrix(row-1,col-1) + mismatchFlag)
            {
                alignment.emplace_back(INSERTION);
                --row;

            } else {

                // test if mismatch
                if (mismatchFlag)
                    alignment.emplace_back(MISMATCH);
                else
                    alignment.emplace_back(MATCHING);
                --row;
                --col;

            }

        } else {

            // we have a deletion in the pattern
            if (dpMatrix(row,col-1) + 1 < dpMatrix(row-1,col-1) + mismatchFlag)
            {
                alignment.emplace_back(DELETION);
                --col;

            } else {

                // test if mismatch
                if (mismatchFlag)
                    alignment.emplace_back(MISMATCH);
                else
                    alignment.emplace_back(MATCHING);
                --row;
                --col;

            }
        }

    }
}
template <typename T, size_t band>
template <typename C>
void LevenshtDP<T, band>::backtrackDPRev(C& comp, std::vector<ERROR_T>& alignment)
{

    // the trace of errors
    // note that insertion means we have an additional character in the PATTERN (rowPat)
    // and deletion analogously
    // we start indexing for the character in the back
    // that is alignment[0] is the error type of the last character
    alignment.reserve(rowPat.size() + band);
    // find minimum in the last part of table
    T minimum = std::numeric_limits<T>::max();
    // stores the index of the cell of the minimum
    size_t col = 0;
    for (long offset = -static_cast<long>(band); offset <= static_cast<long>(band); ++offset)
    {
        const T& valRef = dpMatrix(rowPat.size(), rowPat.size() + offset);
        if (valRef < minimum)
        {

            col = rowPat.size() + offset;
            minimum = valRef;
        }
    }

    // BACKTRACK
    // starting from the minimum cell in the last column
    for (long row = rowPat.size(); row > 0;)
    {

        // check if characters match
        T mismatchFlag = comp(rowPat[row - 1], *(colPat -col + 1));

        // see where does it came from
        if (dpMatrix(row-1,col) < dpMatrix(row,col-1))
        {
            // we have an insertion in the pattern
            if (dpMatrix(row-1,col) + 1 < dpMatrix(row-1,col-1) + mismatchFlag)
            {
                alignment.emplace_back(INSERTION);
                --row;

            } else {

                // test if mismatch
                if (mismatchFlag)
                    alignment.emplace_back(MISMATCH);
                else
                    alignment.emplace_back(MATCHING);
                --row;
                --col;

            }

        } else {

            // we have a deletion in the pattern
            if (dpMatrix(row,col-1) + 1 < dpMatrix(row-1,col-1) + mismatchFlag)
            {
                alignment.emplace_back(DELETION);
                --col;

            } else {

                // test if mismatch
                if (mismatchFlag)
                    alignment.emplace_back(MISMATCH);
                else
                    alignment.emplace_back(MATCHING);
                --row;
                --col;

            }
        }
    }
}

#endif /* LEVENSHTDP_H */
