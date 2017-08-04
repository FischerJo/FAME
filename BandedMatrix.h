#ifndef BANDEDMATRIX_H
#define BANDEDMATRIX_H

#include <vector>

// (2D)Matrix class for values of type T
// only storing the CENTRAL band which is of size band
// in both directions away from the main diagonal
//
// Example: band = 1, 5x5 matrix, T=uint8_t
//
//    bandwidth
//     <-->
//   | 0  1  -  -  - |
//   | 1  0  1  -  - |
//   | -  1  0  0  - |
//   | -  -  0  2  3 |
//   | -  -  -  3  4 |
//
//   - is not represented by matrix
//
template <typename T, size_t band>
class BandedMatrix
{

    public:

        // ------ Ctors ------

        BandedMatrix() = delete;

        // Create BandedMatrix object with nrow rows and ncol columns
        // Second version fills matrix with copies of initVal
        BandedMatrix(size_t nrow, size_t ncol);
        BandedMatrix(size_t nrow, size_t ncol, T initVal);

        // -------------------


        // access element in row rowIdx and using the bandwidth as offset
        // start indexing with 0
        //
        // ARGUMENTS:
        //      rowIdx      row of desired object
        //      offset      offset of object in the row (offset \in [0, 2*band+1])
        //
        // RETURN:
        //      reference to desired object
        inline T& operator()(long rowIdx, long colIdx)
        {
            // in case we are at the borders
            // first term gives the offset for cut off elements to the left of matrix
            // second term gives offset due to the column
            long offset = std::max(static_cast<long>(band) - rowIdx, static_cast<long>(0)) + (colIdx - rowIdx + static_cast<long>(band));
            // compute number of elements that are cut
            return m[(rowIdx*(2*band + 1)) + offset];
        }


    private:

        // dimensions of the matrix (second dimension is implicit by bandwidth)
        size_t rowNum, colNum;
        // underlying matrix object
        // banded rows are stored consecutively
        std::vector<T> m;

};


template <typename T, size_t band>
BandedMatrix<T,band>::BandedMatrix(size_t nrow, size_t ncol) :
        rowNum(nrow)
    ,   colNum(ncol)
        // note +1 for main diagonal
    ,   m(nrow*(2*band + 1))
{
}


template <typename T, size_t band>
BandedMatrix<T,band>::BandedMatrix(size_t nrow, size_t ncol, T initVal) :
        rowNum(nrow)
    ,   colNum(ncol)
        // note +1 for main diagonal
    ,   m(nrow*(2*band + 1), initVal)
{
}


#endif /* BANDEDMATRIX_H */
