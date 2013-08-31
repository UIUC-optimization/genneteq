/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: matrix.h                                                            */
/* Description:                                                              */
/*   Contains design and implementation details for a template matrix class. */
/*****************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

// Required include's
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <vector>
using std::vector;

template <typename T>
class Matrix {
  public:
    Matrix()
    {
        // Do nothing
    }

    ~Matrix()
    {
        // Do nothing
    }

    void initialize(int size, T val)
    {
        initialize(size, size, val);
        return;
    }

    void initialize(int rows, int cols, T val)
    {
        numRows = rows;
        numCols = cols;
        data.resize(numRows);
        for (int i = 0; i < numRows; ++i) {
            data[i].resize(numCols);
            for (int j = 0; j < numCols; ++j) {
                data[i][j] = val;
            }
        }
        return;
    }

    inline
    T get(int i, int j) const
    {
        return data[i][j];
    }

    inline
    void set(int i, int j, T val)
    {
        data[i][j] = val;
        return;
    }

    inline
    void increment(int i, int j, T val)
    {
        data[i][j] += val;
        return;
    }

    inline 
    void decrement(int i, int j, T val)
    {
        data[i][j] -= val;
        return;
    }

    void print() const
    {
        std::cout << "Matrix with dimensions " << numRows << " x " << numCols 
                  << std::endl;;
        for (int i = 0; i < numRows; ++i) {
            std::cout << "|";
            for (int j = 0; j < numCols; ++j) {
                std::cout << data[i][j] << ",";
            }
            std::cout << "|" << std::endl;
        }
        return;
    }

  protected:
    // Variables
    int numRows;
    int numCols;
    vector<vector<T> > data;

  private:
    // Nothing
};

#endif // MATRIX_H

