//
// Created by victor on 11.03.2020.
//
#include <grid.h>

#ifndef CFD_2D_INITIALIZATOR_H
#define CFD_2D_INITIALIZATOR_H

#endif //CFD_2D_INITIALIZATOR_H

class Initializator{
private:
    static int n; // dimension of the matrix in terms of block-units
    static int nnz; // number of non-zero entries

public:
    void static init(Grid* g, int matr_dim, int block_dimx);
    static int get_n(){return n;};
    static int get_nnz(){return nnz;};
};
