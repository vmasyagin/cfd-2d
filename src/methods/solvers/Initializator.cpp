//
// Created by victor on 11.03.2020.
//

#include "Initializator.h"

int Initializator::n;
int Initializator::nnz;

void Initializator::init(Grid* g, int matr_dim, int block_dimx) {
    try {
        int block_num = matr_dim / block_dimx; // blocks per matrix dimension
        n = g->cCount * block_num;
        for (int i = 0; i < g->cCount; i++) {
            Cell& c = g->cells[i];
            nnz+=block_num * block_num;
            for (int j : c.neigh) {
                if (j > -1) nnz+=block_num * block_num;
            }
        }
    }
    catch (Exception e){
        printf(e.getMessage());
    }
}

