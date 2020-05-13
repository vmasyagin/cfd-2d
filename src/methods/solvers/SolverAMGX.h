#pragma once
#include "MatrixSolver.h"
#include "amgx_c.h"
#include <map>
#include <string>
#include <algorithm>

typedef std::map<int, std::vector<double>> block_row;

class SolverAMGX :
	public MatrixSolver
{
public:
	~SolverAMGX();
	virtual void init(Grid* g, int matrDimension, int blockDimension);

	virtual void zero();
    void	assemble();

	virtual void setMatrElement(int i, int j, double** matrDim);
	virtual void setRightElement(int i, double* vectDim);
	virtual void addMatrElement(int i, int j, double** matrDim);
	virtual void addRightElement(int i, double* vectDim);
	virtual void createMatrElement(int i, int j) {}
	virtual void setParameter(const char* name, int val);
	virtual void printToFile(const char* fileName);

protected:
	void initMatrVectors();

protected:
    AMGX_config_handle config;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle solver;
    AMGX_matrix_handle A;
    AMGX_vector_handle rhs;
    AMGX_vector_handle sln;

    int n, nnz, block_dimx, block_dimy, block_size, num_neighbors, matr_dim;

    int *row_ptrs = NULL, *col_indices = NULL, *neighbors = NULL;
    double *values = NULL, *diag = NULL, *dh_x = NULL, *dh_b = NULL;
    int *h_row_ptrs = NULL, *h_col_indices = NULL;
    double *h_values = NULL, *h_diag = NULL, *h_x = NULL, *h_b = NULL;

    block_row		*rows;

	char* solver_name;

	int PRINT_LEVEL = 0;
};

