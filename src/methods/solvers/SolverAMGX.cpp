#include "SolverAMGX.h"
#include "config.h"


void SolverAMGX::initMatrVectors()
{
    AMGX_initialize();
    AMGX_initialize_plugins();
    std::string config_name = CONFIGS_DIR;
    config_name.append(solver_name + 5);
    config_name += ".json";
    AMGX_config_create_from_file(&config, config_name.c_str());
    AMGX_resources_create_simple(&rsrc, config);

    AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, config);
    AMGX_matrix_create(&A, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&sln, rsrc, AMGX_mode_dDDI);
}

void SolverAMGX::zero() {
    AMGX_matrix_destroy(A);
    AMGX_vector_destroy(rhs);
    AMGX_vector_destroy(sln);

    AMGX_matrix_create(&A, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&sln, rsrc, AMGX_mode_dDDI);

    memset(h_values, 0, sizeof(double)*nnz*block_size);
    memset(h_b, 0, sizeof(double)*n*block_dimx);
    memset(h_x, 0, sizeof(double)*n*block_dimx);
    memset(h_col_indices, 0, nnz * sizeof(int));
    memset(h_row_ptrs, 0, (n + 1)*sizeof(int));
    for (int i = 0; i < n; i++) {
        block_row & r = rows[i];
        for (block_row::iterator it = r.begin(); it != r.end(); it++) {
            std::fill(it->second.begin(), it->second.end(), 0);
        }
    }
}

SolverAMGX::~SolverAMGX()
{
    AMGX_matrix_destroy(A);
    AMGX_vector_destroy(rhs);
    AMGX_vector_destroy(sln);
    AMGX_solver_destroy(solver);
    AMGX_resources_destroy(rsrc);

    AMGX_finalize_plugins();
    AMGX_finalize();
	delete[] x;
}

void SolverAMGX::setMatrElement(int i, int j, double** matrDim)
{
    // TODO: реализовать позже при необходимости
}

void SolverAMGX::addMatrElement(int i, int j, double** matrDim)
{
    int block_num = (int) (matr_dim / block_dimx);
    std::vector<double> row;
    for (int i_block = 0; i_block < block_num; i_block++){
        for (int j_block = 0; j_block < block_num; j_block++){
            row = rows[i * block_num + i_block][j * block_num + j_block];
            if (row.empty()) {
                //row.assign(block_size, 0);
                row.resize(block_size, 0);
            }
            for (int ii = 0; ii < block_dimx; ii++) {
                for (int jj = 0; jj < block_dimy; jj++) {
                    row[ii*block_dimx + jj] += matrDim[i_block*block_dimx + ii][j_block*block_dimy + jj];
                }
            }
            rows[i * block_num + i_block][j * block_num + j_block] = row;
        }
    }
}


void SolverAMGX::setRightElement(int i, double* vectDim)
{
    memcpy(&h_b[i*matr_dim], vectDim, matr_dim * sizeof(double));
}


void SolverAMGX::addRightElement(int i, double* vectDim)
{
    for (int ii = 0; ii < matr_dim; ii++) {
        h_b[ii + i * matr_dim] += vectDim[ii];
    }
}

void SolverAMGX::setParameter(const char* name, int val)
{
	if (strcmp(name, "PRINT_LEVEL") == 0) {
		PRINT_LEVEL = val;
	}
}

void SolverAMGX::assemble() {
    int ind = 0;
    for (int i = 0; i < n; i++) {
        h_row_ptrs[i] = ind;
        ind += rows[i].size();
    }
    h_row_ptrs[n] = ind;

    ind = 0;
    for (int i = 0; i < n; i++) {
        block_row & r = rows[i];
        for (block_row::iterator it = r.begin(); it != r.end(); it++) {
            h_col_indices[ind] = it->first;
            double* array = &it->second[0];
            memcpy(&h_values[ind*block_size], array, block_size*sizeof(double));
            ind++;
        }
    }
}

void SolverAMGX::init(Grid* g, int matrDimension, int blockDimension) {
    Initializator::init(g, matrDimension, blockDimension);
    blockDim = blockDimension;
    block_dimx = blockDimension;
    block_dimy = blockDimension;
    block_size = block_dimx * block_dimy;
    matr_dim = matrDimension;
    n = Initializator::get_n();
    nnz = Initializator::get_nnz();

    h_col_indices = new int[nnz];
    h_row_ptrs = new int[n + 1];
    h_values	= new double[nnz * block_size];
    h_b = new double[n * block_dimx];
    h_x = new double[n * block_dimx];

    x = new double[n*block_dimx];

    rows = new block_row[n];

    initMatrVectors();
}

void SolverAMGX::printToFile(const char* fileName)
{
    FILE * fp = fopen(fileName, "w");
	for (int i = 0; i < n*block_dimx; i++) {
	    fprintf(fp, "%25.16e  ", values[h_row_ptrs[i]]);
    }
	fprintf(fp, "\n\n=============================================================================================\n\n\n");
	for (int i = 0; i < n*block_dimx; i++) {
		fprintf(fp, "%25.16e  ", x[i]);
	}

	fclose(fp);
}





