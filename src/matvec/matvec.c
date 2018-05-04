/*
 * 2017.04.23  - Created, based on ideas from https://github.com/sergeyvoronin/LowRankMatrixDecompositionCodes
 * 2017.05.22  - Total rejig to allow backend replace using preprocessor directive/condition compilation (e.g. __SYSMATH_MKL, etc)
 * 2017.06.03  - Introduced support for explicit set of system threads
 * 2017.06.06  - Introduced pseudo random number generator to remove mandatory dependency on Intel MKL
 * 2017.06.24  - Added basic sparse matrix support (CSR format)
 * 2018.04.15  - Improved Sparse matrix support (introduced SSS format)
 */

/* 
 * File:   matvec.c
 * Author: Ilamah, Osho
 * Created on 23 April 2017, 08:42 as matvec.c
 * 
 * Implementation of various matrix and vector support functions. Matrices are assumed to be ROW MAJOR.
 * 
 */
//TODO - csr_from_matrix should simply like sss_from_matrix

#include <stdio.h>


#include "../headers/external/matvec.h"
#include "../headers/external/xalloc.h"

#ifdef  __SYSMATH_OPENBLAS__
#include <omp.h>
#endif

// private function prototypes
uint32_t get_nnz_count(double* data, const uint64_t BUFFER_SIZE);

Vector* create_vector(uint32_t nrows) {
    Vector* v = malloc(sizeof (Vector));
    v->data = (double*) aligned_calloc(nrows, sizeof (double), 64);
    if (v->data == NULL) {
        return NULL;
    }
    v->nrows = nrows;
    return v;
}

/**
 * Creates a new matrix and and initializes all the elements to zero. If symmetric flag is set, only the upper half of the matrix is stored
 * @param nrows
 * @param ncols
 * @param symmetric
 * @return 
 */
Matrix* create_matrix(uint32_t nrows, uint32_t ncols, bool symmetric) {
    Matrix* m = malloc(sizeof (Matrix));
    if (symmetric && (nrows == ncols)) {
        //symmetric matrix  - for symmetric dense matrix, we store lower triangular matrix (including diagonal elements), this requires (N*(N+1))/2:
        //      Bottom half (excluding diagonal elements = [N*N - N]/2
        //      Diagonal elements = N
        //      Total Elements = N + [N*N - N]/2   ==  (N*(N+1))/2
        m->symmetric = symmetric;
        m->data = (double*) aligned_calloc((nrows * (nrows + 1)) / 2, sizeof (double), 64);
    } else {
        m->symmetric = false;
        m->data = (double*) aligned_calloc(nrows * ncols, sizeof (double), 64);
    }
    if (m->data == NULL) {
        return NULL;
    }
    m->nrows = nrows;
    m->ncols = ncols;
    return m;
}

/**
 * Creates a new CSRMatrix and initializes all the elements to zero. 
 * @param nrows
 * @param ncols
 * @param nnz
 * @return 
 */
CSRMatrix* create_csrmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz) {
    CSRMatrix* m = malloc(sizeof (CSRMatrix));
    m->values = (double*) aligned_calloc(nnz, sizeof (double), 64);
    m->col_ind = (uint32_t*) aligned_calloc(nnz, sizeof (uint32_t), 64);
    m->row_ptr = (uint32_t*) aligned_calloc(nrows + 1, sizeof (uint32_t), 64);
    if (m->values == NULL || m->col_ind == NULL || m->row_ptr == NULL) {
        return NULL;
    }
    m->nrows = nrows;
    m->ncols = ncols;
    m->nnz = nnz;
    return m;
}

SSSMatrix* create_sssmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnzl) {
    SSSMatrix* m = malloc(sizeof (SSSMatrix));
    //
    m->dvalues = (double*) aligned_calloc(nrows, sizeof (double), 64);
    m->values = (double*) aligned_calloc(nnzl, sizeof (double), 64);
    m->col_ind = (uint32_t*) aligned_calloc(nnzl, sizeof (uint32_t), 64);
    m->row_ptr = (uint32_t*) aligned_calloc(nrows + 1, sizeof (uint32_t), 64);
    if (m->values == NULL || m->values == NULL || m->col_ind == NULL || m->row_ptr == NULL) {
        return NULL;
    }
    //
    m->nrows = nrows; 
    m->ncols = ncols;
    m->nnzl = nnzl;
    return m;
}

//

Matrix* wrap_as_matrix(uint32_t nrows, uint32_t ncols, bool symmetric, double* data) {
    Matrix* m = malloc(sizeof (Matrix));
    m->data = data;
    m->nrows = nrows;
    m->ncols = ncols;
    m->symmetric = symmetric;
    return m;
}

CSRMatrix* wrap_as_csrmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz, double* data, uint32_t* row_ptr, uint32_t* col_ind) {
    CSRMatrix* m = malloc(sizeof (CSRMatrix));
    m->values = data;
    m->col_ind = col_ind;
    m->row_ptr = row_ptr;
    m->nrows = nrows;
    m->ncols = ncols;
    m->nnz = nnz;
    return m;
}

SSSMatrix* wrap_as_sssmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnzl, double* dvalues, double* values, uint32_t* row_ptr, uint32_t* col_ind) {
    SSSMatrix* m = malloc(sizeof (SSSMatrix));
    m->dvalues = dvalues;
    m->values = values;
    m->col_ind = col_ind;
    m->row_ptr = row_ptr;
    m->nrows = nrows;
    m->ncols = ncols;
    m->nnzl = nnzl;
    return m;
}
//

void delete_vector(Vector* vector) {
    if (vector != 0) {
        aligned_free(vector->data);
        free(vector);
        vector = 0;
    }
}

void delete_matrix(Matrix* matrix) {
    if (matrix != 0) {//prevents seg fault, assuming client set pointer to zer for unalloced
        aligned_free(matrix->data);
        free(matrix);
        matrix = 0;
    }
}

void delete_csrmatrix(CSRMatrix* matrix) {
    if (matrix != 0) {
        aligned_free(matrix->values);
        aligned_free(matrix->row_ptr);
        aligned_free(matrix->col_ind);
        free(matrix);
        matrix = 0;
    }
}

void delete_sssmatrix(SSSMatrix* matrix) {
    if (matrix != 0) {
        aligned_free(matrix->dvalues);
        aligned_free(matrix->values);
        aligned_free(matrix->row_ptr);
        aligned_free(matrix->col_ind);
        free(matrix);
        matrix = 0;
    }
}

/**
 * Sets the value of the provided matrix, at the specified row and column indices. It is assumed that the backing array is organized in ROW MAJOR format from left to right and then top to down.  If the 
 * matrix is symmetric, only a portion is stored - the upper triangular half (+ diagonal) in row major format, which corresponds to the lower triangular half (+ diagonal) in COLUMN major format.
 * 
 * @param M
 * @param row_num
 * @param col_num
 * @param val
 */
//https://stackoverflow.com/questions/3187957/how-to-store-a-symmetric-matrix
//https://www.codeguru.com/cpp/cpp/algorithms/general/article.php/c11211/TIP-Half-Size-Triangular-Matrix.htm
//
//QC'ed

void set_matrix_element(Matrix* M, const uint32_t row_num, const uint32_t col_num, const double val) {
    //#ifdef __DEBUG__
    //    printf("%d,%d:: %lf\n", row_num, col_num, val);
    //#endif
    if (M->symmetric) {
        if (row_num >= col_num) {
            M->data[((col_num * M->nrows) + row_num - ((col_num * (col_num + 1)) / 2))] = val; //LOWER HALF - transpose upper half ROW MAJOR
        } else {
            M->data[((row_num * M->ncols) + col_num - ((row_num * (row_num + 1)) / 2))] = val; //Upper half - normal write as ROW MAJOR  
        }
    } else {
        M->data[(row_num * M->ncols) + col_num] = val; //ROW MAJOR
    }
}

/**
 * Gets the value of the provided matrix, at the specified row and column numbers. Assumes row major format
 * @param M
 * @param row_num
 * @param col_num
 * @return 
 */
//QC'ed

double get_matrix_element(const Matrix* M, const uint32_t row_num, const uint32_t col_num) {
    if (M->symmetric) {
        if (row_num >= col_num) {
            return M->data[((col_num * M->nrows) + row_num - ((col_num * (col_num + 1)) / 2))]; //LOWER HALF - transpose upper half ROW MAJOR        
        } else {
            return M->data[((row_num * M->ncols) + col_num - ((row_num * (row_num + 1)) / 2))]; //Upper half - normal read as ROW MAJOR  
        }
    } else {
        return M->data[(row_num * M->ncols) + col_num];
    }
}

//unpacks a truncatedly (half triangle) stored dense symmetric matrix into an untruncated storage backing

Matrix* unpack_symmetric_matrix(Matrix* M) {
    if (M->symmetric) {
        Matrix* A = create_matrix(M->nrows, M->ncols, 0);
        if (A == NULL) {
            return NULL;
        }
        const uint32_t nrows = M->nrows, ncols = M->ncols;
        for (uint32_t i = 0, j = 0; i < nrows; i++) {
            //read diagonal
            j = i;
            double value = M->data[((i * ncols) + j - ((i * (i + 1)) / 2))]; //copied from get_matrix_element
            A->data[(i * ncols) + j] = value; //copied from set matrix
            for (uint32_t j = i + 1; j < ncols; j++) {
                //set upper
                value = M->data[((i * ncols) + j - ((i * (i + 1)) / 2))];
                A->data[(i * ncols) + j] = value;
                //set lower
                //uint32_t ii = j, jj = i; 
                A->data[(j * ncols) + i] = value;
            }
        }
        return A;
    } else {
        return M;
    }
}


//only support non symetric matrix - for symetric, use SSS sparse matrix plse
//QC'ed

CSRMatrix* csr_from_matrix(const Matrix* M) {
    if (M->symmetric) {
        fprintf(stderr, "MATVEC: SSS Matrix should be used for symmetric cases");
        return NULL;
    }
    const double ZERO = 0;
    const uint64_t BUFFER_SIZE = M->ncols * M->nrows;
    const uint32_t nnz = get_nnz_count(M->data, BUFFER_SIZE);
    const uint32_t ncols = M->ncols, nrows = M->nrows;
    //now we know the number of elements, allocate the sparse CSR matrix
    CSRMatrix* out = create_csrmatrix(M->nrows, M->ncols, nnz);
    if (out == NULL) {
        return NULL;
    }
    //
    uint32_t row = 0, row_ptr = 0, sparse_element_idx = -1, elements_in_row = 0;
    for (uint64_t dense_element_idx = 0; row < nrows; row++) {
        elements_in_row = 0;
        for (uint32_t col = 0; col < ncols; col++) {
            const double value = M->data[dense_element_idx++];
            if (value != ZERO) {
                out->values[++sparse_element_idx] = value;
                elements_in_row++;
                out->col_ind[sparse_element_idx] = col;
            }
        }
        //just completed one row
        row_ptr += elements_in_row;
        out->row_ptr[row + 1] = row_ptr;
    }
    return out;
}

//only support  symetric matrix - for non symetric, use CRS sparse matrix plse
//QC'ed

SSSMatrix* sss_from_matrix(const Matrix* M) {
    if (!M->symmetric) {
        return NULL; ////TODO message - use CSR
    }
    const double ZERO = 0;
    const uint32_t nrows = M->nrows; // nrows == ncols so no need for ncols
    const uint64_t BUFFER_SIZE = ((nrows * (nrows + 1)) / 2);
    uint32_t nnzs = get_nnz_count(M->data, BUFFER_SIZE);//since matrix is symmetric (compact, lower triangle + diagonal), nnzs is for diagonal + lower (OR upper) triangle
    uint32_t nnzl = 0, nnzd = 0; /* non-zero values in lower half and diagonal respectively  */
    //scan the diagonal to detect number of non zero diagonal elements - this is useful in determining nnzl (the non-zero off diagonal elements
    for (uint32_t row = 0; row < nrows; row++) {
        if (get_matrix_element(M, row, row) != ZERO) {
            nnzd++;
        }
    }
    nnzl = (nnzs - nnzd); // our sym matrix will only store the lower half so nnzs must be nnzl + nnzd
    //const uint32_t nnz_sss = nnzl + nnzl /*+ nnzd  we must store all the diagonal, zero or not so use nrows*/ + nrows;
    //
    //now we know the number of elements, allocate the sparse SSS matrix
    SSSMatrix* out = create_sssmatrix(M->nrows, M->ncols, nnzl);
    if (out == NULL) {
        return NULL;
    }
    //
    for (uint32_t row = 0, row_ptr = 0, sparse_element_idx = 0, elements_in_row = 0; row < nrows; row++) {
        elements_in_row = 0;
        for (uint32_t col = 0; col < row /* LOWER TRIANGLE  */; col++) {
            const double value = get_matrix_element(M, row, col);
            if (value != ZERO) {
                out->values[sparse_element_idx] = value;
                elements_in_row++;
                out->col_ind[sparse_element_idx++] = col;
            }
        }
        //just completed one row (lower half)
        row_ptr += elements_in_row;
        out->row_ptr[row + 1] = row_ptr;
        //
        out->dvalues[row] = get_matrix_element(M, row, row); // diagonal will always be stored
        //dense_element_idx = dense_element_idx + m->ncols - row; // advance to start of next row
    }
    return out;
}

//QC'ed

void print_matrix(const Matrix* M, FILE *stream, const uint32_t rows, const uint32_t cols) {
    if (M == NULL || M->data == NULL) {
        return;
    }
    const uint32_t ROW_LIMIT = min(max(rows, 1), M->nrows);
    const uint32_t COL_LIMIT = min(max(cols, 1), M->ncols);
    const char* NEW_LINE = "\n";
    const char* FORMAT = "%.3g ";
    double value;
    fprintf(stream, NEW_LINE);
    fprintf(stream, NEW_LINE);
    for (uint32_t row = 0; row < ROW_LIMIT; row++) {
        for (uint32_t col = 0; col < COL_LIMIT; col++) {
            value = get_matrix_element(M, row, col);
            fprintf(stream, FORMAT, value);
        }
        fprintf(stream, NEW_LINE);
    }
    fprintf(stream, NEW_LINE);
}

//QC'ed

void print_csr_matrix(const CSRMatrix* M, FILE *stream, const uint32_t n) {
    if (M == NULL || M->values == NULL) {
        return;
    }
    const uint32_t entries = min(max(n, 1), M->nnz);
    const char* FORMAT = "\n(%d, %d) -> %.2g";
    fprintf(stream, "\n***********:: CSR Matrix::\n%d rows, %d columns, %d entries\nPrinting first %d entries:\n", M->nrows, M->ncols, M->nnz, entries);
    for (uint32_t i = 0, row = 0; i < entries;) {
        const uint32_t elements_in_row = M->row_ptr[row + 1] - M->row_ptr[row];
        for (uint32_t ii = 0; ii < elements_in_row && i < entries; ii++) {
            fprintf(stream, FORMAT, row, M->col_ind[i], M->values[i]);
            i++;
        }
        row++;
    }
}

//QC'ed

void print_sss_matrix(const SSSMatrix* M, FILE *stream, const uint32_t n) {
    if (M == NULL || M->values == NULL) {
        return;
    }
    const uint32_t nnzl = M->nnzl;
    const uint32_t entries = min(max(n, 1), nnzl + M->nrows);
    const char* FORMAT = "\n(%d, %d) -> %.2g";
    fprintf(stream, "\n***********:: SSS Matrix::\n%d rows, %d columns, Printing first %d entries:\n", M->nrows, M->ncols,  entries);
    for (uint32_t i = 0, row = 0; i < entries;) {
        const uint32_t elements_in_row = M->row_ptr[row + 1] - M->row_ptr[row];
        for (uint32_t ii = M->row_ptr[row]; ii - M->row_ptr[row] < elements_in_row && i < entries; ii++) {
            fprintf(stream, FORMAT, row, M->col_ind[ii], M->values[ii]);
            i++;
        }
        fprintf(stream, FORMAT, row, row, M->dvalues[row]);
        i++;
        row++;
    }
}

uint32_t get_nnz_count(double* data, const uint64_t BUFFER_SIZE) {
    //http://cs.umw.edu/~finlayson/class/fall16/cpsc425/notes/11-parfor.html
    const double ZERO = 0;
    uint32_t nnz = 0;
    //#pragma omp parallel for reduction(+:nnz)  
    for (uint64_t dense_element_idx = 0; dense_element_idx < BUFFER_SIZE; dense_element_idx++) {
        const double value = data[dense_element_idx];
        if (value != ZERO) {
            nnz++;
        }
    }
    return nnz;
}


