/* 
 * File:   xmmio.c
 * Author: Ilamah, Osho
 *
 * Created on 17 April 2018, 09:18
 */

#include <stdio.h>
//#include <omp.h>
//


#include "../headers/external/xmmio.h"
#include "../headers/external/mmio.h"
#include "../headers/external/xalloc.h"

COOMatrix* create_coomatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz);
void coomatrix_delete(COOMatrix* matrix);
//
COOMatrix* load_coo_matrix(char* filename, MM_typecode* typecode);
CSRMatrix* coo_to_csr(const COOMatrix* coo_matrix);

void print_coo_matrix(const COOMatrix* m, const uint32_t n);
void sort_row_by_column(uint32_t* col_idx, double* a, int start, int end);

/**
 * Reads data from matrix market formatted text file into Matrix data structure. Both dense and sparse matrix market files can be read.
 * 
 * @param filename the absolute path to the input file
 * @return a Matrix structure containing the input data if successful, otherwise NULL
 */
//QC'ed

Matrix* load_mm_matrix(const char* filename) {
    FILE* file;
    MM_typecode mmcode;
    int rows, cols, nnz;
    uint32_t read_entries, i, j;
    double dvalue;
    // Open the matrix file
    if ((file = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "MMIO DRIVER: can't open file \"%s\"\n", filename);
        return NULL;
    }
    //
    if (mm_read_banner(file, &mmcode) != 0) {
        fprintf(stderr, "MMIO DRIVER: Could not process Matrix Market banner.\n");
        return NULL;
    }
    //
    //Check if the matrix type is one we support
    if (!mm_is_matrix(mmcode)
            || mm_is_complex(mmcode)/* complex number matrix */
            //|| !mm_is_coordinate(mmcode) /* coordinate - i.e. ascii format */
            || !(mm_is_real(mmcode) || mm_is_integer(mmcode)) /* floating point or integral numbers */
            ) {
        fprintf(stderr, "MMIO DRIVER: Unsupported Matrix Market format: [%s]\n", mm_typecode_to_str(mmcode));
        return NULL;
    }
    if (mm_is_sparse(mmcode)) { //coordinate formatted MM file
        // Read Matrix dimensions
        if ((mm_read_mtx_crd_size(file, &rows, &cols, &nnz)) != 0) {
            fprintf(stderr, "MMIO DRIVER: Could not read matrix size.\n");
            return NULL;
        }
        Matrix* matrix = create_matrix(rows, cols, mm_is_symmetric(mmcode));
        if (matrix == NULL) {
            fprintf(stderr, "MMIO DRIVER: Could not create %dx%d matrix.\n", rows, cols);
            return NULL;
        }
        const int ELEMENTS_PER_RECORD = 3;
        //    const char* error_message = "MMIO DRIVER: Could not read data.\n";
        // Read the matrix data
        for (read_entries = 0, i = 0, j = 0; read_entries < nnz; read_entries++) {
            if (fscanf(file, "%d %d %lg\n", &i, &j, &dvalue) < ELEMENTS_PER_RECORD) {
                break;
            }
            set_matrix_element(matrix, --i, --j, dvalue); //decrement to bring back to 1 based indices
        }
        fclose(file);
        return matrix;
    } else {//array formatted MM file, we expect
        // Read Matrix dimensions
        if ((mm_read_mtx_array_size(file, &rows, &cols)) != 0) {
            fprintf(stderr, "MMIO DRIVER: Could not read matrix size.\n");
            return NULL;
        }
        Matrix* matrix = create_matrix(rows, cols, mm_is_symmetric(mmcode));
        if (matrix == NULL) {
            fprintf(stderr, "MMIO DRIVER: Could not create %dx%d matrix.\n", rows, cols);
            return NULL;
        }
        const int ELEMENTS_PER_RECORD = 1;
        // Read the matrix data
        uint32_t expected = 0;
        if (mm_is_symmetric(mmcode)) {
            //the below assumes there are no comment lines after the control info line
            //integer or double "array" (i.e. not co-ord) format:: Column major of lower triangle
            for (read_entries = 0, j = 0; j < cols; j++) {//column order is how they are stored
                for (i = j /*lower triangle only */; i < rows; i++) {
                    if (fscanf(file, "%lg\n", &dvalue) < ELEMENTS_PER_RECORD) {
                        break;
                    }
                    set_matrix_element(matrix, i, j, dvalue);
                    read_entries++;
                }
            }
            expected = (rows * (rows + 1)) / 2;
        } else {
            //integer or double "array" (i.e. not co-ord) format:: Column major of entire matrix
            for (read_entries = 0, j = 0; j < cols; j++) {//column order is how they are stored
                for (i = 0; i < rows; i++) {
                    if (fscanf(file, "%lg\n", &dvalue) < ELEMENTS_PER_RECORD) {
                        break;
                    }
                    set_matrix_element(matrix, i, j, dvalue);
                    read_entries++;
                }
            }
            expected = rows * cols;
        }
        fclose(file);
        if (read_entries != expected) {
            //array formatted  (i.e. dense) non-symmetric matrix, should have nrows*ncols entries
            fprintf(stderr, "MMIO DRIVER: Inssuficient entries (%d) encountered, %d expected.\n", read_entries, expected);
            delete_matrix(matrix);
            return NULL;
        }
        return matrix;
    }
}

/**
 * Reads data from matrix market formatted text file into CSRMatrix data structure. Only sparse, non-symmetric matrix market input data is supported. If sparse symmetric input is required, please use the Matrix and SSSMatrix structure and associated 
 * read methods instead.
 * 
 * @param filename the absolute path to the input file
 * @return  a CSRMatrix structure containing the input data if successful, otherwise NULL
 */
//QC'ed

CSRMatrix* load_mm_csrmatrix(char* filename) {
    MM_typecode typecode;
    COOMatrix* coo_matrix = load_coo_matrix(filename, &typecode);
    if (coo_matrix == NULL) {
        return NULL;
    }
    if (mm_is_symmetric(typecode)) {
        //Symmetric COO to CSR 
        coomatrix_delete(coo_matrix);
        fprintf(stderr, "MMIO DRIVER: SSS Matrix should be used for symmetric cases");
        return NULL;
    }
    if (!mm_is_sparse(typecode)) {
        coomatrix_delete(coo_matrix);
        fprintf(stderr, "MMIO DRIVER: CSR Matrix is optimized for sparse matrices");
        return NULL;
    }
    CSRMatrix* csr_matrix = coo_to_csr(coo_matrix);
    coomatrix_delete(coo_matrix);
    return csr_matrix;
}

/**
 * Reads data from matrix market formatted text file into SSSMatrix data structure. Only sparse, symmetric matrix market input data is supported. If sparse NON-symmetric input is required, please use the Matrix and CSRMatrix structure and associated 
 * read methods instead.
 * 
 * @param filename the absolute path to the input file
 * @return  a SSSMatrix structure containing the input data if successful, otherwise NULL
 */
SSSMatrix* load_mm_ssmatrix(char* filename) {
    MM_typecode typecode;
    COOMatrix* coo_matrix = load_coo_matrix(filename, &typecode);
    if (coo_matrix == NULL) {
        return NULL;
    }
    if (!mm_is_symmetric(typecode)) {
        fprintf(stderr, "MMIO DRIVER: SSS Matrix must be symmetric");
        return NULL;
    }
    if (!mm_is_sparse(typecode)) {
        coomatrix_delete(coo_matrix);
        fprintf(stderr, "MMIO DRIVER: SSS Matrix is optimized for sparse matrices");
        return NULL;
    }
    //coo_matrix.nnzl will be greater than 0
    const uint32_t nrows = coo_matrix->nrows, ncols = coo_matrix->ncols;
    const uint32_t nnzr = coo_matrix->nnzr;

    if (coo_matrix->nnzr != coo_matrix->nnzl + coo_matrix->nnzd) {
        fprintf(stderr, "MMIO DRIVER: SSS Matrix minimum data elements not satisfied");
        coomatrix_delete(coo_matrix);
        return NULL;
    }
    SSSMatrix* sss_matrix = create_sssmatrix(nrows, ncols, coo_matrix->nnzl);
    if (sss_matrix == NULL) {
        fprintf(stderr, "MMIO DRIVER: SSS Matrix could not be created");
        coomatrix_delete(coo_matrix);
        return NULL;
    }
    /*================================================================================================
     * The implementation below is based off of the coo_to_csr function, adapted for SSS matrix format
     * ===============================================================================================
     */

    // compute number of non-zero, off diagonal lower triangular entries per row of A, this is achieved by scanning all elements of 
    // the COO and incrementing row_ptr for the row containing each element
    // after this stage, each entry in row_ptr  is the number of elements in the corresponding LT row of the COO matrix
    for (uint32_t ii = 0; ii < nnzr; ii++) {
        const uint32_t row = coo_matrix->row_coord[ii];
        const uint32_t col = coo_matrix->col_coord[ii];
        if (row < col) {
            continue;
        }
        if (row > col) {
            //lower triangle element - store
            sss_matrix->row_ptr[row]++;
        } else { //(row == col)
            //diagonal element, save once and for all here
            sss_matrix->dvalues[row] = coo_matrix->values[ii];
        }
    }
    //
    // We now convert the row_ptr to the actual SSS row_ptr definition - running sum which defines the number of stored elements in the previous row
    // this is achieved by computing the running sum, starting at 0
    for (uint32_t row = 0, running_sum = 0; row < nrows; row++) {
        const uint32_t temp = sss_matrix->row_ptr[row]; // temp is the number of stored (i.e. LT) elements in the ith row
        sss_matrix->row_ptr[row] = running_sum; //SSS row_ptr definition
        running_sum += temp; // update running so that next row will take into account the number of stored elements in the ith row
    }
    sss_matrix->row_ptr[nrows] = coo_matrix->nnzl; // fill in the nrows+1 row_ptr value, which must be nnzl (not nnz!!)
    // At this stage, sss_matrix->row_ptr corresponds to the SSS definition, we now generate the col_ind and values arrays from the COO
    for (uint32_t n = 0; n < nnzr; n++) {
        const uint32_t row = coo_matrix->row_coord[n]; // the row address
        const uint32_t col = coo_matrix->col_coord[n];
        if (row <= col) {//Ignore row <= col as eq case has already been handled and < case is UT entry which is not required in SSS
            continue;
        }
        //lower triangle element - store
        // We determine the col_ind address by looking at row_ptr  for the target row, it tells us where the target row should start in col_ind
        const uint32_t slot = sss_matrix->row_ptr[row];
        // we then save the col_coord and values in this address slot 
        sss_matrix->col_ind[slot] = col;
        sss_matrix->values[slot] = coo_matrix->values[n];
        // now, we play a little trick - we shift the row_ptr (i.e. record of start address of this row, so that the next element in this row will be saved in the next slot
        // this screws up row_ptr but we can easily correct this (see for loop below), and it facilitates slot determination for col_ind and values
        sss_matrix->row_ptr[row]++;
    }
    //correct row_ptr
    for (uint32_t i = 0, last = 0; i <= nrows; i++) {
        const uint32_t temp = sss_matrix->row_ptr[i];
        sss_matrix->row_ptr[i] = last;
        last = temp;
    }
    // Sort each row group by column indices to ensure the natural column order in row major format is preserved. 
    for (uint32_t i = 0; i < nrows; i++) {
        sort_row_by_column(sss_matrix->col_ind, sss_matrix->values, (int) sss_matrix->row_ptr[i], (int) sss_matrix->row_ptr[i + 1]);
    }
    coomatrix_delete(coo_matrix);
    return sss_matrix;
}

/* ****************************************************************************
 *                     INTERBAL (PRIVATE) ROUTINES
 *                     INTERBAL (PRIVATE) ROUTINES
 *                     INTERBAL (PRIVATE) ROUTINES 
 * ***************************************************************************
 */


/**
 * Sorts the subset of packed column index buffer, corresponding to elements of a sparse matrix row, and the associated element values in the packed values buffer, in natural order
 * @param col_idx the column index buffer
 * @param a the packed values buffers
 * @param start offset into the column buffer at which the row entries begin
 * @param end offset into the column buffer at which the row entries end
 */
//QC'ed

void sort_row_by_column(uint32_t* col_idx, double* values, int start, int end) {
    int i, j, it;
    double dt;
    for (i = end - 1; i > start; i--) {
        for (j = start; j < i; j++) {
            if (col_idx[j] > col_idx[j + 1]) {
                if (values) {
                    dt = values[j];
                    values[j] = values[j + 1];
                    values[j + 1] = dt;
                }
                it = col_idx[j];
                col_idx[j] = col_idx[j + 1];
                col_idx[j + 1] = it;
            }
        }
    }
}


/**
 * Allocates memory for a co-ordinate format matrix data structure
 * @param nrows the total number of rows in the matrix
 * @param ncols the total number of cols in the matrix
 * @param nnz the total number of non-zero entries in the matrix
 * @return A pointer to a COOMatrix data structure with sufficient backing memory to hold the specified matrix
 */
//QC'ed

COOMatrix* create_coomatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz) {
    COOMatrix* m = malloc(sizeof (COOMatrix));
    m->values = (double*) aligned_calloc(nnz, sizeof (double), 64);
    m->col_coord = (uint32_t*) aligned_calloc(nnz, sizeof (uint32_t), 64);
    m->row_coord = (uint32_t*) aligned_calloc(nnz, sizeof (uint32_t), 64);
    if (m->values == NULL || m->col_coord == NULL || m->row_coord == NULL) {
        return NULL;
    }
    m->nrows = nrows;
    m->ncols = ncols;
    m->nnz = nnz;
    m->nnzd = 0;
    m->nnzl = 0;
    m->nnzr = 0;
    return m;
}


/**
 * Frees memory for a co-ordinate format matrix data structure
 * @param matrix pointer to the COOMatrix to be freed
 */
//QC'ed

void coomatrix_delete(COOMatrix* matrix) {
    if (matrix != 0) {
        aligned_free(matrix->values);
        aligned_free(matrix->row_coord);
        aligned_free(matrix->col_coord);
        free(matrix);
        matrix = 0; //avoids double free's
    }
}


/**
 * Reads data from matrix market formatted text file into a COOMatrix data structure. This function provides an efficient path to load sparse matrices without allocating the full ncols*nrows buffer, instead using COOMatrix of much fewer i, j, value triplet. It 
 * is intended for use with sparse matrices only, for dense matrix, please use the Matrix structure and associated read method instead.
 * 
 * @param filename the absolute path to the input file
 * @param typecode pointer to matrix market typecode variable allowing typecode info to be pushed back to calling client
 * @return NULL a COOMatrix structure containing the input data if successful, otherwise NULL
 */
//QC'ed

COOMatrix* load_coo_matrix(char* filename, MM_typecode* typecode) {
    FILE* file;
    int rows, cols, nnz;
    // Open the matrix file
    if ((file = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "MMIO DRIVER: can't open file \"%s\"\n", filename);
        return NULL;
    }
    //
    if (mm_read_banner(file, typecode) != 0) {
        fprintf(stderr, "MMIO DRIVER: Could not process Matrix Market banner.\n");
        return NULL;
    }
    //
    //Check if the matrix type is one we support
    if (!mm_is_matrix(*typecode)
            || mm_is_complex(*typecode)/* complex number matrix */
            || !mm_is_coordinate(*typecode) /* coordinate - i.e. ascii format */
            || !(mm_is_real(*typecode) || mm_is_integer(*typecode)) /* floating point or integral numbers */
            ) {
        fprintf(stderr, "MMIO DRIVER: Unsupported Matrix Market format: [%s]\n", mm_typecode_to_str(*typecode));
        return NULL;
    }

    bool symmetric = mm_is_symmetric(*typecode);

    // Read Matrix dimensions
    if ((mm_read_mtx_crd_size(file, &rows, &cols, &nnz)) != 0) {
        fprintf(stderr, "MMIO DRIVER: Could not read matrix size.\n");
        return NULL;
    }

    COOMatrix* matrix = create_coomatrix((uint32_t) rows, (uint32_t) cols, (uint32_t) nnz);

    if (matrix == NULL) {
        fprintf(stderr, "MMIO DRIVER: Could not create %dx%d COO Matrix.\n", rows, cols);
        return NULL;
    }

    //const char* error_message = "MMIO DRIVER: Could not read data.\n";
    const int ELEMENTS_PER_RECORD = 3;
    // Read the matrix data
    uint32_t nnzd = 0, nnzl = 0, nnzr = 0;
    for (uint32_t i = 0, j = 0; nnzr < nnz; nnzr++) {
        if (fscanf(file, "%d %d %lg\n", &i, &j, &matrix->values[nnzr]) < ELEMENTS_PER_RECORD) {
            break;
        }
        matrix->row_coord[nnzr] = --i; //decrement to bring back to 1 based indices
        matrix->col_coord[nnzr] = --j;
        if (symmetric) {
            if (i == j) {
                nnzd++;
            } else if (i > j) {
                nnzl++;
            }
        }
    }
    matrix->nnzr = nnzr;
    if (symmetric) {
        matrix->nnzd = nnzd;
        matrix->nnzl = nnzl;
        matrix->nnz = nnzl + nnzl + nnzd;

        fprintf(stdout, "MMIO DRIVER: SSS Matrix with %d zero valued elements along diagonal\n", rows - nnzd);

    } else {
        matrix->nnz = nnzr;
        //negative nnzl implies non-symmetric matrix
        matrix->nnzd = 0;
        matrix->nnzl = 0;
    }
    fclose(file);
    return matrix;
}

/**
 * Generates a CSR Matrix corresponding to the input COO matrix. It is assumed that there are no duplicates in the coo matrix. Row and Column indices of the COO Matrix *are not* assumed to be ordered.
 * This implementation is inspired by https://github.com/scipy/scipy/blob/3b36a57/scipy/sparse/sparsetools/coo.h#L34. This method should <b>NOT</b> be called if the input COO matrix only contains lower or upper triangular data
 * as the resulting CSR matrix will be incomplete.
 * 
 * @param coo_matrix the input COOMatrix
 * @return a CSR Matrix corresponding to the input COO matrix
 */
//QC'ed

CSRMatrix* coo_to_csr(const COOMatrix* coo_matrix) {
    if (coo_matrix->nnzl > 0) {
        fprintf(stderr, "MMIO DRIVER: COO to CSR expects non-symmetric input.\n");
    }
    const uint32_t nnz = coo_matrix->nnz, nrows = coo_matrix->nrows, ncols = coo_matrix->ncols;
    CSRMatrix* csr_matrix = create_csrmatrix(nrows, ncols, nnz);
    if (csr_matrix == NULL) {
        return csr_matrix;
    }
    // compute number of non-zero entries per row of A, this is achieved by scanning all elements of the COO and incrementing row_ptr for the row containing each element
    // after this stage, each entry in row_ptr  is the number of elements in the corresponding row of the COO matrix
    for (uint32_t n = 0; n < nnz; n++) {
        csr_matrix->row_ptr[coo_matrix->row_coord[n]]++;
    }
    //
    // We now convert the row_ptr to the actual CSR row_ptr definition - running sum which defines the number of elements in the previous row
    // this is achieved by computing the running sum, starting at 0
    for (uint32_t row = 0, running_sum = 0; row < nrows; row++) {
        const uint32_t temp = csr_matrix->row_ptr[row]; // temp is the number of elements in the ith row
        csr_matrix->row_ptr[row] = running_sum; //CSR row_ptr definition
        running_sum += temp; // update running so that next row will take into account the number of elements in the ith row
    }
    csr_matrix->row_ptr[nrows] = nnz; // fill in the nrows+1 row_ptr value, which must be nnz
    // At this stage, csr_matrix->row_ptr corresponds to the CSR definition, we now generate the col_ind and values arrays from the COO
    //
    for (uint32_t n = 0; n < nnz; n++) {
        const uint32_t row = coo_matrix->row_coord[n]; // the row address
        // We determine the col_ind address by looking at row_ptr  for the target row, it tells us where the target row should start in col_ind
        const uint32_t slot = csr_matrix->row_ptr[row];
        // we then save the col_coord and values in this address slot 
        csr_matrix->col_ind[slot] = coo_matrix->col_coord[n];
        csr_matrix->values[slot] = coo_matrix->values[n];
        // now, we play a little trick - we shift the row_ptr (i.e. record of start address of this row, so that the next element in this row will be saved in the next slot
        // this screws up row_ptr but we can easily correct this (see for loop below), and it facilitates slot determination for col_ind and values
        csr_matrix->row_ptr[row]++;
    }
    //correct row_ptr
    for (uint32_t i = 0, last = 0; i <= nrows; i++) {
        const uint32_t temp = csr_matrix->row_ptr[i];
        csr_matrix->row_ptr[i] = last;
        last = temp;
    }
    // Sort each row group by column indices to ensure the natural column order in row major format is preserved. 
    for (uint32_t i = 0; i < nrows; i++) {
        sort_row_by_column(csr_matrix->col_ind, csr_matrix->values, (int) csr_matrix->row_ptr[i], (int) csr_matrix->row_ptr[i + 1]);
    }
    //
    return csr_matrix;
}

///**
// * Prints 1 or more entries from a COOMatrix
// * @param m  pointer to the COOMatrix
// * @param n the number of entries to print
// */
//void print_coo_matrix(const COOMatrix* m, const uint32_t n) {
//    const uint32_t entries = min(max(n, 1), m->nnz);
//    const char* FORMAT = "\n(%d, %d) -> %.2g";
//    printf("\n***********:: COO Matrix::\n%u rows, %u columns\nPrinting first %u entries:\n", m->nrows, m->ncols, entries);
//    for (uint32_t i = 0; i < entries; i++) {
//        printf(FORMAT, m->row_coord[i], m->col_coord[i], m->values[i]);
//    }
//    printf("\n***********");
//}
