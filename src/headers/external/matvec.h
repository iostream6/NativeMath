/*
 * 2017.04.23  - Created initial version
 * 2017.05.22  - Reimplemented for portable BLAS/LAPACK backend
 * 2017.06.24  - Added basic sparse matrix support (CSR format)
 * 2018.04.15  - Improved Sparse matrix support (introduced SSS format)
 */

/* 
 * File:   matvec.h
 * Author: Ilamah, Osho
 *
 * Created on 23 April 2017, 08:42 as matvec.h
 */
//https://stackoverflow.com/questions/3187957/how-to-store-a-symmetric-matrix
//https://stackoverflow.com/questions/9039189/make-efficient-the-copy-of-symmetric-matrix-in-c-sharp/9040526#9040526
#ifndef MATVEC_H
#define MATVEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))

    typedef struct {
        uint32_t nrows;
        double* data;
    } Vector;

    typedef struct {
        uint32_t nrows, ncols;
        bool symmetric;
        double* data;
    } Matrix;

    /*
     * A CSR (or CSR3/Yale) sparse matrix representation
     */
    typedef struct {
        uint32_t nrows, ncols; //number of cols and rows
        uint32_t nnz;
        uint32_t* row_ptr;
        uint32_t* col_ind;
        double* values;

    } CSRMatrix;

    /*
     * A SS sparse matrix representation
     */
    typedef struct {
        uint32_t nrows, ncols; //number of cols and rows
        uint32_t nnzl;
        uint32_t* row_ptr;
        uint32_t* col_ind;
        double* values;
        double* dvalues;
    } SSSMatrix;

    /**
     * Allocates memory for a Vector data structure and initializes all the elements to zero
     * 
     * @param nrows the number of rows (i.e. elements) in the vector
     * @return  Pointer to a Vector data structure with sufficient backing memory to hold the specified number of elements
     */
    Vector* create_vector(uint32_t nrows);


    /**
     * Creates a new matrix and initializes all the elements to zero. If symmetric flag is set, only the upper half of the matrix is stored
     * 
     * @param nrows the number of rows contained in the matrix
     * @param ncols the number of columns contained in the matrix
     * @param symmetric if true, the matrix is taken as symmetric and internally only the upper triangle  (+ diagonal) elements are stored, ROW MAJOR
     * @return Pointer to a Matrix data structure with sufficient backing memory 
     */
    Matrix* create_matrix(uint32_t nrows, uint32_t ncols, bool symmetric);

    /**
     * Creates a new CSRMatrix and initializes all the elements to zero. 
     * 
     * @param nrows the number of rows contained in the matrix
     * @param ncols the number of columns contained in the matrix
     * @param nnz the total number of non-zero elements in the represented matrix
     * @return  Pointer to a CSRMatrix data structure with sufficient backing memory
     */
    CSRMatrix* create_csrmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz);

    /**
     * Creates a new SSSMatrix and initializes all the elements to zero. 
     * @param nrows the number of rows contained in the matrix
     * @param ncols the number of columns contained in the matrix
     * @param nnzl the total number of non-zero elements in the lower triangle (excluding the diagonal) of the represented matrix
     * @return Pointer to a SSSMatrix data structure with sufficient backing memory
     */
    SSSMatrix* create_sssmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnzl);

    /**
     * Wraps a double buffer into a Matrix structure such that the buffer may be read/written/used as a Matrix
     * 
     * @param nrows the number of rows in the wrapped matrix
     * @param ncols the number of columns in the wrapped matrix
     * @param symmetric if true, the matrix is taken as symmetric and internally only the upper triangle  (+ diagonal) elements are assumed to be stored, ROW MAJOR 
     * @param data the buffer to be wrapped
     * @return Pointer to a Matrix data structure with the data buffer serving as backing memory 
     */
    Matrix* wrap_as_matrix(uint32_t nrows, uint32_t ncols, bool symmetric, double* data);

    /**
     * Wraps a double buffer and associated column index and row pointer buffers into a CSRMatrix structure such that the buffer may be read/written/used as a CSRMatrix
     * 
     * @param nrows the number of rows in the wrapped matrix
     * @param ncols the number of columns in the wrapped matrix
     * @param nnz the total number of non-zero elements in the represented matrix
     * @param data the data buffer to be wrapped
     * @param row_ptr the CSR row_pointer associated with the wrapped data buffer
     * @param col_ind the CSR col_ind associated with the wrapped data buffer
     * @return Pointer to a CSRMatrix data structure with the data buffer serving as backing memory 
     */
    CSRMatrix* wrap_as_csrmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz, double* data, uint32_t* row_ptr, uint32_t* col_ind);

    /**
     * Wraps relevant buffers into a SSSMatrix structure such that the buffers may be read/written/used as a SSSMatrix
     * 
     * @param nrows the number of rows in the wrapped matrix
     * @param ncols the number of columns in the wrapped matrix
     * @param nnzl the total number of non-zero elements in the lower triangle (excluding the diagonal) of the represented matrix
     * @param dvalues buffer of diagonal values
     * @param values buffer of non-diagonal values
     * @param row_ptr row pointer buffer
     * @param col_ind col ind buffer
     * @return Pointer to a SSSMatrix data structure with the data buffers serving as backing memory 
     */
    SSSMatrix* wrap_as_sssmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnzl, double* dvalues, double* values, uint32_t* row_ptr, uint32_t* col_ind);

    /**
     * Frees memory for a Vector data structure
     * 
     * @param vector pointer to the structure to be freed
     */
    void delete_vector(Vector* vector);

    /**
     * Frees memory for a Matrix data structure
     * 
     * @param matrix pointer to the structure to be freed
     */
    void delete_matrix(Matrix* matrix);

    /**
     * Frees memory for a CSRMatrix data structure
     * 
     * @param matrix pointer to the structure to be freed
     */
    void delete_csrmatrix(CSRMatrix* matrix);

    /**
     * Frees memory for a SSSMatrix data structure
     * 
     * @param matrix pointer to the structure to be freed
     */
    void delete_sssmatrix(SSSMatrix* matrix);



    /**
     * Sets the value of the provided matrix, at the specified row and column indices. It is assumed that the backing array is organized in ROW MAJOR format from left to right and then top to down.  If the 
     * matrix is symmetric, only a portion is stored - the upper triangular half (+ diagonal) in row major format, which corresponds to the lower triangular half (+ diagonal) in COLUMN major format.
     * 
     * @param M
     * @param row_num
     * @param col_num
     * @param val
     */
    void set_matrix_element(Matrix* M, const uint32_t row_num, const uint32_t col_num, const double val);


    /**
     * Gets the value of the provided matrix, at the specified row and column numbers. Assumes row major format
     * @param M
     * @param row_num
     * @param col_num
     * @return 
     */
    double get_matrix_element(const Matrix* M, const uint32_t row_num, const uint32_t col_num);

    /**
     * Unpacks a truncatedly stored (i.e. only half triangle + diag) symmetric Matrix into an untruncatedly stored one.
     * 
     * @param M the input matrix to unpack
     * @return A newly allocated unpacked matrix if the input matrix was symmetric, else the input matrix is returned
     */
    Matrix* unpack_symmetric_matrix(Matrix* M);

    /**
     * Creates a CSRMatrix from the specified Matrix data structure. This function only supports non-symmetric input matrix, for symmetric input matrix, please use the SSSMatrix structure and associated method instead.
     * 
     * @param M the input Matrix
     * @return  a newly allocated CSRMatrix equivalent representation of the input matrix  
     */
    CSRMatrix* csr_from_matrix(const Matrix* M);

    /**
     * Creates a SSSMatrix from the specified Matrix data structure. This function only supports symmetric input matrix, for non-symmetric input matrix, please use the CSRMatrix structure and associated method instead.
     * 
     * @param M the input Matrix
     * @return a newly allocated SSSMatrix equivalent representation of the input matrix  
     */
    SSSMatrix* sss_from_matrix(const Matrix* M);


    /**
     * Prints a portion of a Matrix.
     * 
     * @param M the input matrix
     * @param stream the output stream to which the matrix is printed
     * @param rows the number of rows from the top of the matrix to be included in the printed portion
     * @param cols the number of cols from the left of the matrix to be included in the printed portion
     */
    void print_matrix(const Matrix* M, FILE* stream, const uint32_t rows, const uint32_t cols);

    /**
     * Prints a portion of a CSRMatrix.
     * 
     * @param M the input matrix
     * @param stream the output stream to which the matrix is printed
     * @param n number of entries in the sparse matrix to be included in the printed portion
     */
    void print_csr_matrix(const CSRMatrix* M, FILE* stream, const uint32_t n);

    /**
     * Prints a portion of a SSSMatrix.
     * 
     * @param M the input matrix
     * @param stream the output stream to which the matrix is printed
     * @param n number of entries in the sparse matrix to be included in the printed portion
     */
    void print_sss_matrix(const SSSMatrix* M, FILE* stream, const uint32_t n);



#ifdef __cplusplus
}
#endif

#endif /* MATVEC_H */

