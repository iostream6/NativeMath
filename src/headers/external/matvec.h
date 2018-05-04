/*
 * 2017.04.23  - Created initial version
 * 2017.05.22  - Reimplemented for portable BLAS/LAPACK backend
 * 2017.06.24  - Added basic sparse matrix support
 * 2018.04.15  - Improved Sparse matrix support
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

    Vector* create_vector(uint32_t nrows);

    Matrix* create_matrix(uint32_t nrows, uint32_t ncols, bool symmetric);

    CSRMatrix* create_csrmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz);

    SSSMatrix* create_sssmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz);

    //
    Matrix* wrap_as_matrix(uint32_t nrows, uint32_t ncols, bool symmetric, double* data);

    CSRMatrix* wrap_as_csrmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz, double* data, uint32_t* row_ptr, uint32_t* col_ind);

    SSSMatrix* wrap_as_sssmatrix(uint32_t nrows, uint32_t ncols, uint32_t nnz, double* dvalues, double* values, uint32_t* row_ptr, uint32_t* col_ind);

    //
    void delete_vector(Vector* vector);

    void delete_matrix(Matrix* matrix);

    void delete_csrmatrix(CSRMatrix* matrix);

    void delete_sssmatrix(SSSMatrix* matrix);


    /*
     *
     *
     *
     */

    /**
     * Sets the value of the provided matrix, at the specified row and column numbers. Assumes row major format
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

    Matrix* unpack_symmetric_matrix(Matrix* M);
    
    CSRMatrix* csr_from_matrix(const Matrix* M);

    SSSMatrix* sss_from_matrix(const Matrix* M);

    void print_matrix(const Matrix* M, FILE* stream, const uint32_t rows, const uint32_t cols);

    void print_csr_matrix(const CSRMatrix* M, FILE* stream, const uint32_t n);

    void print_sss_matrix(const SSSMatrix* M, FILE* stream, const uint32_t n);



#ifdef __cplusplus
}
#endif

#endif /* MATVEC_H */

