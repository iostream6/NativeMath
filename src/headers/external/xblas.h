/*
 * 2018.04.15  - Based on reorg of legacy C codes written in 2017
 */

/* 
 * File:   xblas.h
 * Author: Ilamah, Osho
 * 
 * Provides driver/interface to system optimised BLAS functions (OPENBLAS and MKL backends), as well as custom BLAS routines (e.g. spMVx) 
 *
 * Created on 15 April 2018, 18:30
 */

#ifndef XBLAS_H
#define XBLAS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "matvec.h"


    //openblas is the default backend, this is overriden by the  __SYSMATH_MKL__ preprocessor directive
#ifdef __SYSMATH_MKL__
    //Using MKL backend
#include "mkl.h"

#define  mkl_init  sysmath_init  
#define  mkl_initialize_random_matrix initialize_random_matrix

#define SEED    777
#define BRNG    VSL_BRNG_MCG31
#define METHOD  VSL_RNG_METHOD_GAUSSIAN_ICDF
#else
    //Using OPENBLAS backend - bring in the openblas headers and OMP for multi-threading support
#include "cblas.h"
#include "lapacke.h"
#include <omp.h>

#include "zrng.h"

#define  openblas_init  sysmath_init  
#define  zrng_initialize_random_matrix initialize_random_matrix

#endif  //end __SYSMATH_MKL__ ifdef_else_endif block

    /**
     * Specifies the number of threads to be used by BLAS backend 
     * @param num_threads
     */
    void set_system_num_threads(const uint32_t num_threads);

    /**
     * Initializes diagonal matrix from vector data 
     * @param D
     * @param v
     */
    void initialize_diagonal_matrix(Matrix* D, Vector* v);

    /**
     * Computes y = M*x ; row major 
     * @param M
     * @param x
     * @param y
     */
    void matrix_vector_mult(const Matrix* M, const double* x, double* y);

    /**
     * Performs  a sparse matrix  to dense vector mulplication operation: y = M*x ; row major 
     * @param M Input sparse matrix
     * @param x Input dense vector
     * @param y Output vector
     */
    void sparse_matrix_vector_mult(const CSRMatrix* M, double* x, double* y);

    /**
     * Performs  a symmetric sparse matrix to dense vector mulplication operation: y = M*x ; It is assumed that the output vector has been allocated by the client 
     * @param M Input symmetric sparse matrix in SSS format
     * @param x Input dense vector
     * @param y Output vector
     */
    void sparse_sym_matrix_vector_mult(const SSSMatrix* M, const double* x, double* y, const uint32_t threads);


    /**
     * Computes C = A*B ; assuming ROW major format 
     * @param A
     * @param B
     * @param C
     */
    void matrix_matrix_mult(Matrix* A, Matrix* B, Matrix* C);

    /**
     * Computes C = A^T*B ; assuming ROW major format  
     * @param A
     * @param B
     * @param C
     */
    void matrix_transpose_matrix_mult(Matrix* A, Matrix* B, Matrix* C);

    /**
     * Computes C = A*B^T ; assuming ROW major format  
     * @param A
     * @param B
     * @param C
     */
    void matrix_matrix_transpose_mult(Matrix* A, Matrix* B, Matrix* C);

    /**
     * Computes Q from [Q,R] = qr(M,'0') compact QR factorization, assuming ROW major format  
     * @param M mxn matrix
     * @param Q mxn matrix
     */
    void QR_factorization_getQ(Matrix* M, Matrix* Q);

    /**
     * Computes [Q,R] = qr(M,'0'), or  the compact QR factorization 
     * @param M mxn matrix
     * @param Q mxn matrix
     * @param R min(m,n) x min(m,n)
     */
    void compact_QR_factorization(Matrix* M, Matrix* Q, Matrix* R);

    /**
     * Copies elements from the source to the destination matrix
     * @param destination
     * @param source
     */
    void matrix_copy(Matrix* destination, Matrix* source);

    /**
     * Extracts a column of the provided matrix into a vector
     * @param M the matrix
     * @param j the column to extract
     * @param vector the vector to polulate
     */
    void get_matrix_col(Matrix* M, const uint32_t j, Vector* vector);


    /**
     * Computes the SVD of A = U*Σ*V^T for real routines, assuming row major format
     * @param M
     * @param U
     * @param svals vector containing the diagonal elements of Σ, in descending order of magnitude.
     * @param Vt
     */
    void singular_value_decomposition(Matrix* M, Matrix *U, Vector* svals /*Matrix *S*/, Matrix *Vt);

    //int compute_evals_and_evecs_of_symm_matrix(Matrix *S, vec *evals);

    /**
     * Constitutes a matrix P, using the SVD factors U, S and V. This is useful for validating SVD accuracy.
     * @param U
     * @param S
     * @param V
     * @param P
     */
    void form_svd_product_matrix(Matrix* U, Matrix* S, Matrix* V, Matrix* P);

    /**
     * Computes the frobenius norm of a matrix. This is useful for validating decomposition accuracy.
     * @param M
     * @return the frobenius norm
     */
    double get_matrix_frobenius_norm(Matrix* M);

    void sysmath_init();

    void initialize_random_matrix(Matrix* M);
#ifdef __cplusplus
}
#endif

#endif /* XBLAS_H */

