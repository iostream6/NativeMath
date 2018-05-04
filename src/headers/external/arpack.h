/* 
 * File:   arpack_driver.h
 * Author: Ilamah, Osho
 *
 * Created on 06 May 2017, 21:37
 * 2017.06.24  - Added basic sparse matrix support for Eigen routines
 * 2018.04.27  - Added SSS sparse symmetric support
 */
#ifndef ARPACK_DRIVER_H
#define ARPACK_DRIVER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "xblas.h"

    //If we are on linux then we know both arpack64 and app are compiled with GCC. otherwise if we are on windows and we are using openblas backend, 
    //we know both arpack64(openblas flavour) and app are compiled with GCC. 
    // http://stackoverflow.com/a/4605893
    // https://stackoverflow.com/a/2998876
    //#ifdef  __linux__
#if defined(__linux__) || defined(__SYSMATH_OPENBLAS__)   
    //mangles arpack subroutine calls to lower case and underscore for gcc/gfortran compatibility     
#ifndef FC_MANGLE_GCC
#define FC_MANGLE_GCC

#define  DSAUPD  dsaupd_
#define  DSEUPD  dseupd_
#define  DNAUPD  dnaupd_
#define  DNEUPD  dneupd_

#endif /* FC_MANGLE_GCC */

#endif  /* __linux__ || __SYSMATH_OPENBLAS__*/

    /**
     * @brief Ordering to return eigenvalues in.
     */
    typedef enum EigenOrders {
        LA, /**< Largest algebraic eigenvalues first. */
        SA, /**< Smallest algebraic eigenvalues first. */
        LM, /**< Largest eigenvalues in magnitude first. */
        SM, /**< Smallest eigenvalues in magnitude first. */
        BE, /**< Compute nev eigenvalues, half from each end of the spectrum. When nev is odd, compute one more from the high end than from the low end. */
    } EigenOrder;

    /**
     * Mode of operation of the ARPACK backend.
     *
     * @note Note that currently G modes (generalised problems) are not yet implemented.
     */
    typedef enum EigenSolverModes {
        /* For use in solving Av = vd */
        I_REGULAR, /**< For solving Av=vd in regular mode. */
        //   I_SHIFTINVERT, /**< For solving Av=vd in shift-invert mode. */
        /* For use in solving Av = dMv */
        //   G_REGINVERSE,  /**< For solving Av=dMv in regular inverse mode. */
        //   G_SHIFTINVERT, /**< For solving Av=dMv in shift-invert mode. */
        //   G_BUCKLING,    /**< For solving Av=dMv in Buckling mode. */
        //   G_CAYLEY,      /**< For solving Av=dMv in Cayley mode. */
    } EigenSolverMode;

    // Lanczos method (double precision)
    extern void DSAUPD(int* ido, char* bmat, int* n, char* which, int* nev, double* tol, double* resid, int* ncv, double* v, int* ldv, int* iparam,
            int* ipntr, double* workd, double* workl, int* lworkl, int* info);

    // Eigenvalue problem solver for result of dsaupd_
    void DSEUPD(int* rvec, char* howmny, int* select, double* d, double* z, int* ldz, double* sigma, char* bmat, int* n, char* which, int* nev,
            double* tol, double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr, double* workd, double* workl, int* lworkl,
            int* info);

    /**
     * Standard eigenvalue problem solver (A*x = lambda*x) for symmetric matrix of double precision elements. It computes truncated k eigenvalues (and optionally eigenvectors) of a matrix A using
     * ARPACK's Implicitly Restarted Arnoldi (IRA)/Lanczos method in regular mode. BLAS level 2 ops are provided with System optimised BLAS.
     * @param A
     * @param k
     * @param order
     * @param mode
     * @param computeVectors
     * @param symmetric
     * @param eigenValues
     * @param eigenVectors
     * @return 
     */
    int dense_matrix_symmetric_eigen(Matrix* A, int k, EigenOrder order, EigenSolverMode mode, Vector** eigenValues, Matrix** eigenVectors);

    /**
     * Standard eigenvalue problem solver (A*x = lambda*x) for sparse symmetric matrix of double precision elements. It computes truncated k eigenvalues (and optionally eigenvectors) of a symmetric sparse matrix A (in CSR format) using
     * ARPACK's Implicitly Restarted Arnoldi (IRA)/Lanczos method in regular mode. BLAS level 2 ops are provided with System optimised BLAS.
     * @param A
     * @param k
     * @param order
     * @param mode
     * @param computeVectors
     * @param symmetric
     * @param eigenValues
     * @param eigenVectors
     * @return 
     */
    int csr_sparse_matrix_symmetric_eigen(CSRMatrix* A, int k, EigenOrder order, EigenSolverMode mode, Vector** eigenValues, Matrix** eigenVectors);


    /**
     * Standard eigenvalue problem solver (A*x = lambda*x) for symmetric SPARSE matrix of double precision elements. It computes truncated k eigenvalues (and optionally eigenvectors) of a symmetric SPARSE matrix A (in SSS format) using
     * ARPACK's Implicitly Restarted Arnoldi (IRA)/Lanczos method in regular mode. BLAS level 2 ops are provided with System optimised BLAS.
     * 
     * @param A
     * @param k
     * @param order
     * @param mode
     * @param threads
     * @param eigenValues
     * @param eigenVectors
     * @return 
     */
    int sss_sparse_matrix_symmetric_eigen(SSSMatrix* A, int k, EigenOrder order, EigenSolverMode mode, const uint32_t threads, Vector** eigenValues, Matrix** eigenVectors);

    /**
     * Truncated Singular Value Decomposition solver for mxn matrix A. Inspired by arpack-ng/examples/SVD/dsvd.f AND insights gleamed from B. Erichson paper
     * @param A
     * @param k
     * @param requireLeftVectors
     * @param requireRightVectors
     * @param singularValues
     * @param UU
     * @param VV
     * @return 
     */
    int svd(Matrix* A, int k, bool requireLeftVectors, bool requireRightVectors, Vector** singularValues, Matrix** UU, Matrix** VV);

    //    char* get_status(int info);

#ifdef __cplusplus
}
#endif

#endif /* ARPACK_DRIVER_H */
