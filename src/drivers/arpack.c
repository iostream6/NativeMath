/* 
 * File:   arpack_driver.c
 * Author: Ilamah, Osho
 *
 * Created on 06 May 2017, 21:24
 * 
 * 2017.05.06  - Created
 * 2017.05.18  - Introduced ARPACK SVD implementation
 * 2017.05.22  - QC passed
 * 2017.05.24  - Move working buffers from stack to heap to avoid potential stack overflows (see https://stackoverflow.com/a/80113)
 * 2017.06.24  - Added basic sparse matrix support for Eigen routines
 * 2018.04.27  - Added SSS sparse symmetric support
 */

#include <math.h>
#include <string.h>

#include "../headers/external/arpack.h"
#include "../headers/external/xalloc.h"

//Qc'ed (except for CSR)
int sym_eigen_solve(const SSSMatrix* sss_sym, const CSRMatrix* csr_sym, const Matrix* dense_sym, const uint32_t threads, const uint32_t k,
        EigenOrder order, EigenSolverMode mode, Vector** eigenValues, Matrix** eigenVectors) {

    const uint32_t NUM_THREADS = max(1, threads);
    char* solve_message = "\n\tARPACK_DRIVER->SYMMETRIC EIGS with %s input\n";
    //sanity checks - e.g rows = cols
    //
    int ido = 0; // must be zero on first call
    char bmat[] = "I"; // we only support standard eigenvalue problems
    int N = 0;
    if (sss_sym != NULL) {
        N = sss_sym->nrows; //dimensions of the problem
        printf(solve_message, "SSS");
    } else if (dense_sym != NULL) {
        if(dense_sym->symmetric){
            //symmetric internal format not supported by blas backend 
            fprintf(stderr, "Dense symmetric internal format not supported by blas back ends\n");
            return -1;
        }
        N = dense_sym->nrows; //dimensions of the problem
        printf(solve_message, "DENSE");
    } else {
        N = csr_sym->nrows; //dimensions of the problem
        printf(solve_message, "CSR");
    }

    char* which; //which eigen values do you want?
    //
    switch (order) {
        case LA: which = "LA";
            break;
        case SA: which = "SA";
            break;
        case LM: which = "LM";
            break;
        case SM: which = "SM";
            break;
        case BE: which = "BE";
            break;
        default:
            fprintf(stderr, "Invalid eigenvalue order type.\n");
            return -1;
    }
    //
    int nev = max(min(k, N), 0); // number of eigenvalues
    double tol = -1; //use machine precision
    const int DOUBLE_BYTES = sizeof (double);
    //double resid[N]; // will set info = 0, so that this will be randonly initialized by arpack
    double* resid = (double*) aligned_calloc(N, DOUBLE_BYTES, 64); // prefer to do this on the heap as stack might overflow if N is large!!
    const uint32_t ncvGuess = max(20, 3 * nev); // arpack docs says should be >= 2*NEV. For small problems (i.e. k), we want it to be at least reasonable, say 20
    int ncv = ncvGuess > N ? N : ncvGuess; //number of lanczos column vectors
    // More sanity checks
    if (nev > N) {
        fprintf(stderr, "Condition \"nev <= N\" is not met.\n");
        return -1;
    }
    if (nev + 1 > ncv) {
        fprintf(stderr, "Condition \"nev+1 <= ncv\" is not met.\n");
        return -2;
    }
    //
    //double V[N * ncv]; //The NCV columns of V contain the Lanczos basis vectors.
    double* V = (double*) aligned_calloc(N * ncv, DOUBLE_BYTES, 64); // prefer to do this on the heap as stack might overflow if it is large!!
    int ldv = N; // Leading dimension of V exactly as declared in the calling program
    int iparam[11] = {0}; // will initialize all to zero http://stackoverflow.com/a/201116
    // all iparam elements are zero, lets set those that we want to be non zero:
    iparam[0] = 1; /* Specifies the shift strategy (1->exact). */
    iparam[2] = 10 * N; /* Maximum number of iterations. */
    iparam[3] = 1; /* NB blocksize in recurrence. ARPACK mentions only 1 is supported currently. */
    switch (mode) {
        case I_REGULAR:
            iparam[6] = 1;
            break;
            //        case I_SHIFTINVERT:
            //            iparam[6] = 3;
            //            break;
        default:
            fprintf(stderr, "Invalid Eigs mode.\n");
            return -1;
    }
    //    
    int ipntr[11] = {0}; /* Output. Init to zero */
    //
    //double workd[3 * N]; //double workd[3 * N] = {0}; doesnt work in C:: variable-sized object may not be initialized, so use memset (requires string.h) and note that it works by bytes (see http://stackoverflow.com/a/6816467)
    //memset(workd, 0, sizeof workd);
    double* workd = (double*) aligned_calloc(3 * N, DOUBLE_BYTES, 64); // prefer to do this on the heap as stack might overflow. Added bonus is we get faster initialization
    int lworkl = ncv * (ncv + 8);
    //double workl[lworkl];
    //memset(workl, 0, sizeof workl);
    double* workl = (double*) aligned_calloc(lworkl, DOUBLE_BYTES, 64); // prefer to do this on the heap as stack might overflow if buffer is large!!
    int info = 0;
    do {
        DSAUPD(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
        //
        switch (ido) {
            case -1:
            case 1:
                // compute  Y = A*x  where
                //IPNTR(1) is the pointer into WORKD for X, IPNTR(2) is the pointer into WORKD for Y.
                //note that we subtract 1 to get back in line to C zero based array indices c.f. with fortran (1 based)
                if (sss_sym != NULL) {
                    //SS sparse symmetric matrix
                    sparse_sym_matrix_vector_mult(sss_sym, &(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]), NUM_THREADS);
                } else if (dense_sym != NULL) {
                    matrix_vector_mult(dense_sym, &(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
                } else {
                    //CSR sparse  symmetric matrix
                    sparse_matrix_vector_mult(csr_sym, &(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
                }
                break;
                // case 2://IDO =  2: compute  Y = B * X  should not happen as we only support standard eigen problems
                // case 3://IDO =  3: compute the IPARAM(8) shifts should not happen as we dont use user shifts
        }
    } while (ido != 99); /* Finish condition. */
    //print_lanczos_status(info);
    /* Normal exit? */
    if (info != 0) {
        return -3; //abnormal exit
    }

    // No fatal errors or abnormal situation occurred so we can post-process all nev eigenvalues using DSEUPD, without selection. 
    // Corresponding eigenvectors may also be computed now if desired.  (indicated by rvec = .true.)
    int rvec = 1;
    char howmny[] = "A"; //home many of the eigenvalues' corresponding eigenvectors do we want - we want for all eigenvalues that have been computed
    int select[ncv]; // as howmny is A, this will be used as a workspace for reordering the Ritz values
    // ARPACK orders the eigen values back to front so if mode is LM and a >= b >=c >= d are eigenvalues, ARPACK will arrange results as : d, c, b, a, We dont want this, so we will flip as required. 
    Vector *d = create_vector(nev); //fortran space eigenvalues (ARPACK OUTPUT), these are returned by ARPACK in ascending order
    Matrix *Z = create_matrix(N, nev, false); //  fortran space eigenvectors (ARPACK OUTPUT), Vectors are COLUMN MAJOR
    int ldz = N;
    double sigma = 0; //should not be references as iparam[6] = 1;
    DSEUPD(&rvec, howmny, select, d->data, Z->data, &ldz, &sigma, bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    //free memory
    aligned_free((void*) resid);
    aligned_free((void*) V);
    aligned_free((void*) workd);
    aligned_free((void*) workl);
    //
    //print_lanczos_status(info);
    if (info < 0) {
        return -4; //post processing failure
    }
    // extract (and optionally reorder) the ARPACK eigenpairs
    switch (order) {
        case LA: which = "LA"; //results will be same as MATLAB
        case LM: which = "LM"; //results will be same as MATLAB
            //arpack results are always in ascending order, so if the user asked for LA or LM, we flip arpack results back to front
            // reordering of the ARPACK eigenvalues results
            *eigenValues = create_vector(nev); //C space eigenvalues
            int evalIndex = nev - 1; //fortran space element indexes
            for (int vector = 0; vector < nev;) {
                (*eigenValues)->data[vector++] = d->data[evalIndex--];
            }
            //delete & free memory
            delete_vector(d);
            *eigenVectors = create_matrix(N, nev, false); //C space eigenvectors
#pragma omp parallel for
            for (int vector = 0; vector < nev; vector++) {
                int vStart = (nev - vector - 1) * N;
                const int vEnd = vStart + N;
                int ij = vector; //C space matrix element dereferencer
                for (; vStart < vEnd;) {
                    //read fortran column data elements sequentially and place in sequential elements in the c space matrix
                    (*eigenVectors)->data[ij] = Z->data[vStart++];
                    ij += nev;
                }
            }
            //delete & free memory
            delete_matrix(Z);
            break;
        default:
            //        case SA: which = "SA"; //results will be same as MATLAB
            //        case SM: which = "SM"; //results will be FLIPPED c.f. MATLAB
            //        case BE: which = "BE"; //results will be same as MATLAB
            //user requested order and than of ARPACK will match in this case
            (*eigenValues) = d; // *d must not be freed, dangling pointer otherwise!
            //
            *eigenVectors = create_matrix(N, nev, false); //C space eigenvectors
#pragma omp parallel for
            for (int vector = 0; vector < nev; vector++) {
                int vStart = N*vector;
                const int vEnd = vStart + N;
                int ij = vector; //C space matrix element dereferencer
                for (; vStart < vEnd;) {
                    //read fortran column data elements sequentially and place in sequential elements in the c space matrix
                    (*eigenVectors)->data[ij] = Z->data[vStart++];
                    ij += nev;
                }
            }
            //delete & free memory
            delete_matrix(Z);
            break;
    }
    return (EXIT_SUCCESS);
}


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
//Qc'ed
int dense_matrix_symmetric_eigen(Matrix* A, int k, EigenOrder order, EigenSolverMode mode, Vector** eigenValues, Matrix** eigenVectors) {
    return sym_eigen_solve(NULL, NULL, A, 1, k, order, mode, eigenValues, eigenVectors);
}

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
int csr_sparse_matrix_symmetric_eigen(CSRMatrix* A, int k, EigenOrder order, EigenSolverMode mode, Vector** eigenValues, Matrix** eigenVectors) {
    return sym_eigen_solve(NULL, A, NULL, 1, k, order, mode, eigenValues, eigenVectors);
}

/**
 * Standard eigenvalue problem solver (A*x = lambda*x) for symmetric SPARSE matrix of double precision elements. It computes truncated k eigenvalues (and optionally eigenvectors) of a symmetric SPARSE matrix A (in SSS format) using
 * ARPACK's Implicitly Restarted Arnoldi (IRA)/Lanczos method in regular mode. BLAS level 2 ops are provided with System optimised BLAS.
 * 
 * @param A
 * @param k
 * @param order
 * @param mode
 * @param eigenValues
 * @param eigenVectors
 * @return 
 */
//Qc'ed
int sss_sparse_matrix_symmetric_eigen(SSSMatrix* A, int k, EigenOrder order, EigenSolverMode mode, const uint32_t threads, Vector** eigenValues, Matrix** eigenVectors) {
    return sym_eigen_solve(A, NULL, NULL, threads, k, order, mode, eigenValues, eigenVectors);
}

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
//Qc'ed
int svd(Matrix* A, int k, bool requireLeftVectors, bool requireRightVectors, Vector** singularValues, Matrix** UU, Matrix** VV) {
    int N = A->ncols; //dimensions of the problem
    int M = A->nrows;
    bool solvingRight = M >= N;
    EigenOrder order = LA;
    EigenSolverMode mode = I_REGULAR;
    //unlike ARAPACK example where they solve (A'*A)*v = sigma*v by phased multiplication (first Ax, then A'(Ax)), I do premultiplication once
    //this might be a bit costly, or not, but if m != n it will be faster since we get a smaller matrix
    //it also allows us to generalise for m<n, by using a different formulation (A*A')x = sigma*x instead of (A'*A)x = sigma*x
    Matrix *B;
    //B will be symmetric as B = At*A or A*At 
    if (solvingRight) {
        B = create_matrix(N, N, false /* even though it is symmetric, we dont want compaction as blas will not support this */); //(A'*A)
        matrix_transpose_matrix_mult(A, A, B);
    } else {
        B = create_matrix(M, M, false /* even though it is symmetric, we dont want compaction as blas will not support this */); //(A*A')
        matrix_matrix_transpose_mult(A, A, B); //(A*A')
    }
    const bool requireVectors = requireLeftVectors || requireRightVectors;
    int eigenStatus = dense_matrix_symmetric_eigen(B, k, order, mode, singularValues, solvingRight ? VV : UU);
    //free B
    delete_matrix(B);
    if (eigenStatus == 0) {
        const int nev = (*singularValues)->nrows;
        for (int vector = 0; vector < nev; vector++) {
            (*singularValues)->data[vector] = sqrt((*singularValues)->data[vector]);
        }
        if (requireVectors) {
            if (solvingRight && requireLeftVectors) {
                *UU = create_matrix(M, nev, false); //C space singular vectors
                //U = Av/sigma (see ARPACK-NG SVD example), see my notes
                matrix_matrix_mult(A, *VV, *UU);
#pragma omp parallel for        
                for (int vector = 0; vector < nev; vector++) {
                    //printf("Number of Threads::%d", omp_get_num_threads());
                    const double norm = 1.0 / (*singularValues)->data[vector]; //normalize each vector by the singular vals (also may be norm2 of the corresponding U vect)
                    int address = vector;
                    for (int element = 0; element < M; element++) {
                        (*UU)->data[address] = (*UU)->data[address] * norm;
                        address += nev;
                    }
                }
            } else if (!solvingRight && requireRightVectors) {
                *VV = create_matrix(N, nev, false); //C space singular vectors
                matrix_transpose_matrix_mult(A, *UU, *VV);
#pragma omp parallel for        
                for (int vector = 0; vector < nev; vector++) {
                    //printf("Number of Threads::%d", omp_get_num_threads());
                    const double norm = 1.0 / (*singularValues)->data[vector]; //normalize each vector by the singular vals (also may be norm2 of the corresponding V vect)
                    int address = vector;
                    for (int element = 0; element < N; element++) {
                        (*VV)->data[address] = (*VV)->data[address] * norm;
                        address += nev;
                    }
                }
            }
        }
    }
    return eigenStatus;
}

