/*
 * 2018.04.15  - Based on reorg of legacy C codes written in 2017
 * 2018.04.24  - Implemented multithreaded sparse_sym_matrix_vector_mult using "Effective ranges" method in 2013 paper by Gkountouvas et al
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


#include "../headers/external/xblas.h"
#include "../headers/external/xalloc.h"

int SYSMATH_INIT_ = -1;
int sys_num_threads = 0;

void sparse_sym_matrix_vector_mult_serial();

//TODO specialize dense matrix blas functions (e.g. mat-vec mult) for symmetric internal memory saving format

/**
 * Specifies the number of threads to be used by BLAS backend 
 * @param num_threads
 */
void set_system_num_threads(const uint32_t num_threads) {
    SYSMATH_INIT_ = -1; //ensures any changes to numthreads by client has a chance to be pushed to the lower levels
    sys_num_threads = max(0, num_threads);
}

/**
 * Initializes diagonal matrix from vector data 
 * @param D
 * @param v
 */
//QC'ed

void initialize_diagonal_matrix(Matrix* D, Vector* v) {
    //boundary safety measures, OI
    uint32_t l = min(D->nrows, D->ncols);
    l = min(l, v->nrows);
    for (uint32_t i = 0; i < l; i++) {
        set_matrix_element(D, i, i, v->data[i]);
    }
}

/**
 * Computes y = M*x ; row major 
 * @param M
 * @param x
 * @param y
 */
//QC'ed
void matrix_vector_mult(const Matrix* M, const double* x, double* y) {
    ///\/\/\/\if symmetrix do special
    const double alpha = 1.0, beta = 0.0;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, M->nrows, M->ncols, alpha, M->data, M->ncols, x, 1, beta, y, 1);
}

/**
 * A basic routine to compute the sparse matrix - dense vector product:  y = A*x, row major. Inspired by http://www.netlib.org/utk/people/JackDongarra/etemplates/node382.html. This is mem 
 * bandwidth constrained for numerically symmetric cases, better to use SSS in such instances
 * @param M  the input sparse matrix in CSR3 format
 * @param x  the input dense vector
 * @param y  a dense vector of the results
 */
//QC'ed

void sparse_matrix_vector_mult(const CSRMatrix* M, double * x, double * y) {
    //assume client has checked dimensions are sensible for mults
    const uint32_t ROWS = M->nrows;
    const double* values = M->values;
#pragma omp parallel for
    for (uint32_t i = 0; i < ROWS; i++) {
        double rowSum = 0;
        const int start = M->row_ptr[i];
        const int end = M->row_ptr[i + 1];
        for (uint32_t j = start; j < end; j++) {
            rowSum += values[j] * x[M->col_ind[j]];
        }
        y[i] = rowSum;
    }
}

/**
 * Performs  a symmetric sparse matrix to dense vector mulplication operation: y = M*x ; It is assumed that the output vector has been allocated by the client 
 * @param M Input symmetric sparse matrix in SSS format
 * @param x Input dense vector
 * @param y Output vector
 */
void sparse_sym_matrix_vector_mult(const SSSMatrix* M, const double* x, double* y, const uint32_t threads) {
    const uint32_t NUM_THREADS = max(1, threads);
    if (NUM_THREADS == 1) {
        sparse_sym_matrix_vector_mult_serial(M, x, y, M->nrows);
        return;
    }
    double* slave_buffers[NUM_THREADS]; //an array of double pointers  https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/
    uint32_t starts[NUM_THREADS]; //, ends[NUM_THREADS];
    const uint32_t ROWS_PER_THREAD = max((M->nrows / NUM_THREADS), 1);
    const uint32_t NROWS = M->nrows;
#pragma omp parallel for num_threads(NUM_THREADS)
    for (uint32_t thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
        //prepare structures and partition input matrix - TODO here we can inject more intelligent heuristics for partitioning - at the moment, we are dividing rows approximately 
        const uint32_t start = ROWS_PER_THREAD * thread_id;
        const uint32_t end = thread_id == (NUM_THREADS - 1) ? M->nrows : start + ROWS_PER_THREAD; /* the last thread needs to do all remaining ones */
        //
        starts[thread_id] = start;
        //for all but the zeroth thread, create rightly sized buffers for effective range
        if (thread_id != 0) {
            slave_buffers[thread_id] = (double*) aligned_calloc(start /* - 0 */, sizeof (double), 64);
        } /*else {
            // No adverse side effect as zeroth thread does not require thread local buffer 
            slave_buffers[thread_id] = NULL;
        }*/
        //printf("\n Buffer index = %d  :: Thread Index = %d\n", thread_id, omp_get_thread_num());
        double* slave_buffer = slave_buffers[thread_id];
        // perform the thread calcs
        for (uint32_t row = start; row < end; row++) {
            y[row] = M->dvalues[row] * x[row];
            const uint32_t next_row_pointer = M->row_ptr[row + 1];
            for (uint32_t row_pointer = M->row_ptr[row]; row_pointer < next_row_pointer; row_pointer++) {
                const uint32_t col_index = M->col_ind[row_pointer];
                const double value = M->values[row_pointer];
                y[row] += (value * x[col_index]);
                const double transpose_element_product = value * x[row];
                if (col_index < start) {
                    slave_buffer[col_index] += transpose_element_product;
                } else {
                    y[col_index] += transpose_element_product;
                }
            }
        }
    }
#pragma omp parallel for num_threads( NUM_THREADS) //reduction step TODO in PARALLEL
    for (uint32_t row = 0; row < NROWS; row++) {
        //printf("\n XXXX buffer index = 0  :: Thread Index = %d\n",  omp_get_thread_num());
        double partial_sum = 0;
        for (uint32_t thread_id = 1; thread_id < NUM_THREADS; thread_id++) {
            if (starts[thread_id] > row) {
                partial_sum += slave_buffers[thread_id][row];
            }
        }
        y[row] += partial_sum;
    }
}

//QC'ed
void sparse_sym_matrix_vector_mult_serial(const SSSMatrix* M, const double* x, double* y, const uint32_t NROWS) {
    for (uint32_t row = 0; row < NROWS; row++) {
        y[row] = M->dvalues[row] * x[row];
        const uint32_t start = M->row_ptr[row], end = M->row_ptr[row + 1];
        for (uint32_t j = start; j < end; j++) {
            const uint32_t col = M->col_ind[j];
            const double value = M->values[j];
            y[row] += value * x[col];
            y[col] += value * x[row];
        }
    }
}

/**
 * Computes C = A*B ; assuming ROW major format 
 * @param A
 * @param B
 * @param C
 */
//QC'ed
void matrix_matrix_mult(Matrix* A, Matrix* B, Matrix* C) {
    ///\/\/\/\if symmetrix do special
    const double alpha = 1.0, beta = 0.0;
    //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, k, B, n, beta, C, n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->nrows, B->ncols, A->ncols, alpha, A->data, A->ncols, B->data, B->ncols, beta, C->data, B->ncols);
}

/**
 * Computes C = A^T*B ; assuming ROW major format  
 * @param A
 * @param B
 * @param C
 */
//QC'ed
void matrix_transpose_matrix_mult(Matrix* A, Matrix* B, Matrix* C) {
    ///\/\/\/\if symmetrix do special
    const double alpha = 1.0, beta = 0.0;
    //cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, A, k, B, n, beta, C, n);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, A->ncols, B->ncols, A->nrows, alpha, A->data, A->ncols, B->data, B->ncols, beta, C->data, B->ncols);
}

/**
 * Computes C = A*B^T ; assuming ROW major format  
 * @param A
 * @param B
 * @param C
 */
//QC'ed
void matrix_matrix_transpose_mult(Matrix* A, Matrix* B, Matrix* C) {
    ///\/\/\/\if symmetrix do special
    const double alpha = 1.0, beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, A->nrows, B->nrows, A->ncols, alpha, A->data, A->ncols, B->data, B->ncols, beta, C->data, C->ncols);
}

/**
 * Computes Q from [Q,R] = qr(M,'0') compact QR factorization, assuming ROW major format  
 * @param M mxn matrix
 * @param Q mxn matrix
 */
void QR_factorization_getQ(Matrix* M, Matrix* Q) {
    const uint32_t m = M->nrows, n = M->ncols, k = min(m, n);
    matrix_copy(Q, M);
    Vector* tau = create_vector(k);
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, Q->data, n, tau->data);
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, n, Q->data, n, tau->data);
    delete_vector(tau);
}

/**
 * Computes [Q,R] = qr(M,'0'), or  the compact QR factorization 
 * @param M mxn matrix
 * @param Q mxn matrix
 * @param R min(m,n) x min(m,n)
 */
void compact_QR_factorization(Matrix* M, Matrix* Q, Matrix* R) {
    const uint32_t m = M->nrows, n = M->ncols;
    const uint32_t k = min(m, n);
    //printf("->QR with m = %d, n = %d, k = %d\n", m, n, k);
    Matrix* R_full = create_matrix(m, n, M->symmetric);
    matrix_copy(R_full, M);
    Vector* tau = create_vector(k);
    // get R
    //printf("get R..\n");
    //LAPACKE_dgeqrf(CblasColMajor, m, n, R_full->d, n, tau->d);
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, R_full->nrows, R_full->ncols, R_full->data, R_full->ncols, tau->data);
    for (uint32_t i = 0, j = 0; i < k; i++) {
        for (j = 0; j < k; j++) {
            if (j >= i) {
                set_matrix_element(R, i, j, get_matrix_element(R_full, i, j));
            }
        }
    }
    // get Q
    matrix_copy(Q, R_full);
    //printf("dorgqr..\n");
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, Q->nrows, Q->ncols, min(Q->ncols, Q->nrows), Q->data, Q->ncols, tau->data);
    // clean up
    delete_matrix(R_full);
    delete_vector(tau);
}

/**
 * Copies elements from the source to the destination matrix
 * @param destination
 * @param source
 */
//QC'ed

void matrix_copy(Matrix* destination, Matrix* source) {
    int elements = min(((source->nrows)*(source->ncols)), ((destination->nrows)*(destination->ncols))); //safety - hopefully the caller knows what it is doing - haha! 
    for (uint32_t i = 0; i < elements; i++) {
        destination->data[i] = source->data[i];
    }
}

/**
 * Extracts a column of the provided matrix into a vector
 * @param M the matrix
 * @param j the column to extract
 * @param vector the vector to polulate
 */
//QC'ed
void get_matrix_col(Matrix* M, const uint32_t j, Vector* vector) {
    for (uint32_t i = 0; i < M->nrows; i++) {
        vector->data[i] = get_matrix_element(M, i, j);
    }
}

/**
 * Computes the SVD of A = U*Σ*V^T for real routines, assuming row major format
 * @param M
 * @param U
 * @param svals vector containing the diagonal elements of Σ, in descending order of magnitude.
 * @param Vt
 */
void singular_value_decomposition(Matrix* M, Matrix *U, Vector* svals /*Matrix *S*/, Matrix *Vt) {
    const uint32_t m = M->nrows, n = M->ncols;
    const uint32_t k = min(m, n);
    Vector* work = create_vector(2 * max(3 * min(m, n) + max(m, n), 5 * min(m, n)));
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', m, n, M->data, n, svals->data, U->data, k, Vt->data, n, work->data);
    delete_vector(work);
}

//int compute_evals_and_evecs_of_symm_matrix(Matrix *S, vec *evals);

/**
 * Constitutes a matrix P, using the SVD factors U, S and V. This is useful for validating SVD accuracy.
 * @param U
 * @param S
 * @param V
 * @param P
 */
void form_svd_product_matrix(Matrix* U, Matrix* S, Matrix* V, Matrix* P) {
    const uint32_t k = S->nrows, n = P->ncols;
    Matrix * SVt = create_matrix(k, n, false);
    // form SVt = S*V^T
    matrix_matrix_transpose_mult(S, V, SVt);
    // form P = U*S*V^T
    matrix_matrix_mult(U, SVt, P);
    delete_matrix(SVt);
}

/**
 * Computes the frobenius norm of a matrix. This is useful for validating decomposition accuracy.
 * @param M
 * @return the frobenius norm
 */

double get_matrix_frobenius_norm(Matrix* M) {
    const uint32_t M_times_N = (M->nrows)*(M->ncols);
    double aij_squared = 0, aij = 0;
    for (uint32_t i = 0; i < M_times_N; i++) {
        aij = M->data[i];
        aij_squared += aij*aij;
    }
    return sqrt(aij_squared);
}

//double get_matrix_frobenius_norm0(Matrix* M) {
//    int i;
//    double val, normval = 0;
//#pragma omp parallel shared(M,normval) private(i,val) 
//    {
//#pragma omp for reduction(+:normval)
//        for (i = 0; i < ((M->nrows)*(M->ncols)); i++) {
//            val = M->data[i];
//            normval += val*val;
//        }
//    }
//    return sqrt(normval);
//}

#ifdef __SYSMATH_MKL__
//Using MKL backend

void mkl_init() {
    //the code is only really needed for linux
    if (SYSMATH_INIT_ < 0) {
        printf("\n\n\n******************************************\n******************************************\n");
        printf("MKL thread layer initialization performed\n");
#ifdef  __linux__
        mkl_set_threading_layer(MKL_THREADING_GNU);
#endif  /* __linux__ */ 
        if (sys_num_threads > 0) {
            printf("Parallel domains set by calling client\n");
            mkl_set_dynamic(0);
        }
        mkl_set_num_threads(sys_num_threads);
        printf("System Math MKL calls will use %d threads\n", mkl_get_max_threads());
        printf("System Math MKL threading is %s\n\n", mkl_get_dynamic() == 0 ? "FIXED" : "DYNAMIC");
        SYSMATH_INIT_ = 1;
    }
}

/**
 * Initialises the elements of the provided matrix with normally distributed values of zero mean and sd of 1.0.
 * @param M the matrix to randomly initialise
 * //TODO:: OI query, why not just call vdrnggaussian directly onthe matrices buffer - avoids the need for alocating extra membuffer as well as omp for loop
 */
void mkl_initialize_random_matrix(Matrix *M) {
    const uint32_t m = M->nrows, n = M->ncols;
    const float a = 0.0, sigma = 1.0;
    int N = m*n;
    float* r = (float*) mkl_calloc(N, sizeof (double), 64);
    VSLStreamStatePtr stream;
    vslNewStream(&stream, BRNG, time(NULL));
    //vslNewStream( &stream, BRNG,  SEED );
    vsRngGaussian(METHOD, stream, N, r, a, sigma);
    // read and set elements
    for (uint32_t i = 0; i < N; i++) {
        M->data[i] = r[i];
    }
    mkl_free(r);
}

#else
//Using OPENBLAS backend - bring in the openblas headers and OMP for multi-threading support

/**
 * Initialises the elements of the provided matrix with normally distributed values of zero mean and sd of 1.0 using the zrng pseudo random number provider. 
 * This function will be used by non-MKL backend as implementation of initialize_random_matrix(Matrix *M)
 * 
 * @param M the matrix to randomly initialise
 */
void zrng_initialize_random_matrix(Matrix* M) {
    const uint32_t seed = (uint32_t) time(NULL);
    // so time(null) will generate 8 byte integer contaiming number of sceconds since epoch. 
    // The narrowing cast above will be safe until the 22nd century - I would be long dead by then, at the earliest (4 byte unsigned int will be able to take up to 4294967295 seconds
    // see https://stackoverflow.com/questions/471248
    const int N = M->ncols * M->nrows;
    rngDGaussian(seed, N, M->data, 0, 1.0);
}

void openblas_init() {
    if (SYSMATH_INIT_ < 0 && sys_num_threads > 0) {
        printf("\n\n\n******************************************\n******************************************\n");
        printf("OPENBLAS thread layer initialization performed\n");
        omp_set_num_threads(sys_num_threads);
        printf("System Math BLAS calls will use %d threads\n\n", sys_num_threads);
        SYSMATH_INIT_ = 1;
    }
}

#endif  //end __SYSMATH_MKL__ ifdef_else_endif block