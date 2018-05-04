/*
 * 20.04.2018
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "../../src/headers/external/xalloc.h"
#include "../../src/headers/external/matvec.h"
#include "../../src/headers/external/xblas.h"
#include "../../src/headers/external/arpack.h"

#include "../../src/headers/external/xmmio.h"


#define WINDOW_SIZE 10
#define SMALL_MAT_SIZE 20
#define BASE_DIR "X:\\delete\\x\\"


int xblas_test();

int arpack_sym_eigen_test();
int arpack_svd_test();

int main() {
    time_t t;
    /* Intializes random number generator */
    srand((unsigned) time(&t));

    //xblas_test();
    //arpack_sym_eigen_test();
    
    arpack_svd_test();

    printf("\n\nPlease type a character to exit::");
    getchar();
    return EXIT_SUCCESS;

}

int xblas_test() {
    sysmath_init(); //required for GNU/Linux with MKL backend to force use of OPEN MP and avoid shared object RTL conflicts with Intel OPENMP
    printf("\n\n======== initialize diagonal matrix ========\n");
    Vector* x = create_vector(16);
    Matrix* D = create_matrix(x->nrows, x->nrows, false);
    for (uint32_t row = 0; row < x->nrows; row++) {
        //x->data[row] = row + 1/*row / (2.0)*/;
        //printf("%lg ", x->data[row]);
        const int val = (rand() % 200) - 101;
        x->data[row] = val;
    }
    initialize_diagonal_matrix(D, x);
    for (uint32_t row = 0; row < x->nrows; row++) {
        for (uint32_t col = 0; col < x->nrows; col++) {
            if (col == row) {
                if (get_matrix_element(D, row, col) != x->data[row]) {
                    printf("***** !!!!  initialize diagonal matrix FAILED!!!\n");
                    delete_vector(x);
                    delete_matrix(D);
                    return -1;
                }
            } else {
                if (get_matrix_element(D, row, col) != 0) {
                    printf("***** !!!!  initialize diagonal matrix FAILED!!!\n");
                    delete_vector(x);
                    delete_matrix(D);
                    return -1;
                }
            }
        }
    }
    printf("***** !!!!  initialize diagonal matrix SUCCESS!!!\n");

    printf("\n\n======== copy dense matrix ========\n");
    Matrix* Dcopy = create_matrix(x->nrows, x->nrows, false);
    matrix_copy(Dcopy, D);
    for (uint32_t row = 0; row < x->nrows; row++) {
        for (uint32_t col = 0; col < x->nrows; col++) {
            if (col == row) {
                if (get_matrix_element(Dcopy, row, col) != x->data[row]) {
                    printf("***** !!!!  copy dense matrix FAILED!!!\n");
                    delete_vector(x);
                    delete_matrix(D);
                    delete_matrix(Dcopy);
                    return -1;
                }
            } else {
                if (get_matrix_element(Dcopy, row, col) != 0) {
                    printf("***** !!!!  copy dense matrix FAILED!!!\n");
                    delete_vector(x);
                    delete_matrix(D);
                    delete_matrix(Dcopy);
                    return -1;
                }
            }
        }
    }
    printf("***** !!!!  copy dense matrix SUCCESS!!!\n");

    delete_vector(x);
    delete_matrix(D);
    delete_matrix(Dcopy);


    printf("\n\n======== Mat-vect  mult ========\n");
    //See https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-matrix-storage-formats
    // https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-blas-csr-matrix-storage-format
    Matrix* M = create_matrix(5, 5, false);
    double* data = M->data;
    data[0] = 1;
    data[1] = -1;
    data[3] = -3;
    data[5] = -2;
    data[6] = 5;
    data[12] = 4;
    data[13] = 6;
    data[14] = 4;
    data[15] = -4;
    data[17] = 2;
    data[18] = 7;
    data[21] = 8;
    data[24] = -5;

    CSRMatrix* csrM = csr_from_matrix(M);
    print_csr_matrix(csrM, stdout, 15);
    delete_matrix(M);

    double* yy = (double*) aligned_calloc(5, sizeof (double), 64);
    double* xx = (double*) aligned_calloc(5, sizeof (double), 64);
    double* w = (double*) aligned_calloc(5, sizeof (double), 64);
    //
    xx[0] = 1;
    xx[1] = 2;
    xx[2] = 3;
    xx[3] = 4;
    xx[4] = 5;
    //
    yy[0] = -13;
    yy[1] = 8;
    yy[2] = 56;
    yy[3] = 30;
    yy[4] = -9;

    sparse_matrix_vector_mult(csrM, xx, w);
    printf("\n\tSparse non-sym matvec::");
    for (int v = 0; v < 5; v++) {
        if (yy[v] - w[v] > 0.0000001) {
            printf("***** !!!!  Sparse non-sym matrix-vector mult FAILED!!!\n");
            aligned_free(xx);
            aligned_free(yy);
            aligned_free(w);

            delete_csrmatrix(csrM);
            return -1;
        }
        printf("%0.1f  ", w[v]);
    }
    printf(":: SUCCESS!!!");
    aligned_free(xx);
    aligned_free(yy);
    aligned_free(w);

    delete_csrmatrix(csrM);

    printf("\n\tSparse SYM matvec::");
    char* sparse_sym_filename = BASE_DIR "Ssym_mat_100.mtx";
    char* result_filename = BASE_DIR "Y_mat_100.mtx";
    SSSMatrix* S = load_mm_ssmatrix(sparse_sym_filename); //will still work if diagonal elements are zero
    print_sss_matrix(S, stdout, 20);

    M = load_mm_matrix(result_filename);
    if (M == NULL) {
        printf("***** !!!!  Sparse symmetric matrix-vector mult FAILED on load!!!\n");
        return -1;
    }

    x = create_vector(M->nrows);
    get_matrix_col(M, 0, x);

    yy = (double*) aligned_calloc(M->nrows, sizeof (double), 64);

    sparse_sym_matrix_vector_mult(S, x->data, yy, 1);

    for (int v = 0; v < M->nrows; v++) {
        if (yy[v] - get_matrix_element(M, v, 2) > 0.0000001) {
            printf("***** !!!!  Sparse symmetric matrix-vector mult FAILED!!!\n");
            aligned_free(yy);
            delete_vector(x);
            delete_sssmatrix(S);
            delete_matrix(M);
            return -1;
        }
    }
    printf("\n\tSparse symmetric matvec:: SUCCESS!!!");

    aligned_free(yy);
    delete_vector(x);
    delete_sssmatrix(S);
    delete_matrix(M);


    printf("\n\tDense mat-vec & mat-mat::");
    char* denseA_filename = BASE_DIR "A_dns_dense_non_sym_mat_10.mtx";
    char* denseB_filename = BASE_DIR "B_dns_dense_non_sym_mat_10.mtx";
    char* denseC_filename = BASE_DIR "C_dns_dense_non_sym_mat_10.mtx";
    char* denseP_filename = BASE_DIR "P_dns_dense_non_sym_mat_10.mtx";

    Matrix* A = load_mm_matrix(denseA_filename);
    Matrix* B = load_mm_matrix(denseB_filename);
    Matrix* C = load_mm_matrix(denseC_filename);
    Matrix* P = load_mm_matrix(denseP_filename);

    if (A == NULL || B == NULL || C == NULL || P == NULL) {
        printf("***** !!!!  dense non symmetric matrix load FAILED !!!\n");
        return -1;
    }
    //perform matrix vect mult
    yy = (double*) aligned_calloc(A->nrows, sizeof (double), 64);
    matrix_vector_mult(M, C->data, yy);
    for (int v = 0; v < A->nrows; v++) {
        if (yy[v] - get_matrix_element(P, v, 30) > 0.0000001) {
            printf("\n\t***** !!!!  Dense non-symmetric matrix-vector mult FAILED!!!\n");
            aligned_free(yy);
            delete_matrix(A);
            delete_matrix(B);
            delete_matrix(C);
            delete_matrix(P);
            return -1;
        }
    }
    aligned_free(yy);
    printf("\n\tDense non-symmetric matvec:: SUCCESS!!!");


    //perform matrix matrix mult
    M = create_matrix(A->nrows, B->ncols, 0);

    matrix_matrix_mult(A, B, M);
    for (int i = 0; i < A->nrows; i++) {
        for (int j = 0; j < B->ncols; j++) {
            if (get_matrix_element(M, i, j) - get_matrix_element(P, i, j) > 0.0000001) {
                printf("\n\t***** !!!!  Dense non-symmetric matrix-matrix mult FAILED!!!\n");
                delete_matrix(A);
                delete_matrix(B);
                delete_matrix(C);
                delete_matrix(P);
                delete_matrix(M);
                return -1;
            }
        }
    }
    delete_matrix(M);
    printf("\n\tDense non-symmetric matrix-matrix:: SUCCESS!!!");

    //perform matrixT matrix mult
    M = create_matrix(A->ncols, B->ncols, 0);
    matrix_transpose_matrix_mult(A, B, M);
    for (int i = 0; i < A->nrows; i++) {
        for (int j = 0; j < B->ncols; j++) {
            if (get_matrix_element(M, i, j) - get_matrix_element(P, i, (B->ncols + j)) > 0.0000001) {
                printf("\n\t***** !!!!  Dense non-symmetric matrix-transpose-matrix mult FAILED!!!\n");
                delete_matrix(A);
                delete_matrix(B);
                delete_matrix(C);
                delete_matrix(P);
                delete_matrix(M);
                return -1;
            }
        }
    }
    delete_matrix(M);
    printf("\n\tDense non-symmetric matrix-transpose-matrix:: SUCCESS!!!");

    //perform matrix matrixT mult
    M = create_matrix(A->nrows, B->nrows, 0);

    matrix_matrix_transpose_mult(A, B, M);
    for (int i = 0; i < A->nrows; i++) {
        for (int j = 0; j < B->ncols; j++) {
            if (get_matrix_element(M, i, j) - get_matrix_element(P, i, (B->ncols + B->ncols + j)) > 0.0000001) {
                printf("\n\t***** !!!!  Dense non-symmetric matrix-matrix-transpose mult FAILED!!!\n");
                delete_matrix(A);
                delete_matrix(B);
                delete_matrix(C);
                delete_matrix(P);
                delete_matrix(M);
                return -1;
            }
        }
    }
    delete_matrix(A);
    delete_matrix(B);
    delete_matrix(C);
    delete_matrix(P);
    delete_matrix(M);
    printf("\n\tDense non-symmetric matrix-matrix-transpose:: SUCCESS!!!");
    return 0;
}

int arpack_sym_eigen_test() {
    sysmath_init(); //required for GNU/Linux with MKL backend to force use of OPEN MP and avoid shared object RTL conflicts with Intel OPENMP
    printf("\n\n======== Symmetric matrix eigs :: ========\n");
    int k = 10; // the rank we want
    char* sparse_sym_filename = BASE_DIR "Ssym_mat_100.mtx";
    char* sym_filename = BASE_DIR "Asym_mat_100.mtx";
    SSSMatrix* S = load_mm_ssmatrix(sparse_sym_filename); //will still work if diagonal element is zero
    if (S == NULL) {
        printf("***** !!!!  Symmetric matrix eigs FAILED on load!!!\n");
        return -1;
    }

    Matrix* Asym = load_mm_matrix(sym_filename);
    
    Matrix* A = 0;

    if (Asym == NULL) {
        delete_sssmatrix(S);
        printf("***** !!!!  Symmetric matrix eigs FAILED on load!!!\n");
        return -1;
    } else {
        print_matrix(Asym, stdout, k, 20);
        //we need to make "non symmetric" in the sense of remove memory saving compression so that it can be used againts lower level blas
        A = unpack_symmetric_matrix(Asym);
        if (A == NULL) {
            delete_sssmatrix(S);
            delete_matrix(Asym);
            printf("***** !!!!  Symmetric matrix eigs FAILED on load!!!\n");
        } else {
            delete_matrix(Asym);
        }
    }

    //A and S available now
    print_matrix(A, stdout, k, 20);

    //initialize the vectors and matrices so that we can have safe freeing if unallocated
    Matrix* eigen_vectors = 0;
    Vector* eigen_values = 0;

    printf("\n\tLanczos dense matrix decomposition, please wait!!!....");
    EigenOrder order = LA;
    EigenSolverMode mode = I_REGULAR;
    //time(&start_time);
    int dense_return_value = dense_matrix_symmetric_eigen(A, k, order, mode, &eigen_values, &eigen_vectors);
    //time(&end_time);
    if (dense_return_value == 0) {
        printf("returned:: Success\n");
        //quick check, print first j eigen values
        printf("\n\tPrinting dense eigenvalues::");
        int ii;
        for (ii = 0; ii < k; ii++) {
            printf("\t%0.4f  ", eigen_values->data[ii]);
        }
        printf("\n\t");
        printf("\nPrinting eigenvectors::");
        print_matrix(eigen_vectors, stdout, k, k);
    } else {
        printf("returned:: PROBLEMS!! @ %d\n", dense_return_value);
    }

    delete_matrix(eigen_vectors);
    delete_vector(eigen_values);

    eigen_vectors = 0;
    eigen_values = 0;

    printf("\n\tLanczos sparse matrix decomposition, please wait!!!....");
    int sparse_return_value = sss_sparse_matrix_symmetric_eigen(S, k, order, mode, 1, &eigen_values, &eigen_vectors);
    if (sparse_return_value == 0) {
        printf("returned:: Success\n");
        //quick check, print first j eigen values
        printf("\n\tPrinting dense eigenvalues::");
        int ii;
        for (ii = 0; ii < k; ii++) {
            printf("\t%0.4f  ", eigen_values->data[ii]);
        }
        printf("\n\t");
        printf("\nPrinting eigenvectors::");
        print_matrix(eigen_vectors, stdout, k, k);
        //        matrix_print_top_left(V, 10);
    } else {
        printf("returned:: PROBLEMS!! @ %d\n", sparse_return_value);
    }

    delete_matrix(A);
    delete_matrix(eigen_vectors);
    delete_vector(eigen_values);
    delete_sssmatrix(S);

    return 0;
}

int arpack_svd_test() {
    sysmath_init(); //required for GNU/Linux with MKL backend to force use of OPEN MP and avoid shared object RTL conflicts with Intel OPENMP
    printf("\n\n======== Symmetric matrix svds :: ========\n");
    //
    int k = 10; // the rank we want
    char* sym_filename = BASE_DIR "Asym_mat_100.mtx";
    
    Matrix* Sym = 0;
    
    Matrix* Asym = load_mm_matrix(sym_filename);


    if (Asym == NULL) {
        printf("***** !!!!  Svds FAILED on load!!!\n");
        return -1;
    } else {
        //print_matrix(Asym, stdout, k, 20);
        //we need to make "non symmetric" in the sense of remove memory saving compression so that it can be used against lower level blas
        Sym = unpack_symmetric_matrix(Asym);
        if (Sym == NULL) {
            delete_matrix(Asym);
            printf("***** !!!!  Svds FAILED on load!!!\n");
            return -1;
        } else {
            delete_matrix(Asym);
        }
    }

    Matrix *U, *V;
    Vector *sigma;
    //initialize the vectors and matrices so that we can have safe freeing if unallocated
    U = 0;
    V = 0;
    sigma = 0;

    int sym_return_value = svd(Sym, k, true, true, &sigma, &U, &V);
    if (sym_return_value == 0) {
        printf("\nTruncated SVD with SYMMETRIC inputs returned:: Success\n");
        //quick check, print first j eigen values
        printf("\nPrinting symmetric singular values::");
        int ii;
        for (ii = 0; ii < k; ii++) {
            printf("\t%0.4f  ", sigma->data[ii]);
        }
        printf("\n\t");
        printf("\nPrinting right singular vectors::");
        print_matrix(U, stdout, k, k);

        printf("\nPrinting left singular vectors::");
        print_matrix(V, stdout, k, k);
        //
        delete_matrix(U);
        delete_matrix(V);
        delete_vector(sigma);
        //
        U = 0;
        V = 0;
        sigma = 0;
    } else {
        printf("\nTruncated SVD with SYMMETRIC inputs returned:: PROBLEMS!! @ %d\n", sym_return_value);
        //clean up
        delete_matrix(Sym);
        return -1;
    }
    
    printf("\n\n======== Non-symmetric matrix svds :: ========\n");   
    
    char* non_sym_filename = BASE_DIR "A_dns_dense_non_sym_mat_100.mtx";

    Matrix* NonSym = load_mm_matrix(non_sym_filename);

    if (NonSym == NULL) {
        printf("***** !!!!  Svds FAILED on load!!!\n");
        return -1;
    }
    int nonsym_return_value = svd(NonSym, k, true, true, &sigma, &U, &V);
    if (nonsym_return_value == 0) {
        printf("\nTruncated SVD with NON-SYMMETRIC inputs returned:: Success\n");
        //quick check, print first j eigen values
        printf("\n\tPrinting singular values::");
        int ii;
        for (ii = 0; ii < k; ii++) {
            printf("\t%0.4f  ", sigma->data[ii]);
        }
        printf("\n\t");
        printf("\nPrinting right singular vectors::");
        print_matrix(U, stdout, k, k);

        printf("\nPrinting left singular vectors::");
        print_matrix(V, stdout, k, k);
        //
        delete_matrix(U);
        delete_matrix(V);
        delete_vector(sigma);
    } else {
        printf("\nTruncated SVD with SYMMETRIC inputs returned:: PROBLEMS!! @ %d\n", nonsym_return_value);
        //clean up
        delete_matrix(Sym);
        delete_matrix(NonSym);
        return -1;
    }

    delete_matrix(Sym);
    delete_matrix(NonSym);
    return 0;
}