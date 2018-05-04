///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//
//#include <stdlib.h>
//#include <stdio.h>
//#include <time.h>
//
//
//#include "../../src/headers/external/matvec.h"
//#include "../../src/headers/external/xblas.h"
//#include "../../src/headers/external/xmmio.h"
//#include "../../src/headers/external/xalloc.h"
//
//#define WINDOW_SIZE 10
//#define SMALL_MAT_SIZE 20
//
//int allocate_and_delete_memory();
//int create_and_delete_matrices();
//int load_matrices_from_file();
//int csr_from_matrix_file();
//int sss_from_matrix_file();
//int load_matrices_from_file_and_convert();
//
//int spMxV();
//int load_matrices_from_file_and_unpack();
//
//int main() {
//    time_t t;
//    /* Intializes random number generator */
//    srand((unsigned) time(&t));
//
//    //        if (allocate_and_delete_memory() == 0) {
//    //            printf("Allocate and Delete Memory Okay");
//    //        } else {
//    //            printf("Create and Delete Memory FAILED!!!!!");
//    //            return EXIT_FAILURE;
//    //        }
//    //        printf("\n=================================\n");
//    //
//    //        if (create_and_delete_matrices() == 0) {
//    //            printf("Matrix create and delete method Okay");
//    //        } else {
//    //            printf("Matrix create and delete method FAILED!!!!!");
//    //            return EXIT_FAILURE;
//    //        }
//    //
//    //        printf("\n=================================\n");
//    //        if (load_matrices_from_file() == 0) {
//    //            printf("Matrix market matrix load:: Okay");
//    //        } else {
//    //            printf("Matrix market matrix load:: FAILED!!!!!");
//    //            return EXIT_FAILURE;
//    //        }
//    //
//    //
//    //        printf("\n=================================\n");
//    //        if (csr_from_matrix_file() == 0) {
//    //            printf("CSRMatrix from MATRIX:: Okay");
//    //        } else {
//    //            printf("CSRMatrix from MATRIX:: FAILED!!!!!");
//    //            return EXIT_FAILURE;
//    //        }
//    //
//            printf("\n=================================\n");
//            if (sss_from_matrix_file() == 0) {
//                printf("SSSMatrix from MATRIX:: Okay");
//            } else {
//                printf("SSSMatrix from MATRIX:: FAILED!!!!!");
//                return EXIT_FAILURE;
//            }
//
//    //                printf("\n=================================\n");
//    //                if (load_matrices_from_file_and_convert() == 0) {
//    //                    printf("Matrix conversion:: Okay");
//    //                } else {
//    //                    printf("Matrix conversion:: FAILED!!!!!");
//    //                    return EXIT_FAILURE;
//    //                }
////    printf("\n=================================\n");
////    if (spMxV() == 0) {
////        printf("Sparse Matrix-Vector product:: Okay");
////    } else {
////        printf("Sparse Matrix-Vector product:: FAILED!!!!!");
////        return EXIT_FAILURE;
////    }
//
//
//
//
//    //load_matrices_from_file_and_unpack();
//
//    printf("\n\nPlease type a character to exit::");
//    //getchar();
//    return EXIT_SUCCESS;
//
//}
//
//int allocate_and_delete_memory() {
//    const uint32_t elements = 200;
//    const size_t size = sizeof (elements);
//    double* buffer = aligned_calloc(elements, size, 64);
//    if (buffer == NULL) {
//        return -1;
//    }
//    aligned_free(buffer);
//    return 0;
//}
//
//int create_and_delete_matrices() {
//    Matrix* sym = create_matrix(SMALL_MAT_SIZE, SMALL_MAT_SIZE, true);
//    Matrix* nonsym = create_matrix(SMALL_MAT_SIZE, SMALL_MAT_SIZE, false);
//    for (int i = 0; i < SMALL_MAT_SIZE; i++) {
//        for (int j = i; j < SMALL_MAT_SIZE; j++) {
//            const int val = (rand() % 85) + 10;
//            set_matrix_element(sym, i, j, val);
//            //
//            set_matrix_element(nonsym, i, j, val);
//            set_matrix_element(nonsym, j, i, ((rand() % 85) + 10));
//        }
//    }
//    //#ifdef __DEBUG__
//    //        printf("=====================SYM_MATRIX\n\n");
//    //        print_matrix(sym, stdout, SMALL_MAT_SIZE, SMALL_MAT_SIZE);
//    //        printf("=====================NONSYM_MATRIX\n\n");
//    //        print_matrix(nonsym, stdout, SMALL_MAT_SIZE, SMALL_MAT_SIZE);
//    //#endif
//    // delete the matrices
//    delete_matrix(sym);
//    delete_matrix(nonsym);
//    return 0;
//}
//
//int load_matrices_from_file() {
//    //
//    //char* filename = "X:\\delete\\x\\Ragusa16\\test.mtx";
//    char* filename = "X:\\delete\\x\\Ragusa16\\gk.mtx";
//    //char* filename = "X:\\delete\\x\\Ragusa16\\Ragusa16.mtx";
//    //char* filename = "/media/MERCURY/Samba/DEL/Ragusa16.mtx";
//    //    
//    Matrix* m = load_mm_matrix(filename);
//    if (m == NULL) {
//        return -1;
//    }
//    printf("File read");
//
//#ifdef __DEBUG__
//    printf("=====================MARKET MATRIX\n\n");
//    print_matrix(m, stdout, SMALL_MAT_SIZE, SMALL_MAT_SIZE);
//    printf("=====================\n\n");
//#endif
//    delete_matrix(m);
//    return 0;
//}
//
//int csr_from_matrix_file() {
//    char* filename = "X:\\delete\\x\\Ragusa16\\test.mtx";
//    //char* filename = "X:\\delete\\x\\Ragusa16\\Ragusa16.mtx";
//    //char* filename = "/media/MERCURY/Samba/DEL/Ragusa16.mtx";
//    //   
//    //    
//    CSRMatrix* m = load_mm_csrmatrix(filename);
//    if (m == NULL) {
//        return -1;
//    }
//
//
//#ifdef __DEBUG__
//    printf("=====================CSR MARKET MATRIX\n\n");
//    print_csr_matrix(m, stdout, 10);
//    printf("=====================\n\n");
//#endif
//    delete_csrmatrix(m);
//    return 0;
//}
//
//int sss_from_matrix_file() {
//    char* filename = "X:\\delete\\x\\Ragusa16\\LFAT5.mtx";
//    //   
//    //    
//    SSSMatrix* m = load_mm_ssmatrix(filename);
//    if (m == NULL) {
//        return -1;
//    }
//
//#ifdef __DEBUG__
//    printf("=====================SSS MARKET MATRIX\n\n");
//    print_sss_matrix(m, stdout, 30);
//    printf("=====================\n\n");
//#endif
//    delete_sssmatrix(m);
//    return 0;
//}
//
//int load_matrices_from_file_and_convert() {
//    //char* sym_filename = "/media/MERCURY/Samba/DEL/nos2.mtx";
//    //char* sym_filename = "X:\\delete\\x\\Ragusa16\\gk.mtx";
//    char* sym_filename = "X:\\delete\\x\\Ragusa16\\nos2.mtx";
//    char* filename = "X:\\delete\\x\\Ragusa16\\test.mtx";
//
//    Matrix* matrix_file = load_mm_matrix(filename);
//    CSRMatrix* csr_matrix_file = load_mm_csrmatrix(filename);
//
//    Matrix* sym_matrix_file = load_mm_matrix(sym_filename);
//    SSSMatrix* sss_matrix_file = load_mm_ssmatrix(sym_filename);
//    if (matrix_file == NULL || csr_matrix_file == NULL /*|| sss_matrix_file == NULL*/) {
//        return -1;
//    }
//
//    CSRMatrix* csr_matrix_conv = csr_from_matrix(matrix_file);
//    SSSMatrix* sss_matrix_conv = sss_from_matrix(sym_matrix_file);
//
//    const double epsilon = 1e-8;
//
//    for (uint32_t ii = 0; ii < csr_matrix_file->nnz; ii++) {
//        if ((csr_matrix_file->values[ii] - csr_matrix_conv->values[ii]) > epsilon) {
//            printf("CSR FILE vs CONV DELTA on values[%d] is %lg \n", ii, (csr_matrix_file->values[ii] - csr_matrix_conv->values[ii]));
//            return -1;
//        }
//    }
//
//    for (uint32_t ii = 0; ii < sss_matrix_file->nnzl; ii++) {
//        if ((sss_matrix_file->values[ii] - sss_matrix_conv->values[ii]) > epsilon) {
//            printf("SSS FILE vs CONV DELTA on values[%d] is %lg \n", ii, (sss_matrix_file->values[ii] - sss_matrix_conv->values[ii]));
//            return -1;
//        }
//        if ((sss_matrix_file->col_ind[ii] - sss_matrix_conv->col_ind[ii]) > epsilon) {
//            printf("SSS FILE vs CONV DELTA on col_ind[%d] is %d \n", ii, (sss_matrix_file->col_ind[ii] - sss_matrix_conv->col_ind[ii]));
//            return -1;
//        }
//    }
//
//    for (uint32_t ii = 0; ii < sss_matrix_file->nrows; ii++) {
//        if ((sss_matrix_file->dvalues[ii] - sss_matrix_conv->dvalues[ii]) > epsilon) {
//            printf("SSS FILE vs CONV DELTA on dvalues[%d] is %lg \n", ii, (sss_matrix_file->dvalues[ii] - sss_matrix_conv->dvalues[ii]));
//            return -1;
//        }
//        if ((sss_matrix_file->row_ptr[ii] - sss_matrix_conv->row_ptr[ii]) > epsilon) {
//            printf("SSS FILE vs CONV DELTA on row_ptr[%d] is %d \n", ii, (sss_matrix_file->row_ptr[ii] - sss_matrix_conv->row_ptr[ii]));
//            return -1;
//        }
//    }
//
//#ifdef __DEBUG__
//    //    printf("=====================CSR MATRIX == CSR MATRIX CONV \n\n");
//    //    print_csr_matrix(csr_matrix_file, stdout, 8);
//    //    printf("=====================\n\n");
//    //
//    //    printf("=====================SSS MATRIX  == SSS MATRIX CONV\n\n");
//    //    print_sss_matrix(sss_matrix_file, stdout, 30);
//    //    printf("=====================\n\n");
//#endif
//
//    return 0;
//}
//
//int load_matrices_from_file_and_unpack() {
//    char* sym_filename = "X:\\delete\\x\\Ragusa16\\nos2.mtx";
//    Matrix* packed_sym_matrix = load_mm_matrix(sym_filename);
//    if (packed_sym_matrix == NULL) {
//        return -1;
//    }
//
//    Matrix* unpacked_sym_matrix = unpack_symmetric_matrix(packed_sym_matrix);
//    if (unpacked_sym_matrix == NULL) {
//        return -1;
//    }
//
//#ifdef __DEBUG__
//    for (uint32_t i = 0; i < packed_sym_matrix->nrows; i++) {
//        for (uint32_t j = 0; j < packed_sym_matrix->ncols; j++) {
//            if (get_matrix_element(packed_sym_matrix, i, j) != get_matrix_element(unpacked_sym_matrix, i, j)) {
//                printf("\n\n-------------Matrix unpack FAILED!!!\n");
//                return -1;
//            }
//        }
//    }
//    printf("\n\n-------------Matrix unpack OK!!!\n");
//#endif
//
//    return 0;
//}
//
//int spMxV() {
//    //perform/test sparse matrix vector kernels
//    //char* sym_filename = "X:\\delete\\x\\Ragusa16\\nos2.mtx";
//    //char* sym_filename = "X:\\delete\\x\\Ragusa16\\c-73.mtx";
//    char* sym_filename = "X:\\delete\\x\\Ragusa16\\LFAT5.mtx";
//    //char* sym_filename = "/media/MERCURY/Samba/DEL/nos2.mtx";
//    //Matrix* sym_matrix_file = load_mm_matrix(sym_filename);
//    SSSMatrix* sss_matrix = load_mm_ssmatrix(sym_filename);
//    Vector* x = create_vector(sss_matrix->nrows);
//    for (uint32_t row = 0; row < x->nrows; row++) {
//        x->data[row] = row + 1/*row / (2.0)*/;
//        //printf("%lg ", x->data[row]);
//        //const int val = (rand() % 200) - 101;
//        //->data[row] = val;
//    }
//    Vector* y_serial = create_vector(sss_matrix->nrows);
//    Vector* y_parallel = create_vector(sss_matrix->nrows);
//
//    clock_t start = clock(), diff;
//    for (int iii = 0; iii < 1; iii++) {
//        sparse_sym_matrix_vector_mult(sss_matrix, x->data, y_serial->data, 1);
//        sparse_sym_matrix_vector_mult(sss_matrix, x->data, y_parallel->data, 4);
//    }
//    diff = clock() - start;
//    int msec = diff * 1000 / CLOCKS_PER_SEC;
//    printf("\n\nTime taken %d seconds %d milliseconds\n", msec / 1000, msec % 1000);
//    const double epsilon = 1e-8;
//    const int top = 40 > y_serial->nrows ? y_serial->nrows : 40;
//    printf("\n\n-------------Xerial spMxV Kernel Result:: Top %d ->\n\n", top);
//    for (uint32_t row = 0; row < x->nrows; row++) {
//        if ((y_serial->data[row] - y_parallel->data[row]) > epsilon) {
//            printf("DELTA on row %d is %lg \n", row, y_serial->data[row] - y_parallel->data[row]);
//            return -1;
//        }
//        if (row < top) {
//            printf("%lg \n", y_serial->data[row]);
//        }
//    }
//    printf("\n\n-------------Xerial spMxV Kernel All Good!!!\n");
//    delete_sssmatrix(sss_matrix);
//    delete_vector(x);
//    delete_vector(y_serial);
//    delete_vector(y_parallel);
//    printf("\n\n-------------Xerial spMxV Kernel DEALLOCATE All Good!!!\n");
//
//    //    //serial spMxV kernel::  GKountouvas et al. IEEE, 2013
//    //    for (uint32_t row = 0; row < x->nrows; row++) {
//    //        y->data[row] = sss_matrix_conv->dvalues[row] * x->data[row];
//    //        const uint32_t start = sss_matrix_conv->row_ptr[row], end = sss_matrix_conv->row_ptr[row + 1];
//    //        for (uint32_t j = start; j < end; j++) {
//    //            const uint32_t col = sss_matrix_conv->col_ind[j];
//    //            const double value = sss_matrix_conv->values[j];
//    //            y->data[row] += value * x->data[col];
//    //            y->data[col] += value * x->data[row];
//    //        }
//    //    }
//    //    const int top = 40 > y->nrows ? y->nrows : 40;
//    //    //print result
//    //    printf("\n\n-------------Serial spMxV Kernel result vector:: top %d\n\n", top);
//    //    for (uint32_t row = 0; row < top; row++) {
//    //        printf("%lg \n", y->data[row]);
//    //    }
//    //    spMxVNaiive(sss_matrix_conv, x, 7);
//    return 0;
//}