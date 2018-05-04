/*
 * 2018.04.17  - Created
 */

/* 
 * Provides function prototypes for matrix loading routines 
 * 
 * 
 * File:   xmmio.h
 * Author: Ilamah, Osho
 *
 * Created on 17 April 2018, 09:18
 */

#ifndef XMMIO_H
#define XMMIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include "matvec.h"
    /*
     * A "co-ordinate list" sparse matrix structure. This is an unoptimised sparse matrix primitive 
     * intended as an intermediate transfer structure (e.g. when reading Market Matrix files into 
     * other compact sparse matrix structures)
     */
    typedef struct {
        uint32_t nrows, ncols; //number of cols and rows
        uint32_t nnz, nnzr, nnzd, nnzl; 
        uint32_t* row_coord;
        uint32_t* col_coord; // arrays containing row and column co-ordinate
        double* values;
    } COOMatrix;


    //loads data from matrix market format file into matrix structure
    Matrix* load_mm_matrix(const char* filename);
    //loads data from matrix market format file into csr matrix structure
    CSRMatrix* load_mm_csrmatrix(char* filename);
    //loads data from matrix market format file into sss matrix structure
    SSSMatrix* load_mm_ssmatrix(char* filename);

    //
    /**
     * Loads data from a proprietary binary file into a matrix (see matlab scripts that do this; Format = int m, int n, row major double data). 
     * @param fname pointer to the file name
     * @return A matrix containing the data in the file
     */
    //Matrix * matrix_load_from_binary_file(char *fname);

#ifdef __cplusplus
}
#endif

#endif /* XMMIO_H */

