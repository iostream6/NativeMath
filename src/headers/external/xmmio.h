/*
 * 2018.04.17  - Created
 */

/* 
 * Provides function prototypes for matrix input/output routines 
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

    /**
     * Reads data from matrix market formatted text file into Matrix data structure. Both dense and sparse matrix market files can be read.
     * 
     * @param filename the absolute path to the input file
     * @return a Matrix structure containing the input data if successful, otherwise NULL
     */
    Matrix* load_mm_matrix(const char* filename);

    /**
     * Reads data from matrix market formatted text file into CSRMatrix data structure. Only sparse, non-symmetric matrix market input data is supported. If sparse symmetric input is required, please use the Matrix and SSSMatrix structure and associated 
     * read methods instead.
     * 
     * @param filename the absolute path to the input file
     * @return  a CSRMatrix structure containing the input data if successful, otherwise NULL
     */
    CSRMatrix* load_mm_csrmatrix(char* filename);
    
    /**
     * Reads data from matrix market formatted text file into SSSMatrix data structure. Only sparse, symmetric matrix market input data is supported. If sparse NON-symmetric input is required, please use the Matrix and CSRMatrix structure and associated 
     * read methods instead.
     * 
     * @param filename  the absolute path to the input file
     * @return  a SSSMatrix structure containing the input data if successful, otherwise NULL
     */
    SSSMatrix* load_mm_ssmatrix(char* filename);

#ifdef __cplusplus
}
#endif

#endif /* XMMIO_H */

