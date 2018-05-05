/*
 * 2017.05.17  - Created
 * 2017.05.24  - QC passed
 * 2017.06.03  - Introduced Java_munoor_math_system_MatrixUtils_setThreads method
 * 2017.06.15  - Introduced Java_munoor_math_system_MatrixUtils_allocateAlligned and Java_munoor_math_system_MatrixUtils_freeAlligned
 * 2017.06.25  - Added basic sparse matrix support for Eigen routines
 * 2018.04.30  - Restructure and introduce SSS support. Added bulk allocate. Rename methods/Java class. Overall QC required
 */
/* 
 * File:   jni.c
 * Author: Ilamah, Osho
 *
 * Created on 17 May 2017, 21:24
 */

#include <math.h>
#include <stdbool.h>

#include "../headers/external/arpack.h"
#include "../headers/external/xalloc.h"
#include "../headers/external/matvec.h"
#include "../headers/internal/munoor_math_system_NativeMath.h"

//#include "../interfaces/public/utils.h"
//#include "../interfaces/public/rmd.h"

/*
 *  Modified, removing JNIEXPORT to use instead export lists. As otherwise non JNI functions are hidden
 *  http://anadoxin.org/blog/control-over-symbol-exports-in-gcc.html,   https://www.ibm.com/developerworks/aix/library/au-aix-symbol-visibility/,  https://github.com/stack-of-tasks/jrl-doc/wiki/Symbol-visibility-management
 *  http://man7.org/conf/lca2006/shared_libraries/slide18b.html,  http://stackoverflow.com/questions/19422660/when-to-use-jniexport-and-jnicall-in-android-ndk,  https://gcc.gnu.org/wiki/Visibility
 */


/*
 * Class:     munoor_math_system_NativeMath
 * Method:    sssSparseSymmetricEigs
 * Signature: (JJJJIIIIIJ)I
 */
//QC'ed
/*JNIEXPORT*/ jint JNICALL Java_munoor_math_system_NativeMath_sssSparseSymmetricEigs(JNIEnv* env, jclass class, jlong values_pointer, jlong dvalues_pointer, jlong row_pointer, jlong col_index_pointer, jint m, jint n, jint nnzl,
        jint k, jint o, jint threads, jlong outputPtr) {
    sysmath_init();
    //Sanity checks expected to have been done in JVM layer. The below is based on codes from test_eigen_arpack() function in SysMathTest
    SSSMatrix* S = wrap_as_sssmatrix((uint32_t) m, (uint32_t) n, (uint32_t) nnzl, ((double*) dvalues_pointer), ((double*) values_pointer), ((uint32_t*) row_pointer), ((uint32_t*) col_index_pointer));
    const uint32_t NUM_THREADS = threads > 0 ? (uint32_t) threads : 1;
    //create and initialize the vectors and matrices so that we can have safe freeing if unallocated
    Matrix *eigenVectors = 0;
    Vector *eigenValues = 0;
    //
    EigenOrder order;
    switch (o) {
        case 1:
            order = LA;
            break;
        case 2: //SA
            order = SA;
            break;
        case 3: //LM
            order = LM;
            break;
        case 4: //SM
            order = SM;
            break;
        default: //BE
            order = BE;
            break;
    }
    EigenSolverMode mode = I_REGULAR;
    //time_t start_time, end_time;
    //time(&start_time);
    int returnValue = sss_sparse_matrix_symmetric_eigen(S, k, order, mode, NUM_THREADS, &eigenValues, &eigenVectors);
    //time(&end_time);
    free(S); // Just the outer struct, the buffers data memory is still intact. 
    //
    // Return status and optionally results
    // outputPtr comes from jvm as a 2*8BYTE memory location
    // We set from C side as "jlong" (or long long or __int64) to ensure it takes up the 64 bits per entry (on windows long is 32 byte due to LLP, see http://stackoverflow.com/a/7607537 and https://docs.oracle.com/javase/8/docs/technotes/guides/jni/spec/types.html
    if (returnValue == 0) {
#ifdef  __DEBUG__
        printf("\n======================================================================********");
        printf("\n======================================================================********");
        printf("\n======================================================================********");
        printf("\nNATIVE:: Printing Eigenvalues::");
        int ii;
        for (ii = 0; ii < eigenValues->nrows; ii++) {
            printf("\t%f  ", eigenValues->data[ii]);
        }
        printf("\n");
        printf("\nNATIVE:: Printing Eigenvectors::");
        matrix_print_top_left(eigenVectors, k);
        printf("\n======================================================================********");
        printf("\n======================================================================********");
        printf("\n======================================================================********");
#endif
        //
        jlong * outputRegister = (jlong *) outputPtr;
        outputRegister[0] = (jlong) eigenValues->data;
        outputRegister[1] = (jlong) eigenVectors->data;
        //
        free(eigenValues); // Just the outer struct, the data memory is still intact.
        free(eigenVectors);
    } else {
        delete_vector(eigenValues);
        delete_matrix(eigenVectors);
    }
    return returnValue;
}

/*
 * Class:     munoor_math_system_NativeMath
 * Method:    symmetricEigs
 * Signature: (JIIIIJ)I
 */
//QC'ed
/*JNIEXPORT*/ jint JNICALL Java_munoor_math_system_NativeMath_symmetricEigs(JNIEnv* env, jclass class, jlong values_pointer, jint m, jint n, jint k, jint o, jlong outputPtr) {
    sysmath_init();
    //Sanity checks expected to have been done in JVM layer. The below is based on codes from test_eigen_arpack() function in SysMathTest
    Matrix *A = wrap_as_matrix(m, n, false /*even though symmetric, blas does not support compressed format)*/, ((double*) values_pointer));

    //create and initialize the vectors and matrices so that we can have safe freeing if unallocated
    Matrix *eigenVectors = 0;
    Vector *eigenValues = 0;
    //
    EigenOrder order;
    switch (o) {
        case 1:
            order = LA;
            break;
        case 2: //SA
            order = SA;
            break;
        case 3: //LM
            order = LM;
            break;
        case 4: //SM
            order = SM;
            break;
        default: //BE
            order = BE;
            break;
    }
    EigenSolverMode mode = I_REGULAR;
    //time_t start_time, end_time;
    //time(&start_time);
    int returnValue = dense_matrix_symmetric_eigen(A, k, order, mode, &eigenValues, &eigenVectors);
    //time(&end_time);
    free(A); // Just the outer struct, the data memory is still intact. 
    //
    // Return status and optionally results
    // outputPtr comes from jvm as a 2*8BYTE memory location
    // We set from C side as "jlong" (or long long or __int64) to ensure it takes up the 64 bits per entry (on windows long is 32 byte due to LLP, see http://stackoverflow.com/a/7607537 and https://docs.oracle.com/javase/8/docs/technotes/guides/jni/spec/types.html
    if (returnValue == 0) {
#ifdef  __DEBUG__
        printf("\n======================================================================********");
        printf("\n======================================================================********");
        printf("\n======================================================================********");
        printf("\nNATIVE:: Printing Eigenvalues::");
        int ii;
        for (ii = 0; ii < eigenValues->nrows; ii++) {
            printf("\t%f  ", eigenValues->data[ii]);
        }
        printf("\n");
        printf("\nNATIVE:: Printing Eigenvectors::");
        matrix_print_top_left(eigenVectors, k);
        printf("\n======================================================================********");
        printf("\n======================================================================********");
        printf("\n======================================================================********");
#endif
        jlong * outputRegister = (jlong *) outputPtr;
        outputRegister[0] = (jlong) eigenValues->data;
        outputRegister[1] = (jlong) eigenVectors->data;
        //
        free(eigenValues); // Just the outer struct, the data memory is still intact.
        free(eigenVectors);
    } else {
        delete_vector(eigenValues);
        delete_matrix(eigenVectors);
    }
    return returnValue;
}

/*
 * Class:     munoor_math_system_NativeMath
 * Method:    releaseRegisters
 * Signature: (JI)V
 */
//QC'ed
/*JNIEXPORT*/ void JNICALL Java_munoor_math_system_NativeMath_releaseRegisters(JNIEnv* env, jclass class, jlong registersPtr, jint numberofRegisters) {
    jlong * memRegisterPtr = (jlong*) registersPtr; //need this to be java longs to deal with LLP64 and LP64
    for (int index = 0; index < numberofRegisters; index++) {
        aligned_free((void *) memRegisterPtr[index]);
    }
}

/*
 * Class:     munoor_math_system_NativeMath
 * Method:    setKernelThreads
 * Signature: (I)V
 */
//QC'ed
/*JNIEXPORT*/ void JNICALL Java_munoor_math_system_NativeMath_setKernelThreads(JNIEnv* env, jclass class, jint numThreads) {
    set_system_num_threads(numThreads);
}

/*
 * Class:     munoor_math_system_NativeMath
 * Method:    allocateAllignedBuffer
 * Signature: (III)J
 */
//QC'ed
/*JNIEXPORT*/ jlong JNICALL Java_munoor_math_system_NativeMath_allocateAllignedBuffer(JNIEnv* env, jclass class, jint elements, jint bytes_per_element, jint byte_allignment) {
    return (jlong) (aligned_calloc(elements, bytes_per_element, byte_allignment));
}

/*
 * Class:     munoor_math_system_NativeMath
 * Method:    allocateAllignedBuffers
 * Signature: ([I[II)[J
 */
//QC'ed
/*JNIEXPORT*/ jlongArray JNICALL Java_munoor_math_system_NativeMath_allocateAllignedBuffers(JNIEnv* env, jclass class, jintArray elements_per_buffer, jintArray bytes_per_element, jint byte_allignment) {
    const jsize len = (*env)->GetArrayLength(env, elements_per_buffer);
    //allocate C buffer for the inputs (all on stack)
    jint c_elements_per_buffer[len];
    jint c_bytes_per_element[len];
    //allocate C buffer for the output
    jlong c_created_buffer_address[len];
    //
    //copy the JVM array into C - this makes sense as the array is small (see JNI book)
    (*env)->GetIntArrayRegion(env, elements_per_buffer, 0, len, c_elements_per_buffer);
    (*env)->GetIntArrayRegion(env, bytes_per_element, 0, len, c_bytes_per_element);
    //do the allocation
    for(int i = 0; i < len; i++){
        jlong address = (jlong) (aligned_calloc(c_elements_per_buffer[i], c_bytes_per_element[i], byte_allignment));
        c_created_buffer_address[i] = address;
    }
    //create a result array in the JVM
    jlongArray results=(jlongArray)(*env)->NewLongArray(env,len);
    //copy the c results buffer back into JVM
    (*env)->SetLongArrayRegion(env, results, 0, len, c_created_buffer_address);
    return results;
}

/*
 * Class:     munoor_math_system_NativeMath
 * Method:    releaseAllignedBuffer
 * Signature: (J)V
 */
//QC'ed
/*JNIEXPORT*/ void JNICALL Java_munoor_math_system_NativeMath_releaseAllignedBuffer(JNIEnv* env, jclass class, jlong pointer) {
    void* memRegisterPtr = (void*) pointer; //need this to be java longs to deal with LLP64 and LP64
    aligned_free(memRegisterPtr);
}

/*
 * Class:     munoor_math_system_NativeMath
 * Method:    releaseAllignedBuffers
 * Signature: ([J)V
 */
//QC'ed
/*JNIEXPORT*/ void JNICALL Java_munoor_math_system_NativeMath_releaseAllignedBuffers(JNIEnv* env, jclass class, jlongArray pointers){
    const jsize len = (*env)->GetArrayLength(env, pointers);
    //allocate C buffer for the input (all on stack)
    jlong c_pointers[len];
    //
    //copy the JVM array into C - this makes sense as the array is small (see JNI book)
    (*env)->GetLongArrayRegion(env, pointers, 0, len, c_pointers);
    //
    //perform the release operation
    for(int i = 0; i < len; i++){
        void* memRegisterPtr = (void*) c_pointers[i]; //need this to be java longs to deal with LLP64 and LP64
        aligned_free(memRegisterPtr);
    }
}

///*
// * Class:     system_math_NativeMatrixUtils
// * Method:    getTruncatedSVD
// * Signature: (JIIIIZZJ)I
// */
//
///*JNIEXPORT*/ jint JNICALL Java_munoor_math_system_MatrixUtils_truncatedSVD
//(JNIEnv *env, jclass theClass, jlong matrixDataPtr, jint m, jint n, jint k, jboolean solveLeftVectors, jboolean solveRightVectors, jlong outputPtr) {
//    sysmath_init();
//    //Sanity checks expected to have been done in JVM layer. The below is based on codes from test_eigen_arpack() function in SysMathTest
//    Matrix *A = matrix_wrap((uint32_t)m, (uint32_t)n, ((double*) matrixDataPtr));
//    //matrix_print_top_left(A, 10);
//    Matrix *U, *V;
//    vec *sigma;
//    //initialize the vectors and matrices so that we can have safe freeing if unallocated
//    U = 0;
//    V = 0;
//    sigma = 0;
//    //
//    //time_t start_time, end_time;
//    //time(&start_time);
//    int returnValue = svd(A, k, solveLeftVectors, solveRightVectors, &sigma, &U, &V);
//    //time(&end_time);
//    free(A); // Just the outer struct, the data memory is still intact. 
//    //
//    // Return status and optionally results
//    // outputPtr comes from jvm as a 3*8BYTE memory location
//    // We set from C side as "jlong" (or long long or __int64) to ensure it takes up the 64 bits per entry (on windows long is 32 byte due to LLP, see http://stackoverflow.com/a/7607537 and https://docs.oracle.com/javase/8/docs/technotes/guides/jni/spec/types.html
//    if (returnValue == 0) {
//        jlong * outputRegister = (jlong *) outputPtr;
//        outputRegister[0] = (jlong) sigma->data;
//        free(sigma); // Just the outer struct, the data memory is still intact.
//        if (solveLeftVectors) {
//            outputRegister[1] = (jlong) U->data;
//            free(U); // Just the outer struct, the data memory is still intact.
//        } else {
//            matrix_delete(U); //placeholder and the data
//        }
//        if (solveRightVectors) {
//            outputRegister[2] = (jlong) V->data;
//            free(V); // Just the outer struct, the data memory is still intact.
//        } else {
//            matrix_delete(V); //placeholder and the data
//        }
//    } else {
//        vector_delete(sigma);
//        matrix_delete(U);
//        matrix_delete(V);
//    }
//    return returnValue;
//}
//
//



