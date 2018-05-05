/*
 * 2017.05.13  - Created 
 * 2017.06.03  - Introduced setThreads method stub
 * 2017.06.15  - Introduced allocateAlligned and freeAlligned
 * 2017.06.25  - Added sparse matrix client side support for Eigen routines
 * 2018.04.29  - SSS dedicated routine. Bulk allocate and release rotuines
 */
package munoor.math.system;

/**
 *
 * @author Ilamah, Osho
 */
public class NativeMath {

    public static final int LA_EIGEN_ORDER = 1;
    public static final int SA_EIGEN_ORDER = 2;
    public static final int LM_EIGEN_ORDER = 3;
    public static final int SM_EIGEN_ORDER = 4;
    public static final int BE_EIGEN_ORDER = 99;

    /**
     * Native interface for truncated eigen decomposition for sparse symmteric matrix (organised in SSS format). The low level implementation uses ARPACK and a custom spMVx kernel.
     *
     * @param values_pointer pointer to double precision floating point buffer containing lower triangle non-zero elements of the input matrix
     * @param dvalues_pointer pointer to double precision floating point buffer containing the diagonal elements of the input matrix
     * @param row_pointer pointer to integer buffer containing input matrix row start indices (a.k.a row pointers)
     * @param col_index_pointer pointer to integer buffer containing input matrix column indices for the lower triangle non-zero elements
     * @param rows number of rows in the input matrix
     * @param columns number of columns in the input matrix
     * @param nnzl number of lower triangle non-zero element of the input matrix
     * @param rank the number of eigen pairs required
     * @param order specifies the eigen ordering
     * @param threads specifies the number of parallel threads to be used in the spMVx kernel
     * @param register 2 element "long" register to which memory addresses which point to the eigen results will be stored - 1st element is eigen value address, 2nd is eigen vector address
     * @return zero if the native call completed without errors, otherwise a non zero value is retured
     */
    public static native int sssSparseSymmetricEigs(final long values_pointer, final long dvalues_pointer, final long row_pointer, final long col_index_pointer, int rows, int columns, int nnzl, int rank, int order, int threads, long register);

    /**
     * Native interface for truncated eigen decomposition for dense symmteric matrix. The low level implementation uses ARPACK and OPENBLAS (or MKL) MVx routine.
     *
     * @param values_pointer pointer to double precision floating point buffer containing the elements of the input matrix
     * @param rows number of rows in the input matrix
     * @param columns number of columns in the input matrix
     * @param rank the number of eigen pairs required
     * @param order specifies the eigen ordering
     * @param register 2 element "long" register to which memory addresses which point to the eigen results will be stored - 1st element is eigen value address, 2nd is eigen vector address
     * @return
     */
    public static native int symmetricEigs(final long values_pointer, int rows, int columns, int rank, int order, long register);

    /**
     * Releases one or more heap memory buffers. The buffers to be released are assumed to have been previously allocated by one of the alligned allocation methods
     *
     * @param dataPointers n element "long" register whose elements are memory addresses which point to the alligned buffers to be released
     * @param n the number of elements in the valid elements in the dataPointers registers
     */
    public static native void releaseRegisters(long dataPointers, int n);

    /**
     * Sets the concurency level for the core native math kernel routines and associated OPENBLAS or MKL backends
     *
     * @param numberOfThreads the desired concurency level for the core native math kernel routines
     */
    public static native void setKernelThreads(int numberOfThreads);

    /**
     * Allocates a byte alligned memory buffer and clears its contents.
     *
     * @param elements the number of elements expected to be stored in the buffer
     * @param bytesPerElement the bytes required for each element of the buffer
     * @param byteAllignment the byte allignment. if zero or negative, 64 bytes alligments is used
     * @return if sucessful, the address of the byte alligned, allocated and cleared memory area is returned, else a zero is returned
     */
    public static native long allocateAllignedBuffer(int elements, int bytesPerElement, int byteAllignment);

    /**
     * Allocates 1 or more byte alligned memory buffers and clears their contents.
     *
     * @param elementsPerBuffer an array of integers defining the number of elements expected to be stored in each buffer
     * @param bytesPerElement an array of integers defining the bytes required for each element of each of the buffers
     * @param byteAllignment the byte allignment. if zero or negative, 64 bytes alligments is used
     * @return an array of longs where each element is, if sucessful, the address of the byte alligned, allocated and cleared memory buffer, otherwise is a zero value, is returned
     */
    public static native long[] allocateAllignedBuffers(int[] elementsPerBuffer, int[] bytesPerElement, int byteAllignment);

    /**
     * Releases a heap memory buffer. The buffer to be released is assumed to have been previously allocated by one of the alligned allocation methods
     *
     * @param pointer long value which which points to the alligned buffer to be released
     */
    public static native void releaseAllignedBuffer(long pointer);

    /**
     * Releases one or more heap memory buffers. The buffers to be released are assumed to have been previously allocated by one of the alligned allocation methods. This method should be used in
     * preference to the "releaseRegisters" method, although internally, they perform the same function
     *
     * @param pointers "long" register whose elements are memory addresses which point to the alligned buffers to be released
     */
    public static native void releaseAllignedBuffers(long[] pointers);
}
