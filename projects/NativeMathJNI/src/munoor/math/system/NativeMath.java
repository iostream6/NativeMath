/*
 * 2017.05.13  Created 
 * 2017.06.03  Introduced setThreads method stub
 * 2017.06.15  Introduced allocateAlligned and freeAlligned
 * 2017.06.25  - Added sparse matrix client side support for Eigen routines
 * 2018.04.29  - SSS dedicated routine. Bulk allocate and release rotuines
 */
package munoor.math.system;

/**
 *
 * @author Ilamah, Osho
 */
public class NativeMath {

    public static native int sssSparseSymmetricEigs(final long values_pointer, final long dvalues_pointer, final long row_pointer, final long col_index_pointer, int rows, int columns, int nnzl, int rank, int order, int threads, long register);

    public static native int symmetricEigs(final long values_pointer, int rows, int columns, int rank, int order, long register);

    public static native void releaseRegisters(long dataPointers, int registers);

    public static native void setKernelThreads(int numberOfThreads);

    public static native long allocateAllignedBuffer(int elements, int bytesPerElement, int byteAllignment);

    //return contain the origin of each buffer that was successfully allocated or zero
    public static native long[] allocateAllignedBuffers(int[] elementsPerBuffer, int[] bytesPerElement, int byteAllignment);

    public static native void releaseAllignedBuffer(long pointer);
    
    public static native void releaseAllignedBuffers(long[] pointers);
}
