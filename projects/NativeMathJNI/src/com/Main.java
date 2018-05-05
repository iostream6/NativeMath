/*
 * 2018.05.01  - Created
 * 2018.05.05  - Added test methods for symmetric eigs JNI calls
 */
package com;

import munoor.math.system.COOMatrix;
import munoor.math.system.Matrix;
import munoor.math.system.NativeMath;
import munoor.math.tools.MatrixDataProvider;
import sun.misc.Unsafe;

/**
 *
 * @author Ilamah, Osho
 */
public class Main {

    //private final static Logger LOGGER = Logger.getLogger(Matrix.class.getName());
    static {
        try {
            final String os = System.getProperty("os.name").toLowerCase();
            System.out.println("OS::" + os);
            final String libraryName = os.contains("win") ? "libNativeMath" : "NativeMath";
            System.loadLibrary(libraryName);
            //-XshowSettings:properties -version will tell you path
            //http://stackoverflow.com/questions/6092200/how-to-fix-an-unsatisfiedlinkerror-cant-find-dependent-libraries-in-a-jni-pro
            //http://stackoverflow.com/a/25186091
            //MatrixUtils.setThreads(2);
        } catch (UnsatisfiedLinkError e) {
            System.out.println("Native code library failed to load.\n");
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        //testCreate(); passing

        //testReadMatrix();
        testSymmetricEigs();
    }

    public static void testCreate() {
        final Matrix.DenseMatrix dm = Matrix.DenseMatrix.create(5, 5);
        if (dm != null) {
            System.out.println("Create matrix SUCCESS!!\n");
            //https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-blas-csr-matrix-storage-format
            dm.setElementByCoordinates(0, 0, 1);
            dm.setElementByCoordinates(0, 1, -1);
            dm.setElementByCoordinates(0, 3, -3);
            //
            dm.setElementByCoordinates(1, 0, -2);
            dm.setElementByCoordinates(1, 1, 5);
            //
            dm.setElementByCoordinates(2, 2, 4);
            dm.setElementByCoordinates(2, 3, 6);
            dm.setElementByCoordinates(2, 4, 4);
            //
            dm.setElementByCoordinates(3, 0, -4);
            dm.setElementByCoordinates(3, 2, 2);
            dm.setElementByCoordinates(3, 3, 7);
            //
            dm.setElementByCoordinates(4, 1, 8);
            dm.setElementByCoordinates(4, 4, -5);
            //
            dm.print(5, 5);
            //
            dm.destroy();
            System.out.println("Delete matrix SUCCESS!!\n");
        } else {
            System.out.println("Create matrix FAILED!!\n");
        }
        final Matrix.SSSMatrix sm = Matrix.SSSMatrix.create(8, 8, 10);
        //
        if (sm != null) {
            System.out.println("Create SSS matrix SUCCESS!!\n");
            //matrix from Gkountouvas
            final double[] dvalues = {2.7, 5.6, 9.4, 0.7, 2.4, 7.8, 9.8, 4.1};
            final double[] values = {0.5, 6.6, 3.1, 1.2, 9.8, 4.1, 7.2, 4.7, 3.4, 3.3};
            final int[] rowptr = {0, 0, 1, 2, 3, 5, 7, 8, 10};
            final int[] colind = {0, 1, 2, 0, 3, 2, 3, 5, 1, 4};
            //
            long diagonalValuesBaseAddress = sm.getDiagonalValuesBaseAddress();
            long rowPointerBaseAddress = sm.getRowPointerBaseAddress();
            long colIndexBaseAddress = sm.getColIndexBaseAddress();
            long valuesBaseAddress = sm.getValuesBaseAddress();

            try {
                Unsafe unsafe = Matrix.getUnsafe();
                for (int i = 0; i < 10; i++) {
                    unsafe.putDouble(valuesBaseAddress, values[i]);
                    unsafe.putInt(colIndexBaseAddress, colind[i]);

                    valuesBaseAddress += Matrix.DOUBLE_BYTES;
                    colIndexBaseAddress += Matrix.INT_BYTES;
                }
                //
                for (int i = 0; i < 8; i++) {
                    unsafe.putDouble(diagonalValuesBaseAddress, dvalues[i]);
                    unsafe.putInt(rowPointerBaseAddress, rowptr[i]);

                    diagonalValuesBaseAddress += Matrix.DOUBLE_BYTES;
                    rowPointerBaseAddress += Matrix.INT_BYTES;
                }
                //
                //extra row pointer
                unsafe.putInt(rowPointerBaseAddress, rowptr[8]);
                sm.print(18);
            } catch (Exception ex) {
                System.out.println("ERROR IN SSSMatrix set\n");
                ex.printStackTrace();
            }
            sm.destroy();
            System.out.println("Delete SSS matrix SUCCESS!!\n");
        } else {
            System.out.println("Create SSS matrix FAILED!!\n");
        }
    }

    public static void testReadMatrix() {
        //final String filename = "X:\\delete\\x\\Ragusa16\\LFAT5.mtx";
        final String filename = "/media/MERCURY/Samba/DEL/LFAT5.mtx";
        COOMatrix cooMatrix = MatrixDataProvider.readMatrixMarketFile(filename);
        if (cooMatrix == null) {
            System.out.println("Read COO Matrix FAILED!!");
        }
        final Matrix.DenseMatrix dm = MatrixDataProvider.getDenseMatrix(cooMatrix);
        final Matrix.SSSMatrix sm = MatrixDataProvider.getSSSMatrix(cooMatrix);

        dm.print(14, 14);
        sm.print(30);
        dm.destroy();
        sm.destroy();

        System.out.println("readMatrix passed, please check the printed entries!!");
    }

    public static void testSymmetricEigs() {
        //final String filename = "X:\\delete\\x\\Ragusa16\\LFAT5.mtx";
        final String sparseSymFilename = "/media/MERCURY/Samba/DEL/Ssym_mat_100.mtx";
        COOMatrix cooMatrix = MatrixDataProvider.readMatrixMarketFile(sparseSymFilename);
        if (cooMatrix == null) {
            System.out.println("Read COO Matrix FAILED!!");
            return;
        }
        System.out.println("=======SPARSE Symmetric EIGS ========================================================::");
        final Matrix.SSSMatrix sm = MatrixDataProvider.getSSSMatrix(cooMatrix);
        //determine symmetric eigs using native library
        final int k = 10;
        long outputRegisterPointer = NativeMath.allocateAllignedBuffer(2, (int) Matrix.LONG_BYTES, 64);
        final int sparseSymmResult = NativeMath.sssSparseSymmetricEigs(sm.getValuesBaseAddress(), sm.getDiagonalValuesBaseAddress(), sm.getRowPointerBaseAddress(), sm.getColIndexBaseAddress(), sm.getRows(), sm.getColumns(), sm.getLowerTriangleNonZeroCount(),
                k, 1 /* LA */, NativeMath.LA_EIGEN_ORDER, outputRegisterPointer);

        if (sparseSymmResult == 0) {
            System.out.println("Native sparse call completed without errors\nEigen VALUES::");
            //print eigen values
            long eigenValueAddress = sm.getLongElementByAddress(outputRegisterPointer);
            for (int i = 0; i < k; i++) {
                System.out.print(sm.getDoubleElementByAddress(eigenValueAddress) + "  ");
                eigenValueAddress += Matrix.LONG_BYTES;
            }
            System.out.println("\n\nEigen VECTORS (first K rows of the k columns)::");
            long eigenVectorsAddress = sm.getLongElementByAddress(outputRegisterPointer + Matrix.LONG_BYTES);
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    System.out.print(sm.getDoubleElementByAddress(eigenVectorsAddress) + "  ");
                    eigenVectorsAddress += Matrix.DOUBLE_BYTES;
                }
                System.out.println();
            }
            //Clean up the results
            NativeMath.releaseRegisters(outputRegisterPointer, 2);
        } else {
            System.out.println("Native sparse call FAILED!!!");
        }
        //
        System.out.println("=======DENSE Symmetric EIGS ========================================================::");
        final String denseSymFilename = "/media/MERCURY/Samba/DEL/Asym_mat_100.mtx";
        cooMatrix = MatrixDataProvider.readMatrixMarketFile(denseSymFilename);
        if (cooMatrix == null) {
            System.out.println("Read COO Matrix FAILED!!");
            return;
        }
        final Matrix.DenseMatrix dm = MatrixDataProvider.getDenseMatrix(cooMatrix);

        outputRegisterPointer = NativeMath.allocateAllignedBuffer(2, (int) Matrix.LONG_BYTES, 64);
        final int denseSymmResult = NativeMath.symmetricEigs(dm.getValuesBaseAddress(), dm.getRows(), dm.getColumns(), k, NativeMath.LA_EIGEN_ORDER, outputRegisterPointer);
        if (denseSymmResult == 0) {
            System.out.println("Native dense call completed without errors\nEigen VALUES::");
            //print eigen values
            long eigenValueAddress = dm.getLongElementByAddress(outputRegisterPointer);
            for (int i = 0; i < k; i++) {
                System.out.print(dm.getDoubleElementByAddress(eigenValueAddress) + "  ");
                eigenValueAddress += Matrix.LONG_BYTES;
            }
            System.out.println("\n\nEigen VECTORS (first K rows of the k columns)::");
            long eigenVectorsAddress = dm.getLongElementByAddress(outputRegisterPointer + Matrix.LONG_BYTES);
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    System.out.print(dm.getDoubleElementByAddress(eigenVectorsAddress) + "  ");
                    eigenVectorsAddress += Matrix.DOUBLE_BYTES;
                }
                System.out.println();
            }
            //Clean up the results
            NativeMath.releaseRegisters(outputRegisterPointer, 2);
        } else {
            System.out.println("Native dense call FAILED!!!");
        }
        dm.destroy();
        sm.destroy();
        //System.out.println("readMatrix passed, please check the printed entries!!");
    }

    public static void testSVD() {

    }
}
