/*
 * 2017.05.17 - Created
 * 2017.06.15 - Switch to use of native allocate and destroy
 * 2018.04.30 - Introduced dense and SS sparse specializations
 */
package munoor.math.system;

import java.lang.reflect.Field;
import java.text.NumberFormat;
import sun.misc.Unsafe;

/**
 * Java implementation of a matrix of double precision numbers with off-jvm-heap data storage structures. The off-jvm-heap storage structures enable large matrices to be supported and provides more
 * efficient coupling with native math kernels. The backing memory is managed using the Unsafe API. This abstract class provides attributes/methods
 *
 * @author Ilamah, Osho
 */
/*
http://mydailyjava.blogspot.co.uk/2013/12/sunmiscunsafe.html
http://stackoverflow.com/questions/9323416/using-memory-allocated-by-sun-misc-unsafe-allocatememory-in-native-code
http://mishadoff.com/blog/java-magic-part-4-sun-dot-misc-dot-unsafe/
https://www.mkyong.com/java/java-write-directly-to-memory/
//
//
https://github.com/bytedeco/javacpp-presets
 */
public abstract class Matrix {

    //private final static Logger LOGGER = Logger.getLogger(Matrix.class.getName());
    public static final long DOUBLE_BYTES = 8, LONG_BYTES = 8, INT_BYTES = 4;
    private static final int DEFAULT_ALIGNMENT = 64;
    protected final int rows, columns;
    protected final long valuesBaseAddress;

    public final int getRows() {
        return rows;
    }

    public final int getColumns() {
        return columns;
    }

    public final long getValuesBaseAddress() {
        return valuesBaseAddress;
    }

    public static void setDoubleElementByAddress(final long address, final double value) {
        //NO bound checking, I hope the caller knows what they are doing?
        unsafe.putDouble(address, value);
    }

    public static final double getDoubleElementByAddress(final long address) {
        //NO bound checking, I hope the caller knows what they are doing?
        return unsafe.getDouble(address);
    }

    public static void setIntElementByAddress(final long address, final int value) {
        //NO bound checking, I hope the caller knows what they are doing?
        unsafe.putInt(address, value);
    }

    public static final int getIntElementByAddress(final long address) {
        //NO bound checking, I hope the caller knows what they are doing?
        return unsafe.getInt(address);
    }

    public static void setLongElementByAddress(final long address, final long value) {
        //NO bound checking, I hope the caller knows what they are doing?
        unsafe.putLong(address, value);
    }

    public static final long getLongElementByAddress(final long address) {
        //NO bound checking, I hope the caller knows what they are doing?
        return unsafe.getLong(address);
    }

//    public abstract void setElementByCoordinates(final int i, final int j, final double value);
    //public abstract void print(final int maxRow, final int maxCol);
    //package constructor!! blocked from outside world
    Matrix(int rows, int columns, long baseAddress) {
        this.rows = rows;
        this.columns = columns;
        this.valuesBaseAddress = baseAddress;
    }

    public abstract void destroy();

    /**
     * Direct memory (off-heap) implementation of a row-major dense matrix of double precision numbers. The backing memory is managed using the Unsafe API
     *
     * @author Ilamah, Osho
     */
    public static class DenseMatrix extends Matrix {

        private DenseMatrix(int rows, int columns, long baseAddress) {
            super(rows, columns, baseAddress);
        }

        /**
         * Retrieves the matrix element at the specified coordinates. The coordinates are zero based.
         *
         * @param i
         * @param j
         * @return
         */
        public final double getElementByCoordinates(int i, int j) {
            //NO bound checking, I hope the caller knows what they are doing?
            //return M->data[row_num * (M->ncols) + col_num];
            final long address = (((i * columns) + j) * DOUBLE_BYTES) + valuesBaseAddress;
            return unsafe.getDouble(address);
        }

        /**
         * Sets the matrix element at the specified coordinates. The coordinates are zero based.
         *
         * @param i
         * @param j
         * @param value
         */
        public final void setElementByCoordinates(int i, int j, double value) {
            //NO bound checking, I hope the caller knows what they are doing?
            final long address = (((i * columns) + j) * DOUBLE_BYTES) + valuesBaseAddress;
            unsafe.putDouble(address, value);
        }

        public final static DenseMatrix create(int m, int n) {
            try {
                if (m > 0 && n > 0) {
                    getUnsafe();
                    //final long baseAddress = unsafe.allocateMemory(size);
                    final long baseAddress = NativeMath.allocateAllignedBuffer(((long) m) * n /* ensure long multiplication  */, (int) DOUBLE_BYTES, DEFAULT_ALIGNMENT);//    allocateAlligned(size, DEFAULT_ALIGNMENT);
                    if (baseAddress == 0) {
                        throw new RuntimeException("Native Matrix could not allocate sufficient memory!!");
                    }
                    final DenseMatrix nm = new DenseMatrix(m, n, baseAddress);
                    return nm;
                } else {
                    throw new IllegalArgumentException("Rows and columns must be > 0");
                }
            } catch (Exception e) {
                return null;
            }
        }

        public void print(final int maxRow, final int maxCol) {
            final int ROW_LIMIT = Math.min(Math.max(1, maxRow), rows);
            final int COL_LIMIT = Math.min(Math.max(1, maxCol), columns);

            NumberFormat nf = NumberFormat.getNumberInstance();
            nf.setMaximumFractionDigits(3);
            nf.setMinimumFractionDigits(3);
            nf.setGroupingUsed(false);

            System.out.print("\n\n");
            for (int row = 0; row < ROW_LIMIT; row++) {
                for (int col = 0; col < COL_LIMIT; col++) {
                    System.out.print(nf.format(getElementByCoordinates(row, col)) + " ");
                }
                System.out.print("\n");
            }
            System.out.print("\n\n");
        }

        @Override
        public void destroy() {
            NativeMath.releaseAllignedBuffer(valuesBaseAddress);
        }
    }

    public static class SSSMatrix extends Matrix {

        private final int nnzl;
        protected final long diagonalValuesBaseAddress, colIndexBaseAddress, rowPointerBaseAddress;

        private SSSMatrix(final int rows, final int columns, final int nnzl, final long baseAddress, final long diagonalValuesBaseAddress, final long colIndexBaseAddress, final long rowPointerBaseAddress) {
            super(rows, columns, baseAddress);
            this.nnzl = nnzl;
            this.diagonalValuesBaseAddress = diagonalValuesBaseAddress;
            this.colIndexBaseAddress = colIndexBaseAddress;
            this.rowPointerBaseAddress = rowPointerBaseAddress;
        }

        public final static SSSMatrix create(final int m, final int n, final int nnzl) {
            try {
                if (m > 0 && n > 0) {
                    getUnsafe();
                    final int[] bytesPerElement = {(int) DOUBLE_BYTES, (int) DOUBLE_BYTES, (int) INT_BYTES, (int) INT_BYTES};
                    final long[] elementsPerBuffer = {m, nnzl, nnzl, m + 1};
                    //cross layer once!
                    final long[] addresses = NativeMath.allocateAllignedBuffers(elementsPerBuffer, bytesPerElement, DEFAULT_ALIGNMENT);
                    for (int i = 0; i < addresses.length; i++) {
                        if (addresses[i] == 0) {
                            throw new RuntimeException("Native SSSMatrix could not allocate sufficient memory!!");
                        }
                    }
                    final long dValuesBaseAddress = addresses[0];
                    final long baseAddress = addresses[1];
                    final long colIndexBaseAddress = addresses[2];
                    final long rowPointerBaseAddress = addresses[3];
                    final SSSMatrix sss = new SSSMatrix(m, n, nnzl, baseAddress, dValuesBaseAddress, colIndexBaseAddress, rowPointerBaseAddress);
                    return sss;
                } else {
                    throw new IllegalArgumentException("Rows and columns must be > 0");
                }
            } catch (Exception e) {
                return null;
            }
        }

        public final void print(final int n) {
            final int entries = Math.min(Math.max(1, n), nnzl + rows);
            final String message = "\n(%d, %d) -> %g";
            long diagonalValueIndexAddress = diagonalValuesBaseAddress;
            for (int x = 0, row = 0; x < entries; row++) {
                final long rowPointerAddress = (row * INT_BYTES) + rowPointerBaseAddress;
                final int rowPointerStart = unsafe.getInt(rowPointerAddress);
                final int rowPointerEnd = unsafe.getInt(rowPointerAddress + INT_BYTES);
                long colIndexAddress = (rowPointerStart * INT_BYTES) + colIndexBaseAddress;
                long valueIndexAddress = (rowPointerStart * DOUBLE_BYTES) + valuesBaseAddress;
                final int elementsInRow = rowPointerEnd - rowPointerStart;
                for (int ii = rowPointerStart; ii - rowPointerStart < elementsInRow && x < entries; ii++) {
                    final int col = unsafe.getInt(colIndexAddress);
                    final double value = unsafe.getDouble(valueIndexAddress);
                    System.out.print(String.format(message, row, col, value));
                    colIndexAddress += INT_BYTES;
                    valueIndexAddress += DOUBLE_BYTES;
                    x++;
                }
                if (x < entries) {
                    final double value = unsafe.getDouble(diagonalValueIndexAddress);
                    diagonalValueIndexAddress += DOUBLE_BYTES;
                    System.out.print(String.format(message, row, row, value));
                    x++;
                }

            }
            System.out.println("\n\n");
        }

        public final int getLowerTriangleNonZeroCount() {
            return nnzl;
        }

        public final long getDiagonalValuesBaseAddress() {
            return diagonalValuesBaseAddress;
        }

        public final long getColIndexBaseAddress() {
            return colIndexBaseAddress;
        }

        public final long getRowPointerBaseAddress() {
            return rowPointerBaseAddress;
        }

        @Override
        public void destroy() {
            final long[] pointers = {valuesBaseAddress, diagonalValuesBaseAddress, colIndexBaseAddress, rowPointerBaseAddress};
            NativeMath.releaseAllignedBuffers(pointers);
        }

    }

    protected static Unsafe unsafe;

    public static Unsafe getUnsafe() throws Exception {
        if (unsafe == null) {
            // Get the Unsafe object instance
            Field field = sun.misc.Unsafe.class.getDeclaredField("theUnsafe");
            field.setAccessible(true);
            unsafe = (sun.misc.Unsafe) field.get(null);
        }
        return unsafe;
    }
}
