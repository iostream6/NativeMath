/*
 * 2018.05.01  - Created
 */
package munoor.math.tools;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import munoor.math.system.COOMatrix;
import munoor.math.system.Matrix;
import munoor.math.system.Matrix.DenseMatrix;
import munoor.math.system.Matrix.SSSMatrix;

/**
 * Provides read support for matrix data file formats (e.g. Matrix Market .mtx text files) as well as conversion between various matrix structures (dense, sparse, sparse symmetric)
 *
 * @author Ilamah, Osho
 */
public class MatrixDataProvider {

    /**
     * Reads a Matrix Market .mtx file
     *
     * @param inputFilePath the absolute file path to the .mtx file
     * @return if successful, a COOMatrix structure containing the read data, otherwise null
     */
    public static COOMatrix readMatrixMarketFile(final String inputFilePath) {
        final String COMMENT_MARKER = "%";
        final File inputFile = new File(inputFilePath);
        final String INVALID_INPUT_FILE_MESSAGE = "MatrixMarketReader :: invalid input file";
        if (!inputFile.exists() || inputFile.isDirectory()) {
            System.out.println(INVALID_INPUT_FILE_MESSAGE);
            return null;
        }
        try (final Scanner scanner = new Scanner(inputFile)) {
            final boolean symmetric;
            final FORMAT format;
            final DATA_TYPE dataType;
            if (scanner.hasNextLine()) {
                final String firstLine = scanner.nextLine();
                final String[] controlTokens = firstLine.split("\\s+");
                if (controlTokens.length != 5) {
                    System.out.println(INVALID_INPUT_FILE_MESSAGE);//invalid header/control line
                    return null;
                }
                //
                final String HEADER_PREFIX = "%%MatrixMarket";
                final String MATRIX_FLAG = "matrix", NON_SYMMETRIC_FLAG = "general", SYMMETRIC_FLAG = "symmetric";
                //sanity checks
                final String INVALID_CONTROL_MESSAGE = "MatrixMarketReader :: invalid MM control token";
                if (!controlTokens[0].equals(HEADER_PREFIX) || !controlTokens[1].equals(MATRIX_FLAG)) {
                    System.out.println(INVALID_CONTROL_MESSAGE);
                    return null;
                }
                //
                try {
                    format = FORMAT.fromString(controlTokens[2]);
                    dataType = DATA_TYPE.fromString(controlTokens[3]);
                } catch (IllegalArgumentException iae) {
                    System.out.println(INVALID_CONTROL_MESSAGE);
                    return null;
                }
                switch (controlTokens[4]) {
                    case SYMMETRIC_FLAG:
                        symmetric = true;
                        break;
                    case NON_SYMMETRIC_FLAG:
                        symmetric = false;
                        break;
                    default:
                        System.out.println(INVALID_CONTROL_MESSAGE);
                        return null;
                }
            } else {
                System.out.println(INVALID_INPUT_FILE_MESSAGE); // Input file is empty
                return null;
            }
            //
            int ROWS = 0, COLUMNS = 0, MM_DATA_COUNT = 0;
            while (scanner.hasNextLine()) {
                final String secondLine = scanner.nextLine();
                if (secondLine.startsWith(COMMENT_MARKER)) {
                    continue;
                }
                //process second line
                String[] dimensions = secondLine.split("\\s+");
                if (dimensions.length < 2 || (format == FORMAT.COORDINATE && dimensions.length < 3)) {
                    System.out.println(INVALID_INPUT_FILE_MESSAGE);//invalid header/dimens
                    return null;
                }
                try {
                    ROWS = Integer.parseInt(dimensions[0]);
                    COLUMNS = Integer.parseInt(dimensions[1]);
                    if (format == FORMAT.COORDINATE) {
                        MM_DATA_COUNT = Integer.parseInt(dimensions[2]);
                    }
                } catch (Exception e) {
                    // WIll be handled below
                }
                break;
            }
            if (ROWS <= 0 || COLUMNS <= 0 || (format == FORMAT.COORDINATE && MM_DATA_COUNT <= 0)) {
                System.out.println(INVALID_INPUT_FILE_MESSAGE);//invalid header/dimens
                return null;
            }
            //final ArrayList<COOMatrix.COOEntry> mmRecords = new ArrayList<>();
            final COOMatrix cooMatrix = new COOMatrix(ROWS, COLUMNS, format == FORMAT.COORDINATE, symmetric);
            int nnzd = 0, nnzl = 0, nnzr = 0;
            try {
                if (format == FORMAT.COORDINATE) {
                    while (scanner.hasNextLine()) {
                        final String line = scanner.nextLine();
                        if (line.startsWith(COMMENT_MARKER)) {
                            continue;
                        }
                        //process line
                        String[] records = line.split("\\s+");
                        int row = Integer.parseInt(records[0]);
                        int col = Integer.parseInt(records[1]);
                        final double realValue;
                        switch (dataType) {
                            case PATTERN:
                                realValue = 1;
                                break;
                            default: // REAL, INTEGER
                                realValue = Double.parseDouble(records[2]);
                        }
                        final COOMatrix.COOEntry cooEntry = new COOMatrix.COOEntry(--row, --col, realValue);//Note decrement by 1 to get back to Java array index basis                       
                        cooMatrix.getEntries().add(cooEntry);
                        nnzr++;
                        if (symmetric) {
                            if (row == col) {
                                //diagonal element
                                nnzd++;
                            } else if (row > col) {
                                nnzl++;
                            }
                        }
                    }
                } else {
                    //ARRAY format!! TODO
                    double realValue;
                    //the below assumes there are no comment lines after the control info line
                    for (int j = 0; j < COLUMNS; j++) {//column order is how dense MM files are stored
                        for (int i = symmetric ? j : 0; i < ROWS; i++) {/*lower triangle only for symmetric case */
                            final String valueString = scanner.nextLine().trim();
                            realValue = Double.parseDouble(valueString);
                            final COOMatrix.COOEntry cooEntry = new COOMatrix.COOEntry(i, j, realValue);
                            cooMatrix.getEntries().add(cooEntry);
                            nnzr++;
                            if (i == j) {
                                //diagonal element
                                nnzd++;
                            } else if (i > j) {
                                nnzl++;
                            }
                        }
                    }
                    final int expectedEntries = ((ROWS * (ROWS + 1))/2);
                    if((symmetric && (nnzr != expectedEntries)) || (!symmetric && (nnzr != (ROWS * COLUMNS)))){
                        System.out.println("MatrixMarketReader :: insufficient data entries encountered");
                        return null;
                    }
                }
                //sort the COO entries in row major format (i.e. by row and then by column)
                cooMatrix.getEntries().sort((a, b) -> {
                    final int rowPriority = Integer.compare(a.i, b.i);
                    if (rowPriority == 0) {
                        return Integer.compare(a.j, b.j);
                    } else {
                        return rowPriority;
                    }
                });
                if (symmetric) {
                    //
                    cooMatrix.setTotalLowerTriangleNonZeros(nnzl);
                    cooMatrix.setTotalOffDiagonalNonZeros(nnzd);
                } else {
                    cooMatrix.setTotalNumberOfNonZeros(nnzr);
                }
                cooMatrix.setTotalReadRecords(nnzr);
                return cooMatrix;
            } catch (Exception e) {
                System.out.println("MatrixMarketReader :: data parse error");
                return null;
            }
        } catch (IOException ioe) {
            System.out.println("MatrixMarketReader :: data parse error");
            return null;
        }

    }

    public static DenseMatrix getDenseMatrix(final COOMatrix matrix) {
        final DenseMatrix dm = DenseMatrix.create(matrix.getRows(), matrix.getCols());
        if (dm != null) {
            if (matrix.isSymmetric()) {
                matrix.getEntries().stream().forEach(entry -> {
                    dm.setElementByCoordinates(entry.i, entry.j, entry.value);
                    if (entry.i != entry.j) {
                        dm.setElementByCoordinates(entry.j, entry.i, entry.value);
                    }
                });
            } else {
                matrix.getEntries().stream().forEach(entry -> {
                    dm.setElementByCoordinates(entry.i, entry.j, entry.value);
                });
            }
            return dm;
        }
        return null;
    }

    public static SSSMatrix getSSSMatrix(final COOMatrix matrix) {
        if (matrix.isSymmetric() && matrix.isSparse()) {
            final SSSMatrix sm = SSSMatrix.create(matrix.getRows(), matrix.getCols(), matrix.getTotalLowerTriangleNonZeros());
            if (sm != null) {
                // compute number of non-zero, off diagonal lower triangular entries per row of A, this is achieved by scanning all elements of 
                // the COO and incrementing row_ptr for the row containing each element
                // after this stage, each entry in row_ptr  is the number of elements in the corresponding LT row of the COO matrix
                final long diagonalValuesBaseAddress = sm.getDiagonalValuesBaseAddress();
                final long rowPointerBaseAddress = sm.getRowPointerBaseAddress();
                matrix.getEntries().stream().forEach(entry -> {
                    if (entry.i == entry.j) {
                        //diagonal element
                        sm.setDoubleElementByAddress(diagonalValuesBaseAddress + (Matrix.DOUBLE_BYTES * entry.i), entry.value);
                    } else if (entry.i > entry.j) {
                        //lower trianglular element
                        final long address = rowPointerBaseAddress + (Matrix.INT_BYTES * entry.i);
                        sm.setIntElementByAddress(address, sm.getIntElementByAddress(address) + 1);
                    }
                });
                // We now convert the row_ptr to the actual SSS row_ptr definition - running sum which defines the number of stored elements in the previous row
                // this is achieved by computing the running sum, starting at 0
                for (int row = 0, runningSum = 0; row < matrix.getRows(); row++) {
                    final long address = rowPointerBaseAddress + (Matrix.INT_BYTES * row);
                    int temp =  sm.getIntElementByAddress(address);
                    sm.setIntElementByAddress(address, runningSum);
                    runningSum += temp;
                }
                // fill in the nrows+1 row_ptr value, which must be nnzl (not nnz!!)
                sm.setIntElementByAddress(rowPointerBaseAddress + (Matrix.INT_BYTES * matrix.getRows()), matrix.getTotalLowerTriangleNonZeros());
                //
                for (int n = 0; n < matrix.getTotalReadRecords(); n++) {
                    final COOMatrix.COOEntry entry = matrix.getEntries().get(n);
                    //
                    if (entry.i <= entry.j) {
                        //diagonal element
                        continue;
                    }
                    //lower triangle element, please store
                    final long rowPointerAddress = rowPointerBaseAddress + (Matrix.INT_BYTES * entry.i);
//                   sm.setIntElementByAddress(address, sm.getIntElementByAddress(address) + 1);
                    final int slot = sm.getIntElementByAddress(rowPointerAddress);

                    sm.setIntElementByAddress((sm.getColIndexBaseAddress() + (Matrix.INT_BYTES * slot)), entry.j);
                    sm.setDoubleElementByAddress((sm.getValuesBaseAddress() + (Matrix.DOUBLE_BYTES * slot)), entry.value);
//                   
                    sm.setIntElementByAddress(rowPointerAddress, slot + 1);
                }
                for (int row = 0, last = 0; row < matrix.getRows(); row++) {
                    final long address = rowPointerBaseAddress + (Matrix.INT_BYTES * row);
                    int temp =  sm.getIntElementByAddress(address);
                    sm.setIntElementByAddress(address, last);
                    last = temp;
                }
                return sm;
            }
        }
        return null;
    }

    private static enum FORMAT {
        ARRAY, COORDINATE;

        static FORMAT fromString(String token) {
            final String TOKEN = token.toUpperCase();
            for (final FORMAT f : values()) {
                if (f.name().equals(TOKEN)) {
                    return f;
                }
            }
            throw new IllegalArgumentException();
        }
    }

    private static enum DATA_TYPE {
        INTEGER, PATTERN, REAL;

        static DATA_TYPE fromString(String token) {
            final String TOKEN = token.toUpperCase();
            for (final DATA_TYPE dt : values()) {
                if (dt.name().equals(TOKEN)) {
                    return dt;
                }
            }
            throw new IllegalArgumentException();
        }
    }
}
