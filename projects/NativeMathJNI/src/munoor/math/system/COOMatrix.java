/*
 * 2018.05.02  - Created
 */
package munoor.math.system;

import java.util.ArrayList;

/**
 *
 * @author Ilamah, Osho
 */
public class COOMatrix {
    private final int rows, cols;
    private int nnz = -1, nnzd = -1, nnzl = -1, nnzr = -1;
    //
    private final boolean sparse, symmetric;
    
    final ArrayList<COOEntry> entries = new ArrayList<>();

    public COOMatrix(int rows, int cols, boolean sparse, boolean symmetric) {
        this.rows = rows;
        this.cols = cols;
        this.sparse = sparse;
        this.symmetric = symmetric;
    }

    public int getTotalNumberOfNonZeros() {
        return nnz;
    }

    public void setTotalNumberOfNonZeros(int nnz) {
        this.nnz = nnz;
    }

    public int getTotalOffDiagonalNonZeros() {
        return nnzd;
    }

    public void setTotalOffDiagonalNonZeros(int nnzd) {
        this.nnzd = nnzd;
    }

    public int getTotalLowerTriangleNonZeros() {
        return nnzl;
    }

    public void setTotalLowerTriangleNonZeros(int nnzl) {
        this.nnzl = nnzl;
    }

    public int getTotalReadRecords() {
        return nnzr;
    }

    public void setTotalReadRecords(int nnzr) {
        this.nnzr = nnzr;
    }

    public int getRows() {
        return rows;
    }

    public int getCols() {
        return cols;
    }

    public boolean isSymmetric() {
        return symmetric;
    }

    public boolean isSparse() {
        return sparse;
    }
    
    

    public ArrayList<COOEntry> getEntries() {
        return entries;
    }
    
    public static class COOEntry {
        public int i, j;
        public double value;

        public COOEntry(int i, int j, double value) {
            this.i = i;
            this.j = j;
            this.value = value;
        }
        
    }  
}
