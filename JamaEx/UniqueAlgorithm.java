package JamaEx;

/**
 * Jama = Java Matrix class.
 * 
 * This class is trying to implement the algorithm of unique in MatLab.
 * 
 * @author Steven Chang
 * @version 21 June 2015
 */

public class UniqueAlgorithm implements java.io.Serializable {

	/*
	 * ------------------------ Class variables ------------------------
	 */

	/**
	 *  Returns the same data as in A, but with no repetitions.
	 * 
	 */
	private Matrix C;

	/**
	 * Index vectors ia and ic, such that C = A(ia,:) and A = C(ic,:).
	 */
	private Matrix ia, ic;

	/*
	 * ------------------------ Constructors ------------------------
	 */
	
	/**
	 * Construct an empty Unique object.
	 * 
	 */
	public UniqueAlgorithm() {
		C = null;
		ia = null;
		ic = null;
	}
	
	/**
	 * Construct and perform unique method.
	 * As [C,ia,ic] = unique(A) in MatLab.
	 * 
	 * @param A
	 * 		  Source data matrix
	 * @throws Exception 
	 */

	public UniqueAlgorithm(Matrix A) throws Exception {
		unique(A);
	}
	
	

	/*
	 * ------------------------ Public Methods ------------------------
	 */
	
	/**
	 * Perform unique method.
	 * As [C,ia,ic] = unique(A) in MatLab.
	 * 
	 * @param A
	 * 		  Source data matrix
	 * @throws Exception 
	 */

	public void unique(Matrix A) throws Exception {
		
		C = A.sort();
		int num = 1;
		double currVal = A.get(0);
		for(int i=1; i<C.size(); ++i) {
			if(currVal != C.get(i)) {
				currVal = C.get(i);
				C.set(num, currVal);
				++num;
			}
		}
		// Get C
		C = C.getMatrix(0, num);
		// Get ia
		ia = new Matrix(C.size(), 1);
		for(int i=0; i<C.size(); ++i) {
			int index = A.find_first(C.get(i));
			ia.set(i, 0, index);
		}
		// Get ic
		ic = A.copy();
		for(int i=0; i<C.size(); ++i) {
			ic = ic.equalsSubstitute(C.get(i), i);
		}
	}
	
	/**
	 * Perform unique method. It will treat each row of A as a single entity and
	 * return the unique rows of A. The rows of the array C are in sorted order.
	 * As [C,ia,ic] = unique(A,'rows') in MatLab.
	 * 
	 * @param A
	 * 		  Source data matrix
	 * @throws Exception 
	 */
	public void unique_rows(Matrix A) throws Exception {
		
		Matrix sorted = A.sort(1);
		C = new Matrix(1, A.getColumnDimension());
		Matrix currVal = A.getRow(0);
		C = A.getRow(0);
		for(int i=1; i<sorted.getRowDimension(); ++i) {
			for(int j=0; j<currVal.getColumnDimension(); ++j) {
				if(currVal.get(j) != sorted.get(i, j)) {
					currVal = sorted.getRow(i);
					C = Matrix.concatenate(C, currVal, 2);
					break;
				}
			}
		}
		// Get ia
		ia = new Matrix(C.getRowDimension(), 1);
		for(int i=0; i<C.getRowDimension(); ++i) {
			int index = A.find_first_row(C.getRow(i));
			ia.set(i, index);
		}
		// Get ic
		ic = new Matrix(A.getRowDimension(), 1);
		for(int i=0; i<A.getRowDimension(); ++i) {
			int index = C.find_first_row(A.getRow(i));
			ic.set(i, index);
		}
	}
	
	/**
	 * Get C
	 * 
	 * @throws Exception 
	 */
	
	public Matrix getC() {
		if (C==null) {
			throw new RuntimeException("Member matrix is not initialized");
		}
		return C;
	}
	
	/**
	 * Get ia
	 * 
	 * @throws Exception 
	 */
	
	public Matrix getia() {
		if (ia==null) {
			throw new RuntimeException("Member matrix is not initialized");
		}
		return ia;
	}
	
	/**
	 * Get ic
	 * 
	 * @throws Exception 
	 */
	
	public Matrix getic() {
		if (ia==null) {
			throw new RuntimeException("Member matrix is not initialized");
		}
		return ic;
	}
	
	
	/*
	 * ------------------------ Private Methods ------------------------
	 */

	private static final long serialVersionUID = 1;
}