package JamaEx;

import JamaEx.util.Maths;

/**
 * Jama = Java Matrix class.
 * 
 * Here we also provide a class called Matrix3, which has three dimensions of
 * double data. It has no inheritance relationship with Matrix, and it's also
 * not a wise way to convert a Mattrix3 object casted to a Matrix object
 * directly without using ctors.
 * 
 * Matrix3 is also not a set of Matrix. The total size of Matrix3 is m*n*d3.
 * 
 * Not tested yet.
 * 
 * @author Steven Chang
 * @version 20 June 2015
 */

public class Matrix3 implements Cloneable, java.io.Serializable {

	/*
	 * ------------------------ Class variables ------------------------
	 */

	/**
	 * Array for internal storage of elements.
	 * 
	 * @serial internal array storage.
	 */
	private double[][][] A;

	/**
	 * Row, column and the third dimension.
	 * 
	 * @serial row dimension.
	 * @serial column dimension.
	 * @serial the third dimension
	 */
	private int m, n, d3;

	/*
	 * ------------------------ Constructors ------------------------
	 */

	/**
	 * Construct an m-by-n-by-0 matrix of zeros.
	 * 
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of columns.
	 * @param d3
	 *            Number of the third dimension
	 */

	public Matrix3(int m, int n, int d3) {
		this.m = m;
		this.n = n;
		this.d3 = d3;
		A = new double[m][n][d3];
	}

	/**
	 * Construct an m-by-n-by-o constant matrix.
	 * 
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of columns.
	 * @param d3
	 *            Number of the third dimension
	 * @param s
	 *            Fill the matrix with this scalar value.
	 */

	public Matrix3(int m, int n, int o, double s) {
		this.m = m;
		this.n = n;
		A = new double[m][n][o];
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				for (int k = 0; k < d3; ++k) {
					A[i][j][k] = s;
				}
			}
		}
	}

	/**
	 * Construct a Matrix3 object from a 3-D array.
	 * 
	 * @param A
	 *            Three-dimensional array of doubles.
	 * @exception IllegalArgumentException
	 *                All rows must have the same length
	 * @see #constructWithCopy
	 */

	public Matrix3(double[][][] A) {
		m = A.length;
		n = A[0].length;
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException(
						"All rows must have the same length.");
			}
			for (int j = 0; j < n; ++j) {
				if (A[i][j].length != d3) {
					throw new IllegalArgumentException(
							"All rows must have the same length.");
				}
			}
		}
		this.A = A;
	}

	/**
	 * Construct a Matrix3 quickly without checking arguments.
	 * 
	 * @param A
	 *            Three-dimensional array of doubles.
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of columns.
	 * @param d3
	 *            Number of the third dimension
	 */

	public Matrix3(double[][][] A, int m, int n, int d3) {
		this.A = A;
		this.m = m;
		this.d3 = d3;
	}

	/**
	 * Construct a matrix3 from a one-dimensional packed array
	 * 
	 * @param vals
	 *            One-dimensional array of doubles, packed by columns (ala
	 *            Fortran).
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of columns.
	 * @exception IllegalArgumentException
	 *                Array length must be a multiple of m.
	 */

	public Matrix3(double vals[], int m, int n) {
		this.m = m;
		this.n = n;
		d3 = ((m != 0 && n != 0) ? vals.length / m / n : 0);
		if (m * n * d3 != vals.length) {
			throw new IllegalArgumentException(
					"Array length must be a multiple of m.");
		}
		A = new double[m][n][d3];
		for (int k = 0; k < d3; ++k) {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < m; i++) {
					A[i][j][k] = vals[i + j * m + k * m * n];
				}
			}
		}
	}

	/**
	 * Construct a matrix3 from a Matrix array. Matrix data will be set at row,
	 * col dimension.
	 * 
	 * @param matrixArray
	 *            One-dimensional array of Matrix
	 * @exception IllegalArgumentException
	 *                width and length of each element in array must be the
	 *                same.
	 */

	public Matrix3(Matrix[] matrixArray) {

		this.m = matrixArray[0].getRowDimension();
		this.n = matrixArray[0].getColumnDimension();
		for (int i = 1; i < matrixArray.length; ++i) {
			if (m != matrixArray[i].getRowDimension()
					|| n != matrixArray[i].getColumnDimension()) {
				throw new IllegalArgumentException(
						"width and length of each element in array must be the same.");
			}
		}
		// Copy data
		for (int k = 0; k < d3; ++k) {
			Matrix src = matrixArray[k];
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					A[i][j][k] = src.get(i, j);
				}
			}
		}
	}

	/*
	 * ------------------------ Public Methods ------------------------
	 */

	/**
	 * Generate matrix with random elements
	 * 
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of colums.
	 * @param d3
	 *            Number of the third dimension.
	 * @return An m-by-n-by-r3 matrix with uniformly distributed random
	 *         elements.
	 */

	public static Matrix3 random(int m, int n, int d3) {
		Matrix3 A = new Matrix3(m, n, d3);
		double[][][] X = A.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < d3; ++k) {
					X[i][j][k] = Math.random();
				}
			}
		}
		return A;
	}

	/**
	 * Construct a matrix3 from a copy of a 3-D array.
	 * 
	 * @param A
	 *            Three-dimensional array of doubles.
	 * @exception IllegalArgumentException
	 *                All rows must have the same length
	 */

	public static Matrix3 constructWithCopy(double[][][] A) {
		int m = A.length;
		int n = A[0].length;
		int d3 = A[0][0].length;

		// Check
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException(
						"All rows must have the same length.");
			}
			for (int j = 0; j < n; ++j) {
				if (A[i][j].length != d3) {
					throw new IllegalArgumentException(
							"All rows must have the same length.");
				}
			}
		}

		Matrix3 X = new Matrix3(m, n, d3);
		// Copy data
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				for (int k = 0; k < d3; ++k) {
					X.set(i, j, k, A[i][j][k]);
				}
			}
		}
		return X;
	}

	/**
	 * Make a deep copy of a matrix
	 */

	public Matrix3 copy() {
		Matrix3 X = new Matrix3(m, n, d3);
		double[][][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < d3; ++k) {
					C[i][j] = A[i][j];
				}
			}
		}
		return X;
	}

	/**
	 * Clone the Matrix3 object.
	 */

	public Object clone() {
		return this.copy();
	}

	/**
	 * Access the internal three-dimensional array.
	 * 
	 * @return Pointer to the three-dimensional array of matrix elements.
	 */

	public double[][][] getArray() {
		return A;
	}

	/**
	 * Get row dimension.
	 * 
	 * @return m, the number of rows.
	 */

	public int getRowDimension() {
		return m;
	}

	/**
	 * Get column dimension.
	 * 
	 * @return n, the number of columns.
	 */

	public int getColumnDimension() {
		return n;
	}

	/**
	 * Get the third dimension.
	 * 
	 * @return d3, the number of the third dimension.
	 */

	public int getThirdDimension() {
		return d3;
	}

	/**
	 * Get a single element.
	 * 
	 * @param i
	 *            Row index.
	 * @param j
	 *            Column index.
	 * @param k
	 *            The third dimension index
	 * @return A(i, j, k)
	 * @exception ArrayIndexOutOfBoundsException
	 */

	public double get(int i, int j, int k) {
		return A[i][j][k];
	}

	/**
	 * Get a single element from the index. Search order: m-n-d3
	 * 
	 * @param index
	 *            Row index. Start from 0.
	 * @return A(i,j, k)
	 * @exception ArrayIndexOutOfBoundsException
	 */

	public double get(int index) {
		int row = index / m / n;
		int val = (index - row) / m;
		int col = val % n;
		int dim3 = val / n;
		return A[row][col][dim3];
	}

	/**
	 * Get a submatrix.
	 * 
	 * @param i0
	 *            Initial row index
	 * @param i1
	 *            Final row index
	 * @param j0
	 *            Initial column index
	 * @param j1
	 *            Final column index
	 * @return A(i0:i1,j0:j1)
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public Matrix3 getMatrix3(int i0, int i1, int j0, int j1, int k0, int k1) {
		Matrix3 X = new Matrix3(i1 - i0 + 1, j1 - j0 + 1, k1 - k0 + 1);
		double[][][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					for (int k = k0; k <= k1; ++k) {
						B[i - i0][j - j0][k - k0] = A[i][j][k];
					}
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/**
	 * Set a single element.
	 * 
	 * @param i
	 *            Row index.
	 * @param j
	 *            Column index.
	 * @param k
	 *            The third dimension index
	 * @param s
	 *            A(i,j,k).
	 * @exception ArrayIndexOutOfBoundsException
	 */

	public void set(int i, int j, int k, double s) {
		A[i][j][k] = s;
	}

	/**
	 * Set a single element.
	 * 
	 * @param index
	 *            index of elements
	 * @param value
	 *            value to be set A(i,j,k).
	 * @exception ArrayIndexOutOfBoundsException
	 */

	public void set(int index, double value) {
		int row = index / m / n;
		int val = (index - row) / m;
		int col = val % n;
		int dim3 = val / n;
		A[row][col][dim3] = value;
	}

	/**
	 * Get the a Matrix from specified dimension.
	 * 
	 * @param index
	 *            index of the specified dimension
	 * @param dim
	 *            Dimension index, which index searched.
	 */

	public Matrix getMatrix(int index, int dim) throws Exception {

		if (dim == 1) {
			// Row dimension
			Matrix result = Matrix.constructWithCopy(A[index]);
			return result;
		} else if (dim == 2) {
			// Column dimension
			Matrix result = new Matrix(m, d3);
			for (int i = 0; i < m; ++i) {
				for (int k = 0; k < d3; ++k) {
					result.set(i, k, A[i][dim][k]);
				}
			}
			return result;
		} else if (dim == 3) {
			// the third dimension
			Matrix result = new Matrix(m, n);
			// Copy the data
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					result.set(i, j, A[i][j][dim]);
				}
			}
			return result;
		} else {
			throw new IllegalArgumentException("Dimension must be 1, 2 or 3");
		}
	}

	/**
	 * Return the size of elements. m*n*d3.
	 * 
	 */
	public int elementSize() {
		return m * n * d3;
	}

	/**
	 * Perform abs() method for all elements in the matrix.
	 * 
	 */

	public Matrix3 abs() throws Exception {

		double[][][] B = new double[m][n][d3];
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < m; ++j) {
				for (int k = 0; k < d3; ++k) {
					B[i][j][k] = A[i][j][k];
					if (B[i][j][k] < 0) {
						B[i][j][k] *= -1.;
					}
				}
			}
		}
		return new Matrix3(B);
	}

	/**
	 * Unary minus
	 * 
	 * @return -A
	 */

	public Matrix3 uminus() {
		Matrix3 X = this.copy();
		double[][][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < d3; ++k) {
					C[i][j][k] = -C[i][j][k];
				}
			}
		}
		return X;
	}

	/*
	 * ------------------------ Private Methods ------------------------
	 */

	/** Check if size(A) == size(B) **/

	private void checkMatrixDimensions(Matrix3 B) {
		if (B.m != m || B.n != n || B.d3 != d3) {
			throw new IllegalArgumentException("Matrix dimensions must agree.");
		}
	}

	private static final long serialVersionUID = 1;
}