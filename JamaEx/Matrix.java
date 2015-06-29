package JamaEx;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.StreamTokenizer;

import JamaEx.util.Maths;

/**
 * Jama = Java Matrix class.
 * <P>
 * The Java Matrix Class provides the fundamental operations of numerical linear
 * algebra. Various constructors create Matrices from two dimensional arrays of
 * double precision floating point numbers. Various "gets" and "sets" provide
 * access to submatrices and matrix elements. Several methods implement basic
 * matrix arithmetic, including matrix addition and multiplication, matrix
 * norms, and element-by-element array operations. Methods for reading and
 * printing matrices are also included. All the operations in this version of
 * the Matrix Class involve real matrices. Complex matrices may be handled in a
 * future version.
 * <P>
 * Five fundamental matrix decompositions, which consist of pairs or triples of
 * matrices, permutation vectors, and the like, produce Xs in five decomposition
 * classes. These decompositions are accessed by the Matrix class to compute
 * solutions of simultaneous linear equations, determinants, inverses and other
 * matrix functions. The five decompositions are:
 * <P>
 * <UL>
 * <LI>Cholesky Decomposition of symmetric, positive definite matrices.
 * <LI>LU Decomposition of rectangular matrices.
 * <LI>QR Decomposition of rectangular matrices.
 * <LI>Singular Value Decomposition of rectangular matrices.
 * <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square
 * matrices.
 * </UL>
 * <DL>
 * <DT><B>Example of use:</B></DT>
 * <P>
 * <DD>Solve a linear system A x = b and compute the residual norm, ||b - A x||.
 * <P>
 * 
 * <PRE>
 * double[][] vals = { { 1., 2., 3 }, { 4., 5., 6. }, { 7., 8., 10. } };
 * Matrix A = new Matrix(vals);
 * Matrix b = Matrix.random(3, 1);
 * Matrix x = A.solve(b);
 * Matrix r = A.times(x).minus(b);
 * double rnorm = r.normInf();
 * </PRE>
 * 
 * </DD>
 * </DL>
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 5 August 1998
 */

public class Matrix implements Cloneable, java.io.Serializable {

	/*
	 * ------------------------ Class variables ------------------------
	 */

	/**
	 * Array for internal storage of elements.
	 * 
	 * @serial internal array storage.
	 */
	private double[][] A;

	/**
	 * Row and column dimensions.
	 * 
	 * @serial row dimension.
	 * @serial column dimension.
	 */
	private int m, n;

	/*
	 * ------------------------ Constructors ------------------------
	 */

	/**
	 * Construct an m-by-n matrix of zeros.
	 * 
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of columns.
	 */

	public Matrix(int m, int n) {
		this.m = m;
		this.n = n;
		A = new double[m][n];
	}

	/**
	 * Construct an m-by-n constant matrix.
	 * 
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of columns.
	 * @param s
	 *            Fill the matrix with this scalar value.
	 */

	public Matrix(int m, int n, double s) {
		this.m = m;
		this.n = n;
		A = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = s;
			}
		}
	}

	/**
	 * Construct a matrix from a 1-D array. The row dimension will be 1.
	 * 
	 * @param A
	 *            One-dimensional array of doubles.
	 * @exception IllegalArgumentException
	 *                All rows must have the same length
	 * @author Steven Chang
	 */

	public Matrix(double[] B) {
		m = 1;
		n = B.length;
		A = new double[m][n];
		for (int i = 0; i < n; ++i) {
			A[0][i] = B[i];
		}
	}

	/**
	 * Construct a matrix from a 2-D array.
	 * 
	 * @param A
	 *            Two-dimensional array of doubles.
	 * @exception IllegalArgumentException
	 *                All rows must have the same length
	 * @see #constructWithCopy
	 */

	public Matrix(double[][] A) {
		m = A.length;
		n = A[0].length;
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException(
						"All rows must have the same length.");
			}
		}
		this.A = A;
	}

	/**
	 * Construct a matrix quickly without checking arguments.
	 * 
	 * @param A
	 *            Two-dimensional array of doubles.
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of columns.
	 */

	public Matrix(double[][] A, int m, int n) {
		this.A = A;
		this.m = m;
		this.n = n;
	}

	/**
	 * Construct a matrix from a one-dimensional packed array
	 * 
	 * @param vals
	 *            One-dimensional array of doubles, packed by columns (ala
	 *            Fortran).
	 * @param m
	 *            Number of rows.
	 * @exception IllegalArgumentException
	 *                Array length must be a multiple of m.
	 */

	public Matrix(double vals[], int m) {
		this.m = m;
		n = (m != 0 ? vals.length / m : 0);
		if (m * n != vals.length) {
			throw new IllegalArgumentException(
					"Array length must be a multiple of m.");
		}
		A = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = vals[i + j * m];
			}
		}
	}

	/*
	 * ------------------------ Public Methods ------------------------
	 */

	/**
	 * Construct a matrix from a copy of a 2-D array.
	 * 
	 * @param A
	 *            Two-dimensional array of doubles.
	 * @exception IllegalArgumentException
	 *                All rows must have the same length
	 */

	public static Matrix constructWithCopy(double[][] A) {
		int m = A.length;
		int n = A[0].length;
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException(
						"All rows must have the same length.");
			}
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return X;
	}

	/**
	 * Make a deep copy of a matrix
	 */

	public Matrix copy() {
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return X;
	}

	/**
	 * Clone the Matrix object.
	 */

	public Object clone() {
		return this.copy();
	}

	/**
	 * Access the internal two-dimensional array.
	 * 
	 * @return Pointer to the two-dimensional array of matrix elements.
	 */

	public double[][] getArray() {
		return A;
	}

	/**
	 * Copy the internal two-dimensional array.
	 * 
	 * @return Two-dimensional array copy of matrix elements.
	 */

	public double[][] getArrayCopy() {
		double[][] C = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return C;
	}

	/**
	 * Make a one-dimensional column packed copy of the internal array.
	 * 
	 * @return Matrix elements packed in a one-dimensional array by columns.
	 */

	public double[] getColumnPackedCopy() {
		double[] vals = new double[m * n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				vals[i + j * m] = A[i][j];
			}
		}
		return vals;
	}

	/**
	 * Make a one-dimensional row packed copy of the internal array.
	 * 
	 * @return Matrix elements packed in a one-dimensional array by rows.
	 */

	public double[] getRowPackedCopy() {
		double[] vals = new double[m * n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				vals[i * n + j] = A[i][j];
			}
		}
		return vals;
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
	 * Get a single element.
	 * 
	 * @param i
	 *            Row index.
	 * @param j
	 *            Column index.
	 * @return A(i,j)
	 * @exception ArrayIndexOutOfBoundsException
	 */

	public double get(int i, int j) {
		return A[i][j];
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

	public Matrix getMatrix(int i0, int i1, int j0, int j1) {
		Matrix X = new Matrix(i1 - i0 + 1, j1 - j0 + 1);
		double[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					B[i - i0][j - j0] = A[i][j];
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/**
	 * Get a submatrix.
	 * 
	 * @param r
	 *            Array of row indices.
	 * @param c
	 *            Array of column indices.
	 * @return A(r(:),c(:))
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public Matrix getMatrix(int[] r, int[] c) {
		Matrix X = new Matrix(r.length, c.length);
		double[][] B = X.getArray();
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					B[i][j] = A[r[i]][c[j]];
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/**
	 * Get a submatrix.
	 * 
	 * @param i0
	 *            Initial row index
	 * @param i1
	 *            Final row index
	 * @param c
	 *            Array of column indices.
	 * @return A(i0:i1,c(:))
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public Matrix getMatrix(int i0, int i1, int[] c) {
		Matrix X = new Matrix(i1 - i0 + 1, c.length);
		double[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = 0; j < c.length; j++) {
					B[i - i0][j] = A[i][c[j]];
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/**
	 * Get a submatrix.
	 * 
	 * @param r
	 *            Array of row indices.
	 * @param j0
	 *            Initial column index
	 * @param j1
	 *            Final column index
	 * @return A(r(:),j0:j1)
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public Matrix getMatrix(int[] r, int j0, int j1) {
		Matrix X = new Matrix(r.length, j1 - j0 + 1);
		double[][] B = X.getArray();
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					B[i][j - j0] = A[r[i]][j];
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
	 * @param s
	 *            A(i,j).
	 * @exception ArrayIndexOutOfBoundsException
	 */

	public void set(int i, int j, double s) {
		A[i][j] = s;
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param i0
	 *            Initial row index
	 * @param i1
	 *            Final row index
	 * @param j0
	 *            Initial column index
	 * @param j1
	 *            Final column index
	 * @param X
	 *            A(i0:i1,j0:j1)
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int i0, int i1, int j0, int j1, Matrix X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					A[i][j] = X.get(i - i0, j - j0);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param r
	 *            Array of row indices.
	 * @param c
	 *            Array of column indices.
	 * @param X
	 *            A(r(:),c(:))
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int[] r, int[] c, Matrix X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					A[r[i]][c[j]] = X.get(i, j);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param r
	 *            Array of row indices.
	 * @param j0
	 *            Initial column index
	 * @param j1
	 *            Final column index
	 * @param X
	 *            A(r(:),j0:j1)
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int[] r, int j0, int j1, Matrix X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					A[r[i]][j] = X.get(i, j - j0);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Set a submatrix.
	 * 
	 * @param i0
	 *            Initial row index
	 * @param i1
	 *            Final row index
	 * @param c
	 *            Array of column indices.
	 * @param X
	 *            A(i0:i1,c(:))
	 * @exception ArrayIndexOutOfBoundsException
	 *                Submatrix indices
	 */

	public void setMatrix(int i0, int i1, int[] c, Matrix X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = 0; j < c.length; j++) {
					A[i][c[j]] = X.get(i - i0, j);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/**
	 * Matrix transpose.
	 * 
	 * @return A'
	 */

	public Matrix transpose() {
		Matrix X = new Matrix(n, m);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[j][i] = A[i][j];
			}
		}
		return X;
	}

	/**
	 * One norm
	 * 
	 * @return maximum column sum.
	 */

	public double norm1() {
		double f = 0;
		for (int j = 0; j < n; j++) {
			double s = 0;
			for (int i = 0; i < m; i++) {
				s += Math.abs(A[i][j]);
			}
			f = Math.max(f, s);
		}
		return f;
	}

	/**
	 * Two norm
	 * 
	 * @return maximum singular value.
	 */

	public double norm2() {
		return (new SingularValueDecomposition(this).norm2());
	}

	/**
	 * Infinity norm
	 * 
	 * @return maximum row sum.
	 */

	public double normInf() {
		double f = 0;
		for (int i = 0; i < m; i++) {
			double s = 0;
			for (int j = 0; j < n; j++) {
				s += Math.abs(A[i][j]);
			}
			f = Math.max(f, s);
		}
		return f;
	}

	/**
	 * Frobenius norm
	 * 
	 * @return sqrt of sum of squares of all elements.
	 */

	public double normF() {
		double f = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				f = Maths.hypot(f, A[i][j]);
			}
		}
		return f;
	}

	/**
	 * Unary minus
	 * 
	 * @return -A
	 */

	public Matrix uminus() {
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = -A[i][j];
			}
		}
		return X;
	}

	/**
	 * C = A + B
	 * 
	 * @param B
	 *            another matrix
	 * @return A + B
	 */

	public Matrix plus(Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] + B.A[i][j];
			}
		}
		return X;
	}

	/**
	 * A = A + B
	 * 
	 * @param B
	 *            another matrix
	 * @return A + B
	 */

	public Matrix plusEquals(Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] + B.A[i][j];
			}
		}
		return this;
	}

	/**
	 * C = A - B
	 * 
	 * @param B
	 *            another matrix
	 * @return A - B
	 */

	public Matrix minus(Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] - B.A[i][j];
			}
		}
		return X;
	}

	/**
	 * A = A - B
	 * 
	 * @param B
	 *            another matrix
	 * @return A - B
	 */

	public Matrix minusEquals(Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] - B.A[i][j];
			}
		}
		return this;
	}

	/**
	 * Element-by-element multiplication, C = A.*B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.*B
	 */

	public Matrix arrayTimes(Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] * B.A[i][j];
			}
		}
		return X;
	}

	/**
	 * Element-by-element multiplication in place, A = A.*B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.*B
	 */

	public Matrix arrayTimesEquals(Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] * B.A[i][j];
			}
		}
		return this;
	}

	/**
	 * Element-by-element right division, C = A./B
	 * 
	 * @param B
	 *            another matrix
	 * @return A./B
	 */

	public Matrix arrayRightDivide(Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] / B.A[i][j];
			}
		}
		return X;
	}

	/**
	 * Element-by-element right division in place, A = A./B
	 * 
	 * @param B
	 *            another matrix
	 * @return A./B
	 */

	public Matrix arrayRightDivideEquals(Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] / B.A[i][j];
			}
		}
		return this;
	}

	/**
	 * Element-by-element left division, C = A.\B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.\B
	 */

	public Matrix arrayLeftDivide(Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = B.A[i][j] / A[i][j];
			}
		}
		return X;
	}

	/**
	 * Element-by-element left division in place, A = A.\B
	 * 
	 * @param B
	 *            another matrix
	 * @return A.\B
	 */

	public Matrix arrayLeftDivideEquals(Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = B.A[i][j] / A[i][j];
			}
		}
		return this;
	}

	/**
	 * Multiply a matrix by a scalar, C = s*A
	 * 
	 * @param s
	 *            scalar
	 * @return s*A
	 */

	public Matrix times(double s) {
		Matrix X = new Matrix(m, n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = s * A[i][j];
			}
		}
		return X;
	}

	/**
	 * Multiply a matrix by a scalar in place, A = s*A
	 * 
	 * @param s
	 *            scalar
	 * @return replace A by s*A
	 */

	public Matrix timesEquals(double s) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = s * A[i][j];
			}
		}
		return this;
	}

	/**
	 * Linear algebraic matrix multiplication, A * B
	 * 
	 * @param B
	 *            another matrix
	 * @return Matrix product, A * B
	 * @exception IllegalArgumentException
	 *                Matrix inner dimensions must agree.
	 */

	public Matrix times(Matrix B) {
		if (B.m != n) {
			throw new IllegalArgumentException(
					"Matrix inner dimensions must agree.");
		}
		Matrix X = new Matrix(m, B.n);
		double[][] C = X.getArray();
		double[] Bcolj = new double[n];
		for (int j = 0; j < B.n; j++) {
			for (int k = 0; k < n; k++) {
				Bcolj[k] = B.A[k][j];
			}
			for (int i = 0; i < m; i++) {
				double[] Arowi = A[i];
				double s = 0;
				for (int k = 0; k < n; k++) {
					s += Arowi[k] * Bcolj[k];
				}
				C[i][j] = s;
			}
		}
		return X;
	}

	/**
	 * LU Decomposition
	 * 
	 * @return LUDecomposition
	 * @see LUDecomposition
	 */

	public LUDecomposition lu() {
		return new LUDecomposition(this);
	}

	/**
	 * QR Decomposition
	 * 
	 * @return QRDecomposition
	 * @see QRDecomposition
	 */

	public QRDecomposition qr() {
		return new QRDecomposition(this);
	}

	/**
	 * Cholesky Decomposition
	 * 
	 * @return CholeskyDecomposition
	 * @see CholeskyDecomposition
	 */

	public CholeskyDecomposition chol() {
		return new CholeskyDecomposition(this);
	}

	/**
	 * Singular Value Decomposition
	 * 
	 * @return SingularValueDecomposition
	 * @see SingularValueDecomposition
	 */

	public SingularValueDecomposition svd() {
		return new SingularValueDecomposition(this);
	}

	/**
	 * Eigenvalue Decomposition
	 * 
	 * @return EigenvalueDecomposition
	 * @see EigenvalueDecomposition
	 */

	public EigenvalueDecomposition eig() {
		return new EigenvalueDecomposition(this);
	}

	/**
	 * Solve A*X = B
	 * 
	 * @param B
	 *            right hand side
	 * @return solution if A is square, least squares solution otherwise
	 */

	public Matrix solve(Matrix B) {
		return (m == n ? (new LUDecomposition(this)).solve(B)
				: (new QRDecomposition(this)).solve(B));
	}

	/**
	 * Solve X*A = B, which is also A'*X' = B'
	 * 
	 * @param B
	 *            right hand side
	 * @return solution if A is square, least squares solution otherwise.
	 */

	public Matrix solveTranspose(Matrix B) {
		return transpose().solve(B.transpose());
	}

	/**
	 * Matrix inverse or pseudoinverse
	 * 
	 * @return inverse(A) if A is square, pseudoinverse otherwise.
	 */

	public Matrix inverse() {
		return solve(identity(m, m));
	}

	/**
	 * Matrix determinant
	 * 
	 * @return determinant
	 */

	public double det() {
		return new LUDecomposition(this).det();
	}

	/**
	 * Matrix rank
	 * 
	 * @return effective numerical rank, obtained from SVD.
	 */

	public int rank() {
		return new SingularValueDecomposition(this).rank();
	}

	/**
	 * Matrix condition (2 norm)
	 * 
	 * @return ratio of largest to smallest singular value.
	 */

	public double cond() {
		return new SingularValueDecomposition(this).cond();
	}

	/**
	 * Matrix trace.
	 * 
	 * @return sum of the diagonal elements.
	 */

	public double trace() {
		double t = 0;
		for (int i = 0; i < Math.min(m, n); i++) {
			t += A[i][i];
		}
		return t;
	}

	/**
	 * Generate matrix with random elements
	 * 
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of colums.
	 * @return An m-by-n matrix with uniformly distributed random elements.
	 */

	public static Matrix random(int m, int n) {
		Matrix A = new Matrix(m, n);
		double[][] X = A.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				X[i][j] = Math.random();
			}
		}
		return A;
	}

	/**
	 * Generate identity matrix
	 * 
	 * @param m
	 *            Number of rows.
	 * @param n
	 *            Number of colums.
	 * @return An m-by-n matrix with ones on the diagonal and zeros elsewhere.
	 */

	public static Matrix identity(int m, int n) {
		Matrix A = new Matrix(m, n);
		double[][] X = A.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				X[i][j] = (i == j ? 1.0 : 0.0);
			}
		}
		return A;
	}

	/**
	 * Print the matrix to stdout. Line the elements up in columns with a
	 * Fortran-like 'Fw.d' style format.
	 * 
	 * @param w
	 *            Column width.
	 * @param d
	 *            Number of digits after the decimal.
	 */

	public void print(int w, int d) {
		print(new PrintWriter(System.out, true), w, d);
	}

	/**
	 * Print the matrix to the output stream. Line the elements up in columns
	 * with a Fortran-like 'Fw.d' style format.
	 * 
	 * @param output
	 *            Output stream.
	 * @param w
	 *            Column width.
	 * @param d
	 *            Number of digits after the decimal.
	 */

	public void print(PrintWriter output, int w, int d) {
		DecimalFormat format = new DecimalFormat();
		format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
		format.setMinimumIntegerDigits(1);
		format.setMaximumFractionDigits(d);
		format.setMinimumFractionDigits(d);
		format.setGroupingUsed(false);
		print(output, format, w + 2);
	}

	/**
	 * Print the matrix to stdout. Line the elements up in columns. Use the
	 * format object, and right justify within columns of width characters. Note
	 * that is the matrix is to be read back in, you probably will want to use a
	 * NumberFormat that is set to US Locale.
	 * 
	 * @param format
	 *            A Formatting object for individual elements.
	 * @param width
	 *            Field width for each column.
	 * @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */

	public void print(NumberFormat format, int width) {
		print(new PrintWriter(System.out, true), format, width);
	}

	// DecimalFormat is a little disappointing coming from Fortran or C's
	// printf.
	// Since it doesn't pad on the left, the elements will come out different
	// widths. Consequently, we'll pass the desired column width in as an
	// argument and do the extra padding ourselves.

	/**
	 * Print the matrix to the output stream. Line the elements up in columns.
	 * Use the format object, and right justify within columns of width
	 * characters. Note that is the matrix is to be read back in, you probably
	 * will want to use a NumberFormat that is set to US Locale.
	 * 
	 * @param output
	 *            the output stream.
	 * @param format
	 *            A formatting object to format the matrix elements
	 * @param width
	 *            Column width.
	 * @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */

	public void print(PrintWriter output, NumberFormat format, int width) {
		output.println(); // start on new line.
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				String s = format.format(A[i][j]); // format the number
				int padding = Math.max(1, width - s.length()); // At _least_ 1
																// space
				for (int k = 0; k < padding; k++)
					output.print(' ');
				output.print(s);
			}
			output.println();
		}
		output.println(); // end with blank line.
	}

	/**
	 * Read a matrix from a stream. The format is the same the print method, so
	 * printed matrices can be read back in (provided they were printed using US
	 * Locale). Elements are separated by whitespace, all the elements for each
	 * row appear on a single line, the last row is followed by a blank line.
	 * 
	 * @param input
	 *            the input stream.
	 */

	public static Matrix read(BufferedReader input) throws java.io.IOException {
		StreamTokenizer tokenizer = new StreamTokenizer(input);

		// Although StreamTokenizer will parse numbers, it doesn't recognize
		// scientific notation (E or D); however, Double.valueOf does.
		// The strategy here is to disable StreamTokenizer's number parsing.
		// We'll only get whitespace delimited words, EOL's and EOF's.
		// These words should all be numbers, for Double.valueOf to parse.

		tokenizer.resetSyntax();
		tokenizer.wordChars(0, 255);
		tokenizer.whitespaceChars(0, ' ');
		tokenizer.eolIsSignificant(true);
		java.util.Vector<Double> vD = new java.util.Vector<Double>();

		// Ignore initial empty lines
		while (tokenizer.nextToken() == StreamTokenizer.TT_EOL)
			;
		if (tokenizer.ttype == StreamTokenizer.TT_EOF)
			throw new java.io.IOException("Unexpected EOF on matrix read.");
		do {
			vD.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st
															// row.
		} while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

		int n = vD.size(); // Now we've got the number of columns!
		double row[] = new double[n];
		for (int j = 0; j < n; j++)
			// extract the elements of the 1st row.
			row[j] = vD.elementAt(j).doubleValue();
		java.util.Vector<double[]> v = new java.util.Vector<double[]>();
		v.addElement(row); // Start storing rows instead of columns.
		while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
			// While non-empty lines
			v.addElement(row = new double[n]);
			int j = 0;
			do {
				if (j >= n)
					throw new java.io.IOException("Row " + v.size()
							+ " is too long.");
				row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
			} while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
			if (j < n)
				throw new java.io.IOException("Row " + v.size()
						+ " is too short.");
		}
		int m = v.size(); // Now we've got the number of rows.
		double[][] A = new double[m][];
		v.copyInto(A); // copy the rows out of the vector
		return new Matrix(A);
	}

	/**
	 * Return the max value of specific dimension.Read a matrix and specify the
	 * dimension. dim can only be 1 or 2.
	 * 
	 * @param dim
	 *            dimension
	 * @author Steven Chang
	 */
	public Matrix max(int dim) throws Exception {

		// MatLab max simplified function
		if (dim == 1) {
			// row dimension
			Matrix X = new Matrix(1, n);
			for (int j = 0; j < n; ++j) {
				double max = A[0][j];
				for (int i = 0; i < m; ++i) {
					if (max < A[i][j]) {
						max = A[i][j];
					}
				}
				X.set(0, j, max);
			}
			return X;
		} else if (dim == 2) {
			// col dimension
			Matrix X = new Matrix(m, 1);
			for (int i = 0; i < m; ++i) {
				double max = A[i][0];
				for (int j = 0; j < n; ++j) {
					if (max < A[i][j]) {
						max = A[i][j];
					}
				}
				X.set(i, 0, max);
			}
			return X;
		} else {
			throw new IllegalArgumentException("row and col should be the same");
		}
	}

	/**
	 * Return the max value of matrix data.
	 * 
	 * @author Steven Chang
	 */
	public double max() {

		double max = A[0][0];
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				if (max < A[i][j]) {
					max = A[i][j];
				}
			}
		}
		return max;
	}

	/**
	 * Return the min value of specifiec dimension.Read a matrix and specify the
	 * dimension. dim can only be 1 or 2.
	 * 
	 * @param dim
	 *            dimension
	 * @author Steven Chang
	 */
	public Matrix min(int dim) throws Exception {

		// MatLab max simplified funtion
		if (dim == 1) {
			// row dimension
			Matrix X = new Matrix(1, n);
			for (int j = 0; j < n; ++j) {
				double min = A[0][j];
				for (int i = 0; i < m; ++i) {
					if (min > A[i][j]) {
						min = A[i][j];
					}
				}
				X.set(0, j, min);
			}
			return X;
		} else if (dim == 2) {
			// col dimension
			Matrix X = new Matrix(m, 1);
			for (int i = 0; i < m; ++i) {
				double min = A[i][0];
				for (int j = 0; j < n; ++j) {
					if (min > A[i][j]) {
						min = A[i][j];
					}
				}
				X.set(i, 0, min);
			}
			return X;
		} else {
			throw new IllegalArgumentException("row and col should be the same");
		}
	}

	/**
	 * Return the min value of matrix data.
	 * 
	 * @author Steven Chang
	 */
	public double min() {

		double min = A[0][0];
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				if (min > A[i][j]) {
					min = A[i][j];
				}
			}
		}
		return min;
	}

	/**
	 * Return the mean value of specifiec dimension.Read a matrix and specify
	 * the dimension. dim can only be 1 or 2.
	 * 
	 * @param dim
	 *            dimension
	 * @author Steven Chang
	 */

	public Matrix mean(int dim) throws Exception {

		// MatLab mean simplified funtion
		if (dim == 1) {
			// row dimension
			Matrix X = new Matrix(1, n);
			for (int j = 0; j < n; ++j) {
				double means = 0;
				for (int i = 0; i < m; ++i) {
					means += A[i][j];
				}
				X.set(0, j, means / (float) m);
			}
			return X;
		} else if (dim == 2) {
			// col dimension
			Matrix X = new Matrix(m, 1);
			for (int i = 0; i < m; ++i) {
				double means = 0;
				for (int j = 0; j < n; ++j) {
					means += A[i][j];
				}
				X.set(i, 0, means / (float) n);
			}
			return X;
		} else {
			throw new IllegalArgumentException("row and col should be the same");
		}
	}
	
	/**
	 * Return the mean value of the matrix.
	 * 
	 * @author Steven Chang
	 */

	public double mean() throws Exception {

		double mean = 0.;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				mean += A[i][j]/(double)m/(double)n;
			}
		}
		return mean;
	}

	/**
	 * Realize the simplified version of repmat like MatLab. Duplicate the value
	 * to specified size. The size of returned matrix is m*r1, n*r2. Throw
	 * exception if r1 or r2 less than 1. Added by Steven Chang
	 * 
	 * @param r1
	 *            multiplier of the first dimension size
	 * @param r2
	 *            multiplier of the second dimension size
	 * @author Steven Chang
	 */

	public Matrix repmat(int r1, int r2) throws Exception {

		// MatLab repmat simplified function
		if (r1 < 1 || r2 < 1) {
			throw new IllegalArgumentException(
					"r1 and r2 should be larger than 0");
		}
		int nrow = m * r1;
		int ncol = n * r2;
		Matrix X = new Matrix(nrow, ncol);
		for (int i = 0; i < nrow; ++i) {
			for (int j = 0; j < ncol; ++j) {
				X.set(i, j, A[i % m][j % n]);
			}
		}
		return X;
	}

	/**
	 * Realize the simplified version of reshape like MatLab. It will reshape
	 * the matrix with given size. If the size is smaller than before, data at
	 * the end will be lost. If the size is larger, duplicate the former data
	 * and add 0 to fill the matrix. Throw exception if nrow or ncol less than
	 * 1.
	 * 
	 * @param nrow
	 *            new row number of reshaped matrix
	 * @param ncol
	 *            new col number of reshaped matrix
	 * @author Steven Chang
	 */
	public Matrix reshape(int nrow, int ncol) throws Exception {
		return reshape(nrow, ncol, 0.0);
	}

	/**
	 * Realize the simplified version of reshape like MatLab. It will reshape
	 * the matrix with given size. If the size is smaller than before, data at
	 * the end will be lost. If the size is larger, duplicate the former data
	 * and add 0 to fill the matrix. Throw exception if nrow or ncol less than
	 * 1.
	 * 
	 * @param nrow
	 *            new row number of reshaped matrix
	 * @param ncol
	 *            new col number of reshaped matrix
	 * @param fit
	 *            fit is used to initialize those elements without covered.
	 * @author Steven Chang
	 */
	public Matrix reshape(int nrow, int ncol, double fit) throws Exception {

		// Simplified version of reshape method
		if (nrow < 1 || ncol < 1) {
			throw new IllegalArgumentException(
					"nrow and ncol should be larger than 0");
		}
		int tprow = 0;
		int tpcol = 0;
		if (nrow * ncol >= m * n) {
			// larger
			Matrix X = new Matrix(nrow, ncol, fit);
			double[][] array = X.getArray();
			for (int j = 0; j < n; ++j) {
				for (int i = 0; i < m; ++i, ++tprow) {
					if (tprow == nrow) {
						++tpcol;
						tprow = 0;
					}
					array[tprow][tpcol] = A[i][j];
				}
			}
			return X;
		} else {
			// smaller
			Matrix X = new Matrix(nrow, ncol, fit);
			double[][] array = X.getArray();
			for (int j = 0; j < ncol; ++j) {
				for (int i = 0; i < nrow; ++i, ++tprow) {
					if (tprow == m) {
						++tpcol;
						tprow = 0;
					}
					array[i][j] = A[tprow][tpcol];
				}
			}
			return X;
		}
	}

	/**
	 * Perform two matrix concatenation. There are two kinds of cat: if dim==1
	 * then equals to MatLab code [src1, src2] -> cat by row if dim==2 then
	 * equals to MatLab code [src1; src2] -> cat by col For specified
	 * concatenaiton, num of rows or cols must be equal, or there will be an
	 * exception. Also, if dim is not 1 or 2 throw exception.
	 * 
	 * @param src1
	 *            matrix1
	 * @param src2
	 *            matrix2
	 * @param dim
	 *            dimension
	 * @author Steven Chang
	 */
	public static Matrix concatenate(Matrix src1, Matrix src2, int dim)
			throws Exception {

		if (dim == 1) {
			// [src1, src2]
			if (src1.getColumnDimension() != src2.getColumnDimension()) {
				throw new IllegalArgumentException(
						"Number of col must be equal");
			}
			int row1 = src1.getRowDimension();
			int row2 = src2.getRowDimension();
			int col = src1.getColumnDimension();
			Matrix X = new Matrix(row1 + row2, col);
			X.setMatrix(0, row1 - 1, 0, col - 1, src1);
			X.setMatrix(row1, row1 + row2 - 1, 0, col - 1, src2);
			return X;
		} else if (dim == 2) {
			// [src1;src2]
			if (src1.getRowDimension() != src2.getRowDimension()) {
				throw new IllegalArgumentException(
						"Number of rows must be equal");
			}
			int row = src1.getRowDimension();
			int col1 = src1.getColumnDimension();
			int col2 = src2.getColumnDimension();
			Matrix X = new Matrix(row, col1 + col2);
			X.setMatrix(0, row - 1, 0, col1 - 1, src1);
			X.setMatrix(0, row - 1, col1, col1 + col2 - 1, src2);
			return X;
		} else {
			throw new IllegalArgumentException("row and col should be the same");
		}
	}

	/**
	 * Return the size of elements. m*n.
	 * 
	 * @author Steven Chang
	 */
	public int elementSize() {
		return m * n;
	}

	/**
	 * Realize a simplified version of pdist funciton in MatLab. Can only
	 * perform Euclidean dist. Dimension one will be the member vector as the
	 * second dimension will be data vector of each element.
	 * 
	 * @return Matrix, a [m(m-1)/2, 1] row vector
	 * @author Steven Chang
	 */
	public Matrix pdist() {

		// Simplified MatLab version - Euclidean
		int size = m * n;
		int tpIndex = 0;
		Matrix X = new Matrix(size * (size - 1) / 2, 1);
		for (int i = 0; i < m; ++i) {
			for (int j = i + 1; j < m; ++j, ++tpIndex) {
				double dist = 0.0;
				for (int k = 0; k < n; ++k) {
					Maths.hypot(dist, (A[i][k] - A[j][k]));
				}
				X.set(tpIndex, 0, dist);
			}
		}
		return X;
	}

	/**
	 * Reverse the matrix but not change its size.
	 * 
	 * @author Steven Chang
	 */
	public Matrix reverse() {

		Matrix X = new Matrix(m, n);
		for (int index = 0; index < this.elementSize(); ++index) {
			X.set(this.elementSize() - index - 1, this.get(index));
		}
		return X;
	}

	/**
	 * Reverse the matrix but not create a new matrix.
	 * 
	 * @author Steven Chang
	 */
	public void reverseEqual() {

		for (int index = 0; index < this.elementSize() / 2; ++index) {
			double value = this.get(index);
			int rindex = this.elementSize() - index - 1;
			this.set(index, this.get(rindex));
			this.set(rindex, value);
		}
	}

	/**
	 * Simplified implementation of MatLab command squareform. Convert a
	 * distancvector into a square matrix. Data itself should be a vector. Or it
	 * will throw an exception. Data itself should be a [1, m*(m-1)/2] or
	 * [m*(m-1)/2, 1] form.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public Matrix squareform() throws Exception {

		if (m != 1 && n != 1) {
			throw new IllegalArgumentException("Matrix should be a vector.");
		}
		int tpVal = (int) Math.sqrt(2 * m * n);
		if (tpVal * (tpVal + 1) != 2 * m * n) {
			throw new IllegalArgumentException(
					"Matrix should be a dist vector.");
		}
		++tpVal;
		int tpIndex = 0;
		Matrix X = new Matrix(tpVal, tpVal, 0.0);
		if (m == 1) {
			for (int i = 0; i < tpVal; ++i) {
				for (int j = i + 1; j < tpVal; ++j, ++tpIndex) {
					X.set(i, j, A[0][tpIndex]);
					X.set(j, i, A[0][tpIndex]);
				}
			}
		} else {
			for (int i = 0; i < tpVal; ++i) {
				for (int j = i + 1; j < tpVal; ++j, ++tpIndex) {
					X.set(i, j, A[tpIndex][0]);
					X.set(j, i, A[tpIndex][0]);
				}
			}
		}
		return X;
	}

	/**
	 * Get one row.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public Matrix getRow(int index) throws Exception {
		return this.getMatrix(index, index, 0, n-1);
	}

	/**
	 * Get one col.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public Matrix getCol(int index) throws Exception {
		return this.getMatrix(0, m-1, index, index);
	}

	/**
	 * Return a 1 x n matrix containing the number of indices of element
	 * matching the value.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public Matrix find(double value) throws Exception {

		int count = this.find_number(value);
		// Assign the index
		Matrix X = new Matrix(1, count);
		double[][] Xarr = X.getArray();
		int tpCount = 0;
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				if (A[i][j] == value) {
					Xarr[0][tpCount] = i + j * m;
					++tpCount;
				}
			}
		}
		return X;
	}

	/**
	 * Find the first element matching the value, and return its index.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public int find_first(double value) throws Exception {
		int count = this.elementSize();
		for (int i = 0; i < count; ++i) {
			if (this.get(i) == value) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Find the first row matching the value, and return its row index.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public int find_first_row(Matrix row) throws Exception {

		if (row.getRowDimension() != 1) {
			throw new IllegalArgumentException("argument must be a row vector");
		} else if (row.getColumnDimension() != n) {
			throw new IllegalArgumentException(
					"the argument must have the same columns number");
		}
		for (int i = 0; i < m; ++i) {
			int j = 0;
			for (j = 0; j < n; ++j) {
				if (A[i][j] != row.get(j))
					break;
			}
			if (j == n) {
				// Find
				return i;
			}
		}
		return -1;
	}

	/**
	 * Compare with another matrix if every element is the same.
	 * 
	 * @param mat
	 *            matrix to compare
	 * @author Steven Chang
	 * @throws Exception
	 */
	public boolean equals(Matrix mat) throws Exception {

		checkMatrixDimensions(mat);
		double[][] B = mat.getArray();
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				if (A[i][j] - B[i][j] < 0.e-7) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Check each element if equals to a double value. return a matrix with 1
	 * and 0.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public Matrix equals(double value) throws Exception {

		Matrix X = new Matrix(m, n, 0.0);
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				if (A[i][j] == value) {
					X.set(i, j, 1.0);
				}
			}
		}
		return X;
	}

	/**
	 * Check each element if equals to a double value. change the value to
	 * substitute.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public Matrix equalsSubstitute(double value, double substitute)
			throws Exception {

		Matrix X = new Matrix(m, n, 0.0);
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				if (A[i][j] == value) {
					X.set(i, j, substitute);
				}
			}
		}
		return X;
	}

	/**
	 * Find the number of elements matching the value.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */
	public int find_number(double value) throws Exception {

		int count = 0;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				if (A[i][j] == value) {
					++count;
				}
			}
		}
		return count;
	}

	/**
	 * Get a single element from the index. Search order: m-n
	 * 
	 * @param index
	 *            Row index. Start from 0.
	 * @return A(i,j)
	 * @author Steven Chang
	 */

	public double get(int index) {
		int row = index % m;
		int col = index / m;
		return A[row][col];
	}

	/**
	 * Set a single element from the index. Search order: m-n
	 * 
	 * @param index
	 *            Row index. Start from 0.
	 * @param val
	 *            Value
	 * @author Steven Chang
	 */

	public void set(int index, double val) {
		int row = index % m;
		int col = index / m;
		A[row][col] = val;
	}

	/**
	 * Sort all elements in the matrix.
	 * 
	 * @author Steven Chang
	 * @throws Exception
	 */

	public Matrix sort() throws Exception {

		return sort(1);
	}

	/**
	 * Sort the matrix by specified dimension.
	 * 
	 * @param dim
	 *            Specified dimension to perform sort method
	 * @author Steven Chang
	 */

	public Matrix sort(int dim) {

		Matrix X = Matrix.constructWithCopy(A);
		if (dim == 1) {
			// ok
		} else if (dim == 2) {
			// Perform transposition
			X = X.transpose();
		} else {
			throw new IllegalArgumentException("dim must be 1 or 2");
		}
		// For each row, perform quick sort
		for (int i = 0; i < m; ++i) {
			quicksort(0, n - 1, A[i]);
		}
		if (dim == 2) {
			X = X.transpose();
		}
		return X;
	}

	/**
	 * Construct a matrix containing elements from index istart to iend(included)
	 * 
	 * @param istart
	 *            Start of the copied array index.
	 * @param iend
	 *            End of the copied array index. Not included.
	 * @return X A matrix with size (iend - istart) * 1
	 * @author Steven Chang
	 */

	public Matrix getMatrix(int istart, int iend) throws Exception {

		if (iend<=istart) {throw new IllegalArgumentException();}
		Matrix X = new Matrix((iend - istart)+1, 1);
		double[][] array = X.getArray();
		int index = 0;
		for (int i = istart; i < iend+1; ++i, ++index) {
			array[index][0] = this.get(i);
		}
		return X;
	}

	/**
	 * Construct a matrix containing elements mat pointed at in A.
	 * 
	 * @param mat
	 *            storing indices of data in A.
	 * @author Steven Chang
	 */

	public Matrix getMatrix(Matrix mat) throws Exception {

		Matrix X = Matrix.constructWithCopy(mat.getArray());
		double[][] array = X.getArray();
		for (int i = 0; i < mat.getRowDimension(); ++i) {
			for (int j = 0; j < mat.getColumnDimension(); ++j) {
				array[i][j] = this.get((int)mat.get(i, j));
			}
		}
		return X;
	}

	/**
	 * Get a matrix containing rows that mat pointed at in A.
	 * 
	 * @param dim
	 *            row index storage matrix
	 * @author Steven Chang
	 */

	public Matrix getRows(Matrix mat) throws Exception {
		
		Matrix X = new Matrix(mat.elementSize(), n);
		double[][] array = X.getArray();
		for (int k=0; k<mat.elementSize(); ++k) {
			for(int j=0; j<n; ++j) {
				array[k][j] = A[(int)mat.get(k)][j];
			}
		}
		return X;
	}

	/**
	 * GGet a matrix containing rows that mat pointed at in A.
	 * 
	 * @param mat
	 *            col index storage matrix
	 * @author Steven Chang
	 */

	public Matrix getCols(Matrix mat) throws Exception {

		Matrix X = new Matrix(m, mat.elementSize());
		double[][] array = X.getArray();
		for (int k=0; k<mat.elementSize(); ++k) {
			for(int i=0; i<m; ++i) {
				array[i][k] = A[i][(int)mat.get(k)];
			}
		}
		return X;
	}

	/**
	 * Perform abs() method for all elements in the matrix.
	 * 
	 * @author Steven Chang
	 */

	public Matrix abs() throws Exception {

		double[][] B = new double[m][n];
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < m; ++j) {
				if (A[i][j] < 0)
					B[i][j] = -A[i][j];
				else
					B[i][j] = A[i][j];
			}
		}
		return new Matrix(B);
	}

	/**
	 * get the sum of all elements.
	 * 
	 * @author Steven Chang
	 */

	public double sum() throws Exception {

		double val = 0.0;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				val += A[i][j];
			}
		}
		return val / this.elementSize();
	}

	/**
	 * get the sum of assigned dimensions.
	 * 
	 * @param dim
	 *            dim=1: row dim=2: col
	 * 
	 * @author Steven Chang
	 */

	public Matrix sum(int dim) throws Exception {

		if (dim == 1) {
			// row
			Matrix X = new Matrix(m, 1);
			double[][] arr = X.getArray();
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					arr[i][0] += A[i][j] / n;
				}
			}
			return X;
		} else if (dim == 2) {
			// col
			Matrix X = new Matrix(1, n);
			double[][] arr = X.getArray();
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					arr[0][j] += A[i][j] / m;
				}
			}
			return X;
		} else {
			throw new IllegalArgumentException("dim must be 1 or 2");
		}
	}

	/**
	 * Fill the elements given by arguments. It will create an array from start
	 * to end(including end) with the size of m*n, and fill each elements into
	 * A. The insertion order will from row to column.
	 * 
	 * @param start
	 *            start of the array
	 * @param end
	 *            end of the array
	 * @author Steven Chang
	 */

	public void fill(double start, double end) throws Exception {

		double stride = (end - start) / (this.elementSize() - 1);
		double val = start;
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				A[i][j] = val;
				val += stride;
			}
		}
	}

	/**
	 * Simplified version of setdiff method in MatLab. Find
	 * 
	 * @param mat
	 *            Another matrix to compare.
	 * @return X Index of different elements in both matrices.
	 * @author Steven Chang
	 */

	public Matrix setdiff(Matrix mat) throws Exception {

		checkMatrixDimensions(mat);
		double[] arr = new double[m * n];
		int numOfDiff = 0;
		for (int i = 0; i < m * n; ++i) {
			if (this.get(i) != mat.get(i)) {
				arr[numOfDiff] = i;
				++numOfDiff;
			}
		}
		// Copy and construct
		Matrix X = new Matrix(numOfDiff, 1);
		for (int i = 0; i < numOfDiff; ++i) {
			X.set(i, 0, arr[i]);
		}
		return X;
	}

	/**
	 * Create a matrix, containing indices of elements from this pointed at in
	 * mat by specified dimension. Size of two matrices must match. Also for
	 * each row/column, the set of elements must be the same. Similar to this(X)
	 * = mat in MatLab.
	 * 
	 * @param mat
	 *            Another matrix to compare.
	 * @param dim
	 *            Specified dimension
	 * @return X this->mat
	 * @author Steven Chang
	 */

	public Matrix buildBind(Matrix mat, int dim) throws Exception {

		checkMatrixDimensions(mat);
		Matrix X = this.copy();
		double[][] arr = X.getArray();
		if (dim == 1) {
			// row dimension
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					int k = 0;
					for (; k < n; ++k) {
						if (A[i][j] == mat.get(i, k))
							break;
					}
					if (k == n) {
						throw new IllegalArgumentException(
								"Data set of each row must be the same");
					}
					arr[i][k] = j;
				}
			}
		} else if (dim == 2) {
			// column dimension
			for (int j = 0; j < n; ++j) {
				for (int i = 0; i < m; ++i) {
					int k = 0;
					for (; k < m; ++k) {
						if (A[i][j] == mat.get(k, j))
							break;
					}
					if (k == n) {
						throw new IllegalArgumentException(
								"Data set of each row must be the same");
					}
					arr[k][j] = i;
				}
			}
		} else {
			throw new IllegalArgumentException("dim must be 1 or 2");
		}
		return X;
	}

	/*
	 * ------------------------ Private Methods ------------------------
	 */

	/** Quick sort method **/
	private void quicksort(int s, int t, double[] a2) {
		int i = s, j = t;
		double x = a2[(int) (i + j) / 2], y;
		do {
			while (a2[i] < x)
				i++;
			while (a2[i] > x)
				j--;
			if (i <= j) {
				y = a2[i];
				a2[j] = a2[i];
				a2[i] = y;
				i++;
				j--;
			}
		} while (i < j);
		if (j > s) {
			quicksort(s, j, a2);
		}
		if (i < t) {
			quicksort(i, t, a2);
		}
	}

	/** Check if size(A) == size(B) **/

	private void checkMatrixDimensions(Matrix B) {
		if (B.m != m || B.n != n) {
			throw new IllegalArgumentException("Matrix dimensions must agree.");
		}
	}

	private static final long serialVersionUID = 1;
}
