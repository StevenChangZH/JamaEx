# JamaEx
Some extension for Jama (Java Matrix Class). You can find source code [here](http://math.nist.gov/javanumerics/jama/).

Notice the order of index in Matrix class follows as how MatLab does. That is, mat(1) is the point A(1, 0) but not A(0, 1). In fact, you may find Jama was built former by people in Mathworks, so in some way it is a MatLab in Java.

My advice: use it as MatLab. Do not raise high expectations on performance.

## Extensions
Extension for Matrix class:

Add methods:

```
public Matrix(double[] B);

public Matrix abs() throws Exception;
public Matrix buildBind(Matrix mat, int dim) throws Exception;
public static Matrix concatenate(Matrix src1, Matrix src2, int dim) throws Exception;
public Matrix equals(double value) throws Exception;
public boolean equals(Matrix mat) throws Exception;
public Matrix equalsSustitute(double value, double substitute) throws Exception;
public void fill(double start, double end) throws Exception;
public Matrix find(double value) throws Exception;
public int find_first(double value) throws Exception;
public int find_number(double value) throws Exception;
public double get(int index);
public Matrix getCol(int index) throws Exception;
public Matrix getCols(Matrix mat) throws Exception;
public Matrix getMatrix(int istart, int iend) throws Exception;
public Matrix getMatrix(Matrix mat) throws Exception;
public Matrix getRow(int index) throws Exception;
public Matrix getRows(Matrix mat) throws Exception;
public Matrix max(int dim) throws Exception;
public Matrix max();
public Matrix mean(int dim) throws Exception;
public Matrix min(int dim) throws Exception;
public Matrix min();
public Matrix repmat(int r1, int r2) throws Exception;
public Matrix reshape(int nrow, int ncol) throws Exception;
public Matrix reverse();
public void reverseEqual();
public Matrix pdist(Matrix matrix);
public void set(int index, double val);
public Matrix setdiff(Matrix mat) throws Exception;
public int size();
public Matrix sort() throws Exception;
public Matrix sort(int dim);
public Matrix squareform() throws Exception;
public double sum() throws Exception;
public Matrix sum(int dim) throws Exception;

private void quicksort(int s, int t, double[] a2);
```

Add another class called Matrix3. More methods will be updated later.

Not tested yet. Later I'll add test codes.