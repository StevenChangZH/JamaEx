# JamaEx
Some extension for Jama (Java Matrix Library). You can find source code [here](http://math.nist.gov/javanumerics/jama/)

Extension for Matrix class:

Add methods:

```
public static Matrix concatenate(Matrix src1, Matrix src2, int dim) throws Exception;
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
public int size();
public Matrix squareform() throws Exception;
```

Not tested yet. Later I'll add test codes.