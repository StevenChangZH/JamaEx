# JamaEx
Some extension for Jama (Java Matrix Library). You can find source code [here](http://math.nist.gov/javanumerics/jama/)

Extension for Matrix class:

Add methods:

```
public Matrix max(int dim) throws Exception;
public Matrix min(int dim) throws Exception;
public Matrix mean(int dim) throws Exception;
public Matrix repmat(int r1, int r2) throws Exception;
public Matrix reshape(int nrow, int ncol) throws Exception;
public static Matrix concatenate(Matrix src1, Matrix src2, int dim) throws Exception;
public int size();
```

Not tested yet. Later I'll add test codes.