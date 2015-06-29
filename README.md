# JamaEx
Some extension for Jama (Java Matrix Class). You can find source code [here](http://math.nist.gov/javanumerics/jama/).

Notice the order of index in Matrix class follows as how MatLab does. That is, mat(1) is the point A(1, 0) but not A(0, 1). In fact, you may find Jama was built former by people in Mathworks, so in some way it is a MatLab in Java.

My advice: use it as MatLab. Do not raise high expectations on performance.

## Extensions
Extensions for Matrix class:

Add methods:

```
Matrix(double[] B);

Matrix abs();
Matrix buildBind(Matrix mat, int dim);
static Matrix concatenate(Matrix src1, Matrix src2, int dim);
int elementSize();
Matrix equals(double value);
boolean equals(Matrix mat);
Matrix equalsSustitute(double value, double substitute);
void fill(double start, double end);
Matrix find(double value);
int find_first(double value);
int find_first_row(Matrix row);
int find_number(double value);
double get(int index);
Matrix getCol(int index);
Matrix getCols(Matrix mat);
Matrix getMatrix(int istart, int iend);
Matrix getMatrix(Matrix mat);
Matrix getRow(int index);
Matrix getRows(Matrix mat);
Matrix max(int dim);
double max();
Matrix mean(int dim);
double mean();
Matrix min(int dim);
double min();
Matrix repmat(int r1, int r2);
Matrix reshape(int nrow, int ncol);
Matrix reshape(int nrow, int ncol, double fit);
Matrix reverse();
void reverseEqual();
Matrix pdist(Matrix matrix);
void set(int index, double val);
Matrix setdiff(Matrix mat);
Matrix sort();
Matrix sort(int dim);
Matrix squareform();
double sum();
Matrix sum(int dim);

private void quicksort(int s, int t, double[] a2);
```

Add another class called Matrix3. More methods will be updated later.

Tests do not cover all for now. Later I'll add test codes.