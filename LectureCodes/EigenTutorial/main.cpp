// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

// Functions demonstrating the initialization and use of matrices and vectors
// with Eigen.

/** Demo for use of matrix and vector types in Eigen */
/* SAM_LISTING_BEGIN_1 */
template <typename Scalar> void eigenTypeDemo(unsigned int dim) {
  using dynMat_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using dynColVec_t = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using dynRowVec_t = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;
  using index_t = typename dynMat_t::Index;
  using entry_t = typename dynMat_t::Scalar;

  dynColVec_t colvec(dim);
  dynRowVec_t rowvec(dim);
  for (index_t i = 0; i < colvec.size(); ++i)
    colvec[i] = (entry_t)i;
  for (index_t i = 0; i < rowvec.size(); ++i)
    rowvec[i] = (entry_t)1 / (i + 1);
  colvec[0] = (entry_t)3.14;
  rowvec[dim - 1] = (entry_t)2.718;
  cout << "colvec = " << colvec << endl << "rwovec = " << rowvec << endl;
  dynMat_t vecprod = colvec * rowvec;
  const int nrows = vecprod.rows();
  const int ncols = vecprod.cols();
  cout << "vecprod is " << nrows << " x " << ncols << "-matrix" << endl;
}
/* SAM_LISTING_END_1 */

/** Demo for initialization of matrices */
/* SAM_LISTING_BEGIN_2 */
void blockInit(int size = 6) {
  size = 2 * (size / 2);
  cout << "Matrix block initialization" << endl;

  MatrixXd mat(size, size);
  mat << MatrixXd::Zero(size / 2, size / 2),
      MatrixXd::Identity(size / 2, size / 2),
      MatrixXd::Identity(size / 2, size / 2),
      MatrixXd::Zero(size / 2, size / 2);
  std::cout << "#1 " << mat << std::endl;

  MatrixXd botRows(size / 2, size);
  for (int l = 0; l < botRows.size(); ++l)
    botRows(l) = l;
  mat << MatrixXd::Zero(size / 2, size / 2),
      MatrixXd::Identity(size / 2, size / 2), botRows;
  std::cout << "#2 " << mat << std::endl;

  mat << MatrixXd::Constant(size, size - 4, 1.5),
      MatrixXd::Constant(size, 4, 3.5);
  std::cout << "#3 " << mat << std::endl;
  mat << MatrixXd::Constant(size - 2, size - 4, 1.5), // top row, first block
      MatrixXd::Constant(size - 2, 3, 3.5),           // top row, second block
      MatrixXd::Constant(size - 2, 1, 7.5),           //  top row, third block
      MatrixXd::Constant(2, size - 2, 2.5),           // bottom row, left block
      MatrixXd::Constant(2, 2, 4.5);                  // bottom row, right block
  std::cout << "#4 " << mat << std::endl;
}
/* SAM_LISTING_END_2 */

/** Demo for accessing sub-matrices */
/* SAM_LISTING_BEGIN_3 */
void accessColRow(int nrows = 6, int ncols = 7) {
  using index_t = typename Eigen::MatrixXd::Index;
  cout << "Column and row access" << endl;
  // Allocate dynamic \eigen matrix of type double
  MatrixXd m(nrows, ncols);
  // Initialization by direct component access
  for (index_t i = 0; i < m.rows(); i++)
    for (index_t j = 0; j < m.cols(); j++)
      m(i, j) = i + j;
  // Print matrix to standard output
  cout << "Matrix m = " << endl << m << endl;
  // Print rows and columns
  for (index_t l = 0; l < m.rows(); l++)
    cout << "Row " << l << " = " << m.row(l) << endl;
  for (index_t l = 0; l < m.cols(); l++)
    cout << "Col " << l << " = " << m.col(l) << endl;
  // Access rows and columns as vectors
  RowVectorXd row1(m.row(1));
  VectorXd col1(m.col(1));

  cout << "Tensor product" << endl << col1 * row1 << endl;

  m.col(2).setZero();
  m.row(2) = RowVectorXd::Constant(ncols, -1);
  m.row(4).tail(3) = RowVectorXd::Constant(3, 1.5);
  cout << "Modified matrix m = " << m << endl;
}
/* SAM_LISTING_END_3 */

template <class Matrix> void blockAccess(Eigen::MatrixBase<Matrix> &);

void blockOps(int nrows = 6, int ncols = 7) {
  // Allocate uninitialised matrix
  MatrixXd M(nrows, ncols);
  // Fill matrix by accessing entries directly
  for (int i = 0; i < M.rows(); i++)
    for (int j = 0; j < M.cols(); j++)
      M(i, j) = i - j;
  blockAccess(M);
}

/** @brief Access to a matrix block
    @tparam Matrix any Eigen matrix compatible type
    @param M reference to a matrix-like object.

    Since Eigen uses expression templates, the expressions often evaluate to
    temporary objects. This function demonstrates a way to accept such objects
    as arguments of functions
*/
/* SAM_LISTING_BEGIN_4 */
template <typename Matrix> void blockAccess(Eigen::MatrixBase<Matrix> &M) {
  using index_t = typename Eigen::MatrixBase<Matrix>::Index;
  using entry_t = typename Eigen::MatrixBase<Matrix>::Scalar;
  const index_t nrows(M.rows());
  const index_t ncols(M.cols());

  cout << "Accessing matrix blocks" << endl;
  // Print matrix
  cout << "Matrix M = " << M << endl;
  // Set block size
  int p = nrows / 2, q = ncols / 2;
  // Output submatrix with left upper entry at position \texttt{(i,i)}
  for (index_t i = 0; i < min(p, q); i++)
    cout << "Block (" << i << ',' << i << ',' << p << ',' << q
         << ") = " << M.block(i, i, p, q) << endl;
  // Modify sub-matrix by adding a constant
  M.block(1, 1, p, q) += MatrixXd::Constant(p, q, 1.0);
  cout << "With modified block (1,1," << p << ',' << q << "): M = " << M
       << endl;
  // Extract sub-matrix
  MatrixXd B(M.block(1, 1, p, q));
  cout << "Isolated modified block = " << B << endl;
  // Special sub-matrices
  cout << p << " top rows of m = " << M.topRows(p) << endl;
  cout << p << " bottom rows of m = " << M.bottomRows(p) << endl;
  cout << q << " left cols of m = " << M.leftCols(q) << endl;
  cout << q << " right cols of m = " << M.rightCols(p) << endl;
  // Not possible: requires a vector
  // cout << (p+q) << " head elements of m = " << m.head(p+q) << endl;
  // cout << (p+q) << " tail elements of m = " << m.tail(p+q) << endl;
  const MatrixXd T = M.template triangularView<Upper>();
  cout << "Upper triangular part = " << endl << T << endl;
  M.template triangularView<Lower>() *= -1.5;
  cout << "Matrix M = " << endl << M << endl;
}
/* SAM_LISTING_END_4 */

/** Explores internal matrix storage format of Eigen */
/* SAM_LISTING_BEGIN_5 */
void storageOrder(int nrows = 6, int ncols = 7) {
  cout << "Demonstrating different storage orders for Eigen matrices" << endl;
  // Template parameter \texttt{ColMajor} selects column major data layout
  Matrix<double, Dynamic, Dynamic, ColMajor> mcm(nrows, ncols);
  // Template parameter \texttt{RowMajor} selects row major data layout
  Matrix<double, Dynamic, Dynamic, RowMajor> mrm(nrows, ncols);
  // Direct initialization; lazy option: use \texttt{int} as index type
  for (int l = 1, i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++, l++)
      mcm(i, j) = mrm(i, j) = l;

  cout << "Matrix mrm = " << endl << mrm << endl;
  cout << "mcm linear = ";
  for (int l = 0; l < mcm.size(); l++)
    cout << mcm(l) << ',';
  cout << endl;

  cout << "mrm linear = ";
  for (int l = 0; l < mrm.size(); l++)
    cout << mrm(l) << ',';
  cout << endl;

  // Retrieve pointer to raw matrix data
  double *mdat = mcm.data();
  for (int l = 0; l < mcm.size(); l++)
    cout << mdat[l] << ',';
  cout << endl;

  Map<MatrixXd> mrs(mdat, nrows / 2, ncols * 2);
  cout << "reshaped = " << endl << mrs << endl;
  mrs *= -1.5;
  cout << "mrm = " << endl << mrm << endl;
  cout << "mrs = " << endl << mrs << endl;
}
/* SAM_LISTING_END_5 */

// test of [] access operatror
void veccompaccess(int len = 6) {
  cout << "Access to vector compoents by means of []" << endl;
  VectorXd v(VectorXd::Random(len));
  RowVectorXd rv(RowVectorXd::Random(len));

  // Access through [] operator
  for (int l = 0; l < v.size(); l++)
    cout << v[l] - rv[l] << ',';
  cout << endl;

  /* This triggers a compilation error. [] component access is available only
     for vectors:
  MatrixXd M(MatrixXd::Random(len,len));
  double x = M[2];

  /opt/local/include/eigen3/Eigen/src/Core/DenseCoeffsBase.h:375:7: error:
      static_assert failed
      "THE_BRACKET_OPERATOR_IS_ONLY_FOR_VECTORS__USE_THE_PARENTHESIS_OPERATOR_INSTEAD"
      EIGEN_STATIC_ASSERT(Derived::IsVectorAtCompileTime,
  */
}

void matArray(int nrows = 6, int ncols = 7) {
  MatrixXd m1(nrows, ncols), m2(nrows, ncols);
  for (int i = 0; i < m1.rows(); i++)
    for (int j = 0; j < m1.cols(); j++) {
      m1(i, j) = (double)(i + 1) / (j + 1);
      m2(i, j) = (double)(j + 1) / (i + 1);
    }
  cout << "Matrix m1 = " << endl << m1 << endl;
  cout << "Matrix m2 = " << endl << m2 << endl;
  // Compare two matrices entrywise
  cout << "m1 > m2 in " << (m1.array() > m2.array()).count() << " entries"
       << endl;
  // Selective access to entries of matrices
  // MatrixXd mx(m1.array() > 1.0);
  // cout << "mx = " << mx << endl;
  // Entrywise product
  MatrixXd m3 = (m1.array() * m2.array()).matrix();
  cout << "Entrywise product m3 = " << endl << m3 << endl;
  cout << "cwiseproduct: |diff| = " << (m3 - m1.cwiseProduct(m2)).norm()
       << endl;

  // Logarithm
  cout << "Log(m1) = " << endl << log(m1.array()) << endl;
  // Count entries larger than 3
  cout << (m1.array() > 3).count() << " entries of m1 > 3" << endl;
  // Apply a (lambda) function to all entries of a matrix
  auto fnct = [](double x) { return (x + 1.0 / x); };
  m1 = m1.unaryExpr(fnct);
  cout << "f(m1) = " << endl << m1 << endl;
}

template <typename MatType> void reshape(const MatType &M, int k, int l) {
  using index_t = typename MatType::Index;
  using entry_t = typename MatType::Scalar;
  const index_t nrows(M.rows());
  const index_t ncols(M.rows());
  const index_t nsize(M.size());

  eigen_assert(k * l == nsize);
  Matrix<entry_t, Dynamic, Dynamic> R;
  for (index_t l = 0; l < nsize; l++)
    R(l) = M(l);
}

template <typename MatType> void reshapetest(MatType &M) {
  using index_t = typename MatType::Index;
  using entry_t = typename MatType::Scalar;
  const index_t nrows(M.rows());
  const index_t ncols(M.rows());
  const index_t nsize(M.size());

  // reshaping possible only for matrices with non-prime dimensions
  if ((nsize % 2) == 0) {
    entry_t *Mdat = M.data(); // raw data array for M
    // Reinterpretation of data of M
    Map<Eigen::Matrix<entry_t, Dynamic, Dynamic>> R(Mdat, 2, nsize / 2);
    // (Deep) copy data of M into matrix of different size
    Eigen::Matrix<entry_t, Dynamic, Dynamic> S =
        Map<Eigen::Matrix<entry_t, Dynamic, Dynamic>>(Mdat, 2, nsize / 2);

    cout << "Matrix M = " << endl << M << endl;
    cout << "reshaped to " << R.rows() << 'x' << R.cols() << " = " << endl
         << R << endl;
    // Modifying R affects M, because they share the data space !
    R *= -1.5;
    cout << "Scaled (!) matrix M = " << endl << M << endl;
    // Matrix S is not affected, because of deep copy
    cout << "Matrix S = " << endl << S << endl;
  }
}

void rsCall(int nrows = 6, int ncols = 7) {
  // Allocate uninitialised matrix
  MatrixXd M(nrows, ncols);
  // Fill matrix by accessing entries directly
  for (int i = 0; i < M.rows(); i++)
    for (int j = 0; j < M.cols(); j++)
      M(i, j) = i - j;
  reshapetest(M);
  // reshape(2*M,2,(nrows*ncols)/2);
}

void rowConcat(void) {
  RowVectorXd x0(77);
  x0 << RowVectorXd::LinSpaced(26, 0.0, 1.0),
      RowVectorXd::LinSpaced(51, 2.0, 3.0);
  cout << "x0 = " << x0 << endl;
}

int main(int argc, char **argv) {
  cout << "EIGEN TUTORIAL CODES" << endl;
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <selection>" << endl;
    cerr << "1: blockInit()" << std::endl;
    cerr << "2: accessColRow()" << std::endl;
    cerr << "3: blockOps()" << std::endl;
    cerr << "4: storageOrder(4, 3)" << std::endl;
    cerr << "5:  matArray()" << std::endl;
    cerr << "6: eigenTypeDemo<float>(7)" << std::endl;
    cerr << "7: rtsCall()" << std::endl;
    cerr << "8: rowConcat()" << std::endl;
    exit(-1L);
  } else {
    const int sel = atoi(argv[1]);
    switch (sel) {
    case 1: {
      blockInit();
      break;
    }
    case 2: {
      accessColRow();
      break;
    }
    case 3: {
      blockOps();
      break;
    }
    case 4: {
      storageOrder(4, 3);
      break;
    }
    case 5: {
      matArray();
      break;
    }
    case 6: {
      eigenTypeDemo<float>(7);
      break;
    }
    case 7: {
      rsCall();
      break;
    }
    case 8: {
      rowConcat();
      break;
    }
    default: {
      cerr << "Invalid selection" << endl;
      exit(-1L);
    }
    }
    exit(0);
  }

  {
    cout << "Output of a small dense dynamic matrix" << endl;
    MatrixXd m(2, 2);
    //
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    cout << m << endl;
  }
  {
    cout << "Matrix x vector" << endl;
    MatrixXd m = MatrixXd::Random(3, 3);
    m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;
    cout << "m =" << endl << m << endl;
    VectorXd v(3);
    v << 1, 2, 3;
    cout << "m * v =" << endl << m * v << endl;
  }
  {
    cout << "Matrix x Vector with fixed size matrices" << endl;
    Matrix3d m = Matrix3d::Random();
    m = (m + Matrix3d::Constant(1.2)) * 50;
    cout << "m =" << endl << m << endl;
    Vector3d v(1, 2, 3);

    cout << "m * v =" << endl << m * v << endl;
  }
  {
    // Declaration of various matrix types, calling constructors with specifying
    // the size
    int n = 3;
    Matrix<double, 3, 3> matFixed33;
    Matrix<double, Dynamic, Dynamic> matDynamic(n, n);
    using mat3rows_t = Matrix<double, 3, Dynamic>;
    mat3rows_t matDynCols(3, n);
    Matrix<double, Dynamic, 3> matDynRows(n, 3);

    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        matDynamic(i, j) = 1.0 / (i + j + 1.0);
      }
    cout << matDynamic.rows() << 'x' << matDynamic.cols()
         << "-matrix: " << matDynamic << endl;
    double s = 0.0;
    for (int l = 0; l < matDynamic.size(); l++)
      s += matDynamic(l);
    cout << "Sum of all matrix entries = " << s << endl;
  }
  {
    const int worldDim = 2;
    using pointCoords_t = Eigen::Matrix<double, worldDim, 4>;
    pointCoords_t pts(void); // creates uninitialized $2\times 4$-matrix
    int n = 7;
    using dynMat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    dynMat_t mat(n, 2 * n); //
    using varCoords_t = Eigen::Matrix<double, worldDim, Eigen::Dynamic>;
    varCoords_t npts(worldDim, n + 1);
  }

  return 0;
}
