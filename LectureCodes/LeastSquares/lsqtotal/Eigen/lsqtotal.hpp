// computes only solution \Blue{$\Vx$} of fitted consistent LSE

#include <Eigen/Dense>
#include <iostream>

template<typename VecType, typename MatType>
VecType lsqtotal(const Eigen::MatrixBase<MatType> &A, const VecType &b)
{
    using namespace Eigen;
    using index_t = typename MatrixBase<MatType>::Index;
    using entry_t = typename MatrixBase<MatType>::Scalar;
    const index_t m(A.rows()); // No. of rows
    const index_t n(A.cols()); // No. of columns

    Matrix<entry_t, Dynamic, Dynamic> Ab(m, n + 1);
    Ab << A, b; // Ab = [A,b]
    auto V = Ab.jacobiSvd(ComputeThinU | ComputeThinV).matrixV(); // see \eqref{tlsq:1}

    entry_t s = V(n, n);
    if (s == 0)
    {
        std::cerr << "No solution" << std::endl;
        exit(1);
    }

    return -V.col(n).head(n) / s; // see \eqref{tlsq:3};
}