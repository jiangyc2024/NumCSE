// Solves constrained linear least squares problem \eqref{clsq} with \texttt{dim} passing \Blue{$d$}

#include <Eigen/Dense>
#include <iostream>

template<typename VecType, typename MatType, typename Index, typename Scalar>
int clsq(const Eigen::MatrixBase<MatType> &A, const Index dim, Scalar &c, VecType &n)
{
    using namespace Eigen;
    using namespace std;

    const Index p = A.cols();
    Index m = A.rows();

    if (p < dim + 1)
    {
        cerr << "not enough unknowns" << endl;
        return 1;
    }
    if (m < dim)
    {
        cerr << "not enough equations" << endl;
        return 2;
    }

    m = min(m, p);
    Matrix<Scalar, Dynamic, Dynamic> R = A.householderQr()
        .matrixQR()
        .template triangularView<Upper>(); // First step: orthogonal transformation, see Code~\ref{mc:qrlsqsolve}

    auto V = R.block(p - dim, p - dim, m + dim - p, dim)
        .jacobiSvd(ComputeThinU | ComputeThinV)
        .matrixV(); // Solve \eqref{eq:HRmin2}

    n = V.col(dim - 1);

    auto R_topleft = R.topLeftCorner(p - dim, p - dim).template triangularView<Upper>();
    c = -(R_topleft.solve(R.block(0, p - dim, p - dim, dim)) * n)(0);

    return 0;
}