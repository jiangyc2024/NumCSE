// Solution of linear least squares problem \eqref{eq:LSQ1} by means of QR-decomposition
// Note: \Blue{$\VA\in\bbR^{m,n}$} with \Blue{$m>n$}, \Blue{$\operatorname{rank}(\VA) = n$} is assumed

#include <Eigen/Dense>

template<typename VecType, typename MatType>
void qrlsqsolve(const Eigen::MatrixBase<MatType> &A, const VecType &b,
                VecType &x, double &res)
{
    using namespace Eigen;
    using index_t = typename MatrixBase<MatType>::Index;
    using entry_t = typename MatrixBase<MatType>::Scalar;
    const index_t m(A.rows()); // No. of rows
    const index_t n(A.cols()); // No. of columns

    Matrix<entry_t, Dynamic, Dynamic> Ab(m, n + 1);
    Ab << A, b; // Ab = [A,b]

    decltype(Ab) R = Ab.householderQr()
        .matrixQR()
        .template triangularView<Upper>(); // R = triu(qr(Ab,0)), QR-decomposition of extended matrix \label{qrl:1}

    auto R_nn = R.block(0, 0, n, n).template triangularView<Upper>(); // R\_nn = R(1:n,1:n)
    // for the keyword "template", see http://eigen.tuxfamily.org/dox/TopicTemplateKeyword.html

    x = R_nn.solve(R.block(0, n, n, 1)); // R\_nn \ R(1:n,n+1), \Blue{$\wh{\Vx} = (\VR)_{1:n,1:n}^{-1}(\VQ^T\Vb)_{1:n}$}
    res = R(n, n); // \Blue{$= \N{\VA\wh{\Vx}-\Vb}_2$} (why ?) \label{qrl:2}
}