#include <Eigen/Dense>

template<typename VecType, typename MatType>
VecType lsqsvd(const Eigen::MatrixBase<MatType> &A, const VecType &b)
{
    using namespace Eigen;
    using entry_t = typename MatType::Scalar;
    JacobiSVD<Matrix<entry_t, Dynamic, Dynamic>> svd(A, ComputeThinU | ComputeThinV);

    auto sv = svd.singularValues();
    auto U = svd.matrixU();
    auto V = svd.matrixV();
    auto r = svd.nonzeroSingularValues(); // default threshold is NumTraits<Scalar>::epsilon()

    return V.leftCols(r)
        * (
            sv.head(r).cwiseInverse().asDiagonal()
                * (U.leftCols(r).adjoint() * b)
        );
}