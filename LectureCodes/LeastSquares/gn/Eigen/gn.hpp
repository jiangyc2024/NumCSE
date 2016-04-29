#include <Eigen/Dense>
#include <functional>

template<typename VecType, typename MatType>
VecType gn(VecType x,
           std::function<VecType(const VecType &)> F,
           std::function<MatType(const VecType &)> J,
           double tol) {
    using namespace Eigen;

    VecType s = J(x).householderQr().solve(F(x)); // \label{gn:2}
    x = x - s;
    while (s.norm() > tol * x.norm()) // \label{gn:term}
    {
        s = J(x).householderQr().solve(F(x)); // \label{gn:5}
        x = x - s;
    }

    return x;
}
