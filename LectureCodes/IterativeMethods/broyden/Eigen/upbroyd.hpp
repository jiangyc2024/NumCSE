#ifndef UPBROYD_HPP
#define UPBROYD_HPP

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <tuple>
#include <utility>

using namespace Eigen;

// convenience typedef
template <typename T, int N>
using Vector = Eigen::Matrix<T, N, 1>;

template <typename Scalar>
using upbroyd_history_t = std::vector<std::tuple<unsigned, Scalar, Scalar, Scalar>>;

/**
 * \brief Good Broyden rank-1-update quasi-Newton method
 * Straightforward implementation for small problems
 * \param F Non-linear mapping in n dimensions
 * \param x initial guess
 * \param J initial guess for Jacobi matrix at x0
 * \param tol tolerance for termination
 */
template <typename FuncType, typename JacType, typename Scalar=double, int N=Dynamic>
std::pair<Vector<Scalar, N>, upbroyd_history_t<Scalar>> upbroyd(const FuncType &F, Vector<Scalar, N> x, JacType J, 
                                                                const Scalar tol, const unsigned maxit=20)
{
    // calculate LU factorization
    auto fac = J.lu();

    unsigned k = 1;
    Vector2d s = fac.solve(F(x));
    x -= s;
    auto sn = s.squaredNorm();
    std::cout << "Iteration(UPD) " << std::scientific << k << ": |s| = " << s.norm()
              << ", |F(x)| = " << F(x).norm() << std::endl;


    // keeping a record of the convergence history
    std::vector<std::tuple<unsigned, Scalar, Scalar, Scalar>> history;
    history.emplace_back(k, s.norm(), F(x).norm(), 1);

    // update vector
    std::vector<VectorXd> dx = {s};
    std::vector<Scalar> dxn = {sn};

    while ((sn > tol*tol) && (k < maxit)) {
        Vector<Scalar, N> w = fac.solve(F(x));
        for (unsigned l=1; l<k; ++l) {
            w = w + dx[l]*(dx[l-1].dot(w))/dxn[l-1];
        }
        auto z = s.dot(w);
        s = (1+z/(sn-z))*w;
        sn = s.squaredNorm();
        dx.push_back(s);
        dxn.push_back(sn);
        x -= s;
        history.emplace_back(k, std::sqrt(sn), F(x).norm(), w.norm()/sqrt(dxn[k-1]));
        std::cout << "Iteration(UPD) " << std::scientific << k << ": |s| = " << s.norm()
                  << ", |F(x)| = " << F(x).norm()
                  << ", theta = " << std::fixed << w.norm()/sqrt(dxn[k-1]) << std::endl;
        ++k;
    }
    return std::make_pair(x, history);
}

#endif
