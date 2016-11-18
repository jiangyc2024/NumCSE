//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/QR>

#include <mgl2/mgl.h>

#include <figure/figure.hpp>

using namespace Eigen;

std::vector<size_t> order(const VectorXd &values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));
    std::sort(begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

VectorXd polyfit(const VectorXd &x, const VectorXd &y, size_t order)
{
    Eigen::MatrixXd A(x.size(), order+1);
    Eigen::VectorXd result;

    assert(x.size() == y.size());
    assert(x.size() >= order + 1);

    // Create matrix
    for (size_t i=0; i<x.size(); ++i) {
        for (size_t j=0; j<order+1; ++j) {
            A(i, j) = pow(x(i), j);
        }
    }

    // Solve for linear least squares fit
    VectorXd coeff = A.householderQr().solve(y);
    coeff.conservativeResize(order + 1);

    return coeff;
}

/* @brief Compute values of interpolant in knots $\Vx$ from $(t_i,y_i)$
 * @param[in] x Vector of knots
 * @param[in] t Vector of nodes
 * @param[in] y Vector of values of interpolant in nodes $\Vt$
 * @param[out] s Vector of values of interpolant in knots $\Vx$
 */
VectorXd PwLineIntp(const VectorXd &x, const VectorXd &t,
                    const VectorXd &y)

{
    assert(t.size() == y.size() &&
          "t and y must have same size!");

    // Initialization
    size_t n = t.size();
    auto t_indices = order(t);
    size_t m = x.size();
    auto x_indices = order(x);
    // You can also implement a solution which does not need
    // sorted vectors and e.g. for each knot $x_j$ looks
    // for the closest node $t_{i1}$ and the next closest node $t_{i2}$.
    // However, such solution will not become more efficient
    // if you give as input already sorted vectors: for each knot $x_j$
    // you will always have to iterate along the sorted vector $t$
    // to find the included node $t_i$.

    VectorXd s(m);

    // TODO:

    return s;
}



int main() {
// Compute convergence rate for interpolation by piecewise linear polyn.
// Uniform mesh in $[0,1]$, singular $f(t) = t^\alpha$, h-convergence
    {
        // Initialization
        size_t NumAlph = 15;
        size_t NumN = 50;
        VectorXd alphas = VectorXd::LinSpaced(NumAlph,0.1,2.9);
        VectorXd nn = VectorXd::LinSpaced(NumN,1,50); // Used nodes

        // Evaluation points
        VectorXd x = VectorXd::LinSpaced(1000,0,1);

        MatrixXd Err(NumAlph,NumN); // Error with max norm
        MatrixXd LocErr(NumAlph,NumN); // Location of maximal error
        for(size_t i=0; i<NumN; ++i) {
            size_t n = nn(i);

            // Nodes
            VectorXd t = VectorXd::LinSpaced(n+1,0,1);

            for(size_t j=0; j<NumAlph; ++j) {
                VectorXd s = x.array().pow(alphas(j));
                VectorXd y = t.array().pow(alphas(j));

                VectorXd p = PwLineIntp(x, t, y); // Interpolation

                size_t PosErr;
                Err(j,i) = (s - p).cwiseAbs().maxCoeff(&PosErr);

                // PosErr is index of point in $x$ with max error
                // LocErr is index of subinterval  with max error
                std::vector<double> tmp(t.size());
                VectorXd::Map(&tmp.front(), t.size()) = t -
                        x(PosErr)*VectorXd::Ones(t.size());
                LocErr(j,i) = count_if(tmp.begin(), tmp.end(),
                                       [] (double val) {return val <= 0;}) - 1;
                // "count\_if" only works when $t$ are already sorted!

            }
        }


        // Estimate convergence rate
        VectorXd rates(NumAlph);
        for(size_t i=0; i<NumAlph; ++i) {
            VectorXd coeff = polyfit(nn.array().log(), Err.row(i).array().log(), 1);
            rates(i) = -coeff(0);
        }

    }
{
// Compute convergence rate for interpolation by piecewise linear polyn.
// Beta-graded mesh in $[0,1]$, singular $f(t) = t^\alpha$, h-convergence
// for different betas

    // Initialization
    size_t NumAlph = 3;
    size_t NumBeta = 9;
    size_t NumN = 50;
    VectorXd alphas(NumAlph); alphas << 0.5, 0.75, 4/3;
    VectorXd betas = VectorXd::LinSpaced(NumBeta,1,50);
    VectorXd nn = VectorXd::LinSpaced(NumN,1,50); // Used nodes

    // Evaluation points
    VectorXd x = VectorXd::LinSpaced(1000,0,1);
    std::vector<MatrixXd> Err(NumAlph); // Error with max norm (3D data)
    std::vector<MatrixXd> LocErr(NumAlph); // Location of maximal error (3D data)
    for(size_t i=0; i<NumAlph; ++i) {
        Err[i] = MatrixXd(NumBeta,NumN);
        LocErr[i] = MatrixXd(NumBeta,NumN);
    }

    for(size_t i=0; i<NumN; ++i) {
        size_t n = nn(i);

        for(size_t j=0; j<NumBeta; ++j) {
            // Nodes
            VectorXd t = VectorXd::LinSpaced(n+1,0,1).array().pow(betas(j));

            for(size_t k=0; k<NumAlph; ++k) {
                VectorXd s = x.array().pow(alphas(k));
                VectorXd y = t.array().pow(alphas(k));

                VectorXd p = PwLineIntp(x, t, y); // Interpolation

                size_t PosErr;
                Err[k](j,i) = (s - p).cwiseAbs().maxCoeff(&PosErr);

                // PosErr is index of point in $x$ with max error
                // LocErr is index of subinterval  with max error
                std::vector<double> tmp(t.size());
                VectorXd::Map(&tmp.front(), t.size()) = t -
                               x(PosErr)*VectorXd::Ones(t.size());
                LocErr[k](j,i) = count_if(tmp.begin(), tmp.end(),
                               [] (double val) {return val <= 0;}) - 1;
                // "count\_if" only works when $t$ are already sorted!
            }
        }
    }

    VectorXd rates(NumAlph,NumBeta);
    for(size_t i=0; i<NumAlph; ++i) {

        // Estimate convergence rate
        for(size_t j=0; j<NumBeta; ++j) {
            size_t skip = 4;
            VectorXd row = Err[i].row(j);
            VectorXd coeff = polyfit(nn.tail(nn.size()-skip).array().log(),
                             row.tail(row.size()-skip).array().log(),1);
            rates(i) = -coeff(0);
        }

    }
}
{
// Plot of algebraically graded mesh

    // Initialization
    size_t n = 10;
    int beta = 2;

    VectorXd t  = VectorXd::LinSpaced(n+1,0,n)/n;
    VectorXd y  = t.array().pow(beta);
    VectorXd t_ = VectorXd::LinSpaced(101,0,1);
    VectorXd y_ = t_.array().pow(beta);

    mgl::Figure fig;
    fig.plot(t_, y_, "r");
    for(size_t i=0; i<n+1; ++i) {
        VectorXd t_tmp(3); t_tmp << t(i), t(i), 0;
        VectorXd y_tmp(3); y_tmp << 0, y(i), y(i);
        fig.plot(t_tmp, y_tmp, "b");
        fig.plot(t_tmp.tail(2), y_tmp.head(2), " r+");
    }
    fig.xlabel("uniform mesh");
    fig.ylabel("algeb. graded mesh, beta=2");
    fig.save("GradedMesh_cpp.eps");
}
}
