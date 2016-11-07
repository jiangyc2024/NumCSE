//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <Eigen/Dense>

#include "polyfit.hpp"
#include "polyval.hpp"

#include <figure.hpp>

using namespace Eigen;

/*!
 * \brief dipoleval
 * \param t
 * \param y
 * \param x
 * \return
 */
VectorXd dipoleval(const VectorXd & t,
                   const VectorXd & y,
                   const VectorXd & x) {
    assert(t.size() = y.size()
           && "t and y must have same size!");

    VectorXd ret;
    ret.resizeLike(x);

    // TODO
    return ret;
}

/*!
 * \brief dipoleval_alt
 * \param t
 * \param y
 * \param x
 * \return
 */
VectorXd dipoleval_alt(const VectorXd & t,
                       const VectorXd & y,
                       const VectorXd & x) {
    assert(t.size() = y.size()
           && "t and y must have same size!");

    int n = y.size();

    VectorXd P;
    // TODO
    return P;
}

int main() {

    int n = 10;
    int N = 100;

    VectorXd t = VectorXd::LinSpaced(n,0,3);
    VectorXd y = t.unaryExpr([] (double x) { return std::sin(x);});

    VectorXd x = VectorXd::LinSpaced(N,0,3);

    VectorXd P = polyval(polyfit(t,y,n-1), x);
    VectorXd dp = dipoleval(t,y,x);
    VectorXd dp_alt = dipoleval_alt(t,y,x);

    std::cout << "Error: "
              << (dp-dp_alt).norm()
              << std::endl;

    mgl::Figure fig;
    fig.title("Derivative of interpolation polynomial.");
    fig.plot(x, dp, "r").label("p'(t) (AN)");
    fig.plot(x, dp_alt, "g").label("p'(t) (polyfit)");
    fig.plot(t, y, " b+").label("interpolation pts.");
    fig.plot(x, P, "b").label("p(t)");
    fig.xlabel("$x$");
    fig.ylabel("$y$");
    fig.legend(0, 1);
    fig.save("eval_deriv.eps");
    fig.save("eval_deriv.png");
}
