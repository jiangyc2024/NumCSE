//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
# include <Eigen/Dense>
# include <cmath>
# include <iostream>
# include <mgl2/mgl.h>
# include "meshgrid.hpp"

typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;

class F {
    public:
    Vec operator()(Vec& x)
    {
        Vec fx(2);
        fx << x(0)*x(0)*x(0) - 3*x(0)*x(1)*x(1) - 1,
                3*x(0)*x(0)*x(1) - x(1)*x(1)*x(1);
        return fx;
    }
};

class DF {
    public:
    Mat operator()(Vec& x)
    {
        Mat dfx(2,2);
        dfx << 3*x(0)*x(0) - 3*x(1)*x(1),
                -6*x(0)*x(1),
                6*x(0)*x(1),
                3*x(0)*x(0) - 3*x(1)*x(1);
        return dfx;
    }
};

int main()
{
    // Exact roots of f(z) = z^3 - 1, z \in C
    Vec z1(2), z2(2), z3(2);
    z1 << 1, 0;
    z2 << -0.5, 0.5*std::sqrt(3);
    z3 << -0.5, -0.5*std::sqrt(3);
    const unsigned int maxit = 20;
    const double tol = 1e-4;
    const unsigned int N = 1000;
    Vec x = Vec::LinSpaced(N, -2, 2);

    Mat X, Y;
    meshgrid(x, x, X, Y);

    Vec C = Vec::Ones(X.size());

    F Func; DF Jac;
    for (int i = 0; i < X.size(); ++i){
        Vec v(2); v << *(X.data() + i), *(Y.data() + i);

        // Newton's iteration
        for (unsigned int k = 1; k <= maxit; ++k){
            v -= Jac(v).lu().solve(Func(v));

            // Termination criterium: stop when close to one of the roots
            if ((v - z1).norm() < tol){
                C(i) = 1 + k;
                break;
            } else if ((v - z2).norm() < tol){
                C(i) = 1 + k + maxit;
                break;
            } else if ((v - z3).norm() < tol){
                C(i) = 1 + k + 2*maxit;
                break;
            }
        }
    }

    // Normalize results for plot
    C = (C.array()/double(C.maxCoeff())).matrix();

    mglData Xd(X.rows(), X.cols(), X.data());
    mglData Yd(Y.rows(), Y.cols(), Y.data());
    mglData Cd(X.rows(), X.cols(), C.data());

    mglGraph gr;
    gr.SubPlot(1,1,0,"<_");
    gr.SetRanges(-2,2,-2,2);
    gr.SetRange('c', 0, 1);
    gr.Colorbar("bcwyr");
    gr.Title("Juliaset for \\z^3 = 1");
    gr.Axis();
    gr.Tile(Xd, Yd, Cd, "bcwyr");
    gr.WriteEPS("juliaset.eps");
}
