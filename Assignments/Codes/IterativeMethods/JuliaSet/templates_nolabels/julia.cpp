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
// TODO: implement F
        return fx;
    }
};

class DF {
    public:
    Mat operator()(Vec& x)
    {
        Mat dfx(2,2);
// TODO: implement the Jacobian of F
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

// TODO: analyze Newton iteration

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
