//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
# include <Eigen/Dense>
# include <cmath>
# include <iostream>
# include <string>
# include <mgl2/mgl.h>
# include "meshgrid.hpp"

typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;

class Func {
    public:
    Mat operator()(Mat& eps, Mat tau, double C, double p)
    {
        Mat kmin;

// TODO: implement $k_{\min}$

        return kmin;
    }
};

int main()
{
    // Initialization
    double C = 2;
    double p = 1.5;
    double eps_max = std::pow(C, 1/(1-p));
    unsigned ngp = 100; // Number of grid points

    Vec eps_lin = Vec::LinSpaced(ngp, 0, eps_max);
    Vec tau_lin = Vec::LinSpaced(ngp, 0, eps_max);
    Vec eps_lin_seg = eps_lin.segment(1,eps_lin.size()-2);
    Vec tau_lin_seg = tau_lin.segment(1,tau_lin.size()-2);
    Mat eps_msh, tau_msh;
    meshgrid(eps_lin_seg, tau_lin_seg, eps_msh, tau_msh);

    Func kmin;

// TODO: compute $k_{\min}$
    Mat k_ = Mat::Zero(eps_msh.rows(), eps_msh.cols());

    // Normalize results for plot
    Mat k_norm = (k_.array()/double(k_.maxCoeff())).matrix();

    mglData eps(eps_msh.rows(), eps_msh.cols(), eps_msh.data());
    mglData tau(tau_msh.rows(), tau_msh.cols(), tau_msh.data());
    mglData k(k_norm.rows(), k_norm.cols(), k_norm.data());

    mglGraph gr;
    gr.SetRanges(eps.Minimal(), eps.Maximal(), tau.Minimal(), tau.Maximal());
    gr.Title("Minimal number of iterations for error < tau");
    gr.Axis(); gr.Label('x',"epsilon_0",0); gr.Label('y',"tau",0);
    gr.Tile(eps, tau, k, "kRryw");
    gr.SetRange('c', k_.minCoeff(), k_.maxCoeff());
    gr.SetTicksVal('c', " ");
    for(unsigned i=k_.minCoeff(); i<=k_.maxCoeff(); ++i) {
        gr.AddTick('c', i, std::to_string(i).c_str());
    }
    gr.Colorbar("kRryw");
    gr.WriteEPS("k_min_cpp.eps");
}
