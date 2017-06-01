///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <tuple>
#include <utility>
#include <figure/figure.hpp>

#include "broyd.hpp"
#include "upbroyd.hpp"
#include "../../newton/Eigen/newton.hpp"
#include "../../simpnewton/Eigen/simpnewton.hpp"

using namespace Eigen;

int main() {
    auto F = [] (const Eigen::Vector2d& x) -> Vector2d {
        return Vector2d(pow(x[0], 2)-pow(x[1], 4), 
                x[0] - pow(x[1], 3));
    };
    auto DF = [] (const Eigen::Vector2d& x) -> Matrix2d {
        Matrix2d res;
        res <<  2.*x[0],    -4.*pow(x[1], 3), 
                1      ,    -3.*pow(x[1], 2);
        return res;
    };

    Vector2d x0(0.7, 0.7);

    /*
     * Broyden
     */
    {
        // this function is called in every iteration and prints the progress
        auto print_progress_cb = [] (unsigned k, Vector2d x, Vector2d f, Vector2d s) {
            std::cout << "Iteration " << k << ": |s| " << s.norm() << " " << f.norm() << std::endl;
        };
        // run the algorithm (we ignore the result)
        broyd(F, x0, DF(x0), 0.000001, 20, print_progress_cb);
    }

    /*
     * Upbroyd
     */
    {
        // this function is called in every iteration and prints the progress
        auto print_progress_cb = [] (unsigned k, Vector2d x, Vector2d f, Vector2d s, Vector2d w, std::vector<double>& dxn) {
            std::cout << "Iteration(UPD) " << std::scientific << k << ": |s| = " << s.norm()
                  << ", |F(x)| = " << f.norm();
            if (dxn.size() > 1)
                std::cout << ", theta = " << std::fixed << w.norm()/std::sqrt(dxn[k-1]);
            std::cout << std::endl;
        };
        // run the algorithm (we ignore the result)
        upbroyd(F, x0, DF(x0), 0.000001, 20, print_progress_cb);
    }

    return 0;
}
