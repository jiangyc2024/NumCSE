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

    Vector2d x;
    broyd_history_t<double> broyd_history;
    upbroyd_history_t<double> upbroyd_history;
    std::tie(x, upbroyd_history) = upbroyd(F, x0, DF(x0), 0.000001, 20);
    std::tie(x, broyd_history) = broyd(F, x0, DF(x0), 0.000001, 20);

    return 0;
}
