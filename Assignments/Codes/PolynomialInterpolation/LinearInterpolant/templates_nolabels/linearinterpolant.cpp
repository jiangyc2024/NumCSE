//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

#include <Eigen/Dense>

#include <figure.hpp>

using namespace Eigen;

/*!
 * \brief The LinearInterpolant class
 */
class LinearInterpolant {
    public:
        //! A "pair" holds the pair $(t_i, y_i)$ with the value of
        //! the interpolant $y_i$ at the node $t_i$
        using pair = std::pair<double, double>;
        //! "data" holds the list of "pair", may be ordered or not,
        //! but *cannot* contain duplicate $t_i$'s
        using data = std::vector<pair>;

        /*!
         * \brief LinearInterpolant Build interpolant from data
         * Sort the array for the first time:
         * the data is not assumed to be sorted
         * sorting is necessary for binary search
         * \param i_points_
         */
        LinearInterpolant(const data & i_points_);

        /*!
         * \brief operator () Evaluation operator.
         * Return the value of $I$ at $x$, i.e. $I(x)$.
         * Performs bound checks (i.e. if $x < t_0$ or $x >= t_n$ )
         * \param x
         * \return
         */
        double operator() (double x);
    private:
        // TODO: your data there
};

LinearInterpolant::LinearInterpolant(const data & i_points_) {
    // TODO: construct your data here
}

double LinearInterpolant::operator() (double x) {
    // TODO: your data there
}

int main(void) {
    // Test the class with the basis with nodes (-1,1,2,4)
    // and interpolant with values (-1,2,3,4).
    LinearInterpolant I(
        {
            {1, 2},
            {2, -1},
            {4, 4},
            {-1, -1}
        }
    );

    int n = 100;
    VectorXd x = VectorXd::LinSpaced(n, -1, 4);
    VectorXd y = x.unaryExpr(
                [&I] (double x) { return I(x); }
    );

    mgl::Figure fig;
    fig.title("Piecewise linear interpolation polynomial.");
    fig.plot(x, y, "b").label("I(x)");
    fig.xlabel("$x$");
    fig.ylabel("$y$");
    fig.legend(0, 1);
    fig.save("plp.eps");
    fig.save("plp.png");
}
