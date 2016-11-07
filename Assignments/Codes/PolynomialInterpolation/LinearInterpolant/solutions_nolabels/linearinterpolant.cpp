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
        data   i_points; //!< vector of i_points (t_i, y_i)
};

LinearInterpolant::LinearInterpolant(const data & i_points_) {
    assert( i_points_.size() > 1
            && "Must specify at least two nodes");
    // assignment of std::vectors in C++ performs deep copy
    i_points = i_points_;

    // This is needed in std::sort
    auto Ordering = [] (const pair & P, const pair & Q) -> bool {
        return P.first < Q.first;
    };
    // Sort the array
    std::sort(i_points.begin(), i_points.end(), Ordering);
}

double LinearInterpolant::operator() (double x) {
    // Lambda for compasiron of pair type,
    // we want to compare only the
    // first element of the pair, i.e. $t_i$
    auto Compare = [] (const pair & P, double V) -> bool { return P.first < V; };

    // Find the place $i$ (as iterator), wehre $t_i >= x$, notice that i_points
    // *must* be sorted (this is needed in std::lower_bound)
    // IMPORTANT: "lower_bound" performs a binary search
    // on the data, provided the data has a random access iterator
    // http://www.cplusplus.com/reference/iterator/RandomAccessIterator/
    // this is crucial here, since evaluation operators must be fast
    // this allows to reduce the complexity,
    // from linear complexity to a logarithmic
    auto it = std::lower_bound(i_points.begin(), i_points.end(), x, Compare);

    // Bound checks, if i = 0 and t_0 != x (we are before the first node)
    // or if x > t_n (we are after the last node).
    // In such cases return 0 (*it and *(it-1) would be undefined in such cases)
    if( ( it == i_points.begin() && it->first != x )
            || it == i_points.end() )
        return 0;

    // Actually compute the interpolated value, dist_rato contains the distance from t_{i-1} to x (as a ratio)
    double dist_ratio = (x - (it-1)->first) / (it->first - (it-1)->first);
    return (it-1)->second * (1 - dist_ratio) + it->second * dist_ratio;
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
