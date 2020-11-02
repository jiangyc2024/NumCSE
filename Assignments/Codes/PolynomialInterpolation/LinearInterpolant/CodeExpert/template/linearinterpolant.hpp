#ifndef LINEARINTERPOLANT_HPP
#define LINEARINTERPOLANT_HPP

#include <vector>
#include <algorithm>
#include <cassert>

#include <Eigen/Core>

/*!
 * \brief The LinearInterpolant class
 */
/* SAM_LISTING_BEGIN_0 */
class LinearInterpolant {
public:
    //! A "pair" holds the pair $(t_i, y_i)$ with the value of
    //! the interpolant $y_i$ at the node $t_i$
    using pair = std::pair<double, double>;
    //! "data" holds the list of "pair", may be ordered or not,
    //! but *cannot* contain duplicate $t_i$'s
    using data = std::vector<pair>;

    /*!
     * \brief LinearInterpolant builds interpolant from data.
     * Sort the array for the first time:
     * the data is not assumed to be sorted
     * sorting is necessary for binary search
     * \param t x-values of interpolation points
     * \param y y-values of interpolation points
	 */
    LinearInterpolant(const Eigen::VectorXd &t, const Eigen::VectorXd &y);

    /*!
     * \brief operator () Evaluation operator.
     * Return the value of $I$ at $x$, i.e. $I(x)$.
     * Performs bound checks (i.e. if $x < t_0$ or $x >= t_n$ ).
     * \param x Value $x \in \mathbb{R}$.
     * \return Value $I(x)$.
     */
    double operator() (double x);
private:
    // TODO: (6-4.a) Put necessary members of the class here.
    // START
    
    // END
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
LinearInterpolant::LinearInterpolant(const Eigen::VectorXd &t, const Eigen::VectorXd &y) {
    assert(t.size() == y.size() && "Input sizes must match");
    assert(t.size() > 1
            && "Must specify at least two nodes");
  
    // TODO: (6-4.a) Implement the constructor of the class. t and y may not be sorted.
    // START
    
    // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double LinearInterpolant::operator() (double x) {
    double I_x = 0.;
    // TODO: (6-4.a) Implement an efficient evaluation operator.
    // START
    
    // END
    return I_x;
}
/* SAM_LISTING_END_2 */

#endif
