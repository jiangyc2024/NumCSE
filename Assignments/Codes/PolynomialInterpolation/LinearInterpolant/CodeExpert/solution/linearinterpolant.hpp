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
    data   i_points; //!< vector of i_points (t_i, y_i)
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
  
    // helper lambda to transform two Eigen::VectorXd into one std::vector<std::pair<double, double>>
    auto VecToVecPair = [](const double &t_i, const double &y_i) -> pair {
      return std::make_pair(t_i, y_i);
    };
    
    // perform transformation towards std::vector<std::pair<double, double>>
    i_points.resize(t.size());
    std::transform(t.data(), t.data() + t.size(), y.data(), i_points.begin(), VecToVecPair);

    // This is needed in std::sort
    auto Ordering = [] (const pair & P, const pair & Q) -> bool {
        return P.first < Q.first;
    };
    // Sort the array
    std::sort(i_points.begin(), i_points.end(), Ordering);
    // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double LinearInterpolant::operator() (double x) {
    double I_x = 0.;
    // TODO: (6-4.a) Implement an efficient evaluation operator.
    // START
    
    // Lambda for comparison of pair type,
    // we want to compare only the
    // first element of the pair, i.e. $t_i$
    auto Compare = [] (const pair & P, double V) -> bool { return P.first < V; };

    // Find the place $i$ (as iterator), wehre $t_i >= x$, notice that i\_points
    // *must* be sorted (this is needed in std::lower\_bound)
    // IMPORTANT: "lower\_bound" performs a binary search
    // on the data, provided the data has a random access iterator
    // http://www.cplusplus.com/reference/iterator/RandomAccessIterator/
    // this is crucial here, since evaluation operators must be fast
    // this allows to reduce the complexity,
    // from linear complexity to a logarithmic
    auto it = std::lower_bound(i_points.begin(), i_points.end(), x, Compare);

    // Bound checks, if $i = 0$ and $t_0 != x$ (we are before the first node)
    // or if $x > t_n$ (we are after the last node).
    // In such cases return 0 (*it and *(it-1) would be undefined in such cases)
    if( ( it == i_points.begin() && it->first != x )
            || it == i_points.end() )
        return 0;

    // Actually compute the interpolated value, dist\_ratio
    // contains the distance from $t_{i-1}$ to $x$ (as a ratio)
    double dist_ratio = (x - (it-1)->first) / (it->first - (it-1)->first);
    I_x = (it-1)->second * (1 - dist_ratio) + it->second * dist_ratio;
    // END
    return I_x;
}
/* SAM_LISTING_END_2 */

#endif
