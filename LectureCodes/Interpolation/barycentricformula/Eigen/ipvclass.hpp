///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair based on code by J. Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

# include <Eigen/Dense>

/* SAM_LISTING_BEGIN_0 */
template <typename NODESCALAR = double>
class BarycPolyInterp {
private:
  using nodeVec_t = Eigen::Matrix<NODESCALAR,Eigen::Dynamic,1>; 
  using idx_t = typename nodeVec_t::Index;
  // Number \Blue{$n$} of interpolation points, deg polynomial +1
  const idx_t n;
  // Locations of \Blue{$n$} interpolation points
  nodeVec_t t;
  // Precomputed values \Blue{$\lambda_i$}, \Blue{$i=0,\ldots,n-1$}
  nodeVec_t lambda;
public:
  // Constructors taking node vector \Blue{$[t_{0},\ldots,t_{n}]^{\top}$} as argument
  BarycPolyInterp(const nodeVec_t &_t);
  // The interpolation points may also be passed in an STL container
  template <typename SeqContainer> 
    BarycPolyInterp(const SeqContainer &v); 
  // Computation of \Blue{$p(x_{k})$} for data values \Blue{$(y_{0},\ldots,y_{n})$} and
  // evaluation points \Blue{$x_{k}$}
  template <typename RESVEC,typename DATAVEC>
  RESVEC eval(const DATAVEC &y,const nodeVec_t &x) const;
private:
  void init_lambda(void);
 };  
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename NODESCALAR>
BarycPolyInterp<NODESCALAR>::BarycPolyInterp(const nodeVec_t &_t):n(_t.size()),t(_t),lambda(n) {
  init_lambda(); }

template <typename NODESCALAR>
template <typename SeqContainer> 
BarycPolyInterp<NODESCALAR>::BarycPolyInterp(const SeqContainer &v):n(v.size()),t(n),lambda(n) {
  idx_t ti = 0; for(auto tp: v) t(ti++) = tp;
  init_lambda();
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename NODESCALAR>
template <typename RESVEC,typename DATAVEC>
RESVEC BarycPolyInterp<NODESCALAR>::eval
  (const DATAVEC &y,const nodeVec_t &x) const {
  const idx_t N = x.size(); // No. of evaluation points
  RESVEC p(N); // Ouput vector
  // Compute quotient of weighted sums  of \Blue{$\frac{\lambda_i}{t - t_i}$}, effort \Blue{$O(n)$}
  for (unsigned i = 0; i < N; ++i) {
    nodeVec_t z = (x(i)*nodeVec_t::Ones(n) - t);

    // check if we want to evaluate at a node <-> avoid division by zero
    NODESCALAR* ptr = std::find(z.data(), z.data() + n, NODESCALAR(0)); // \Label{pv:1}
    if (ptr != z.data() + n) { // if ptr = z.data + n = z.end no zero was found
      p(i) = y(ptr - z.data()); // ptr - z.data gives the position of the zero
    }
    else {
      const nodeVec_t mu = lambda.cwiseQuotient(z);
      p(i) = (mu.cwiseProduct(y)).sum()/mu.sum();
    }
  } // end for
  return p;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename NODESCALAR>
void BarycPolyInterp<NODESCALAR>::init_lambda(void) {  
  // Precompute the weights \Blue{$\lambda_i$} with effort \Blue{$O(n^2)$}
  for (unsigned k = 0; k < n; ++k) {
    // little workaround: in \eigen cannot subtract a vector
    // from a scalar; multiply scalar by vector of ones
    lambda(k) = 1./((t(k)*nodeVec_t::Ones(k)-t.head(k)).prod()* 
            (t(k)*nodeVec_t::Ones(n-k-1)-t.tail(n-k-1)).prod());
  }
}
/* SAM_LISTING_END_3 */


