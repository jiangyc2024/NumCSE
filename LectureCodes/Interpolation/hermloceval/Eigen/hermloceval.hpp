///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): J. Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Multiple point evaluation of Hermite polynomial
// \Blue{$y_1$}, \Blue{$y_2$}: data values
// \Blue{$c_1$}, \Blue{$c_2$}: slopes
VectorXd hermloceval(VectorXd t, double t1, double t2,
		     double y1, double y2,
		     double c1, double c2) {
  const double h = t2-t1,a1 = y2-y1, a2 = a1-h*c1, a3 = h*c2-a1-a2;
  t = ((t.array()-t1)/h).matrix();
  return (y1+(a1+(a2+a3*t.array())*(t.array()-1))*t.array()).matrix(); 
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
// Multiple point evaluation of Hermite polynomial
// \Blue{$y_1$}, \Blue{$y_2$}: data values
// \Blue{$c_1$}, \Blue{$c_2$}: slopes
void hermloceval(VectorXd t, double t1, double t2,
		 double y1, double y2,
		 double c1, double c2, VectorXd& p) {
  const double h = t2-t1,a1 = y2-y1, a2 = a1-h*c1, a3 = h*c2-a1-a2;
  t = ((t.array()-t1)/h).matrix();
  p = (y1+(a1+(a2+a3*t.array())*(t.array() - 1))*t.array() ).matrix(); 
}
/* SAM_LISTING_END_1 */
