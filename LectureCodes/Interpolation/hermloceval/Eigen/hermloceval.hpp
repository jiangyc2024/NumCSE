# include <Eigen/Dense>

using Eigen::VectorXd;

void hermloceval(VectorXd t, double t1, double t2, double y1, double y2, double c1, double c2, VectorXd& p) {
  // \Blue{$y_1$}, \Blue{$y_2$}: data values
  // \Blue{$c_1$}, \Blue{$c_2$}: slopes
  const double h = t2 - t1,
               a1 = y2 - y1,
               a2 = a1 - h*c1,
               a3 = h*c2 - a1 - a2;
  t = ( (t.array() - t1)/h ).matrix();
  p = ( y1+(a1+(a2+a3*t.array())*(t.array() - 1))*t.array() ).matrix(); 
}
