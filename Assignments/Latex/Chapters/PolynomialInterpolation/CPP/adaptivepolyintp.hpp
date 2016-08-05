# include <iostream>
# include <cmath>
# include <Eigen/Dense>
# include "intpolyval.hpp"

template <class Function>
void adaptivepolyintp(const Function& f, const double a, const double b, const double tol, const unsigned N, Eigen::VectorXd& tRes) {
  Eigen::VectorXd sp = Eigen::VectorXd::LinSpaced(N, a, b), // sampling points
                  fvals = sp.unaryExpr(f); // evaluate f at sampling points
  double fmax = fvals.cwiseAbs().maxCoeff(); // approximate max |f(x)|
  std::vector<double> t { (a+b)/2. },     // set of interpolation nodes
                      y { f((a+b)/2.) };  // values in interpolation nodes

  for (int i = 0; i < N; ++i) {
    // (i) interpolate with current nodes 
    // need to convert std::vector to Eigen::VectorXd to use the function interpolate
    Eigen::Map<Eigen::VectorXd> te(t.data(), t.size());
    Eigen::Map<Eigen::VectorXd> ye(y.data(), y.size());
    Eigen::VectorXd spvals; 
    intpolyval(te, ye, sp, spvals);

    // (ii) find node where error is the largest
    Eigen::VectorXd err = (fvals - spvals).cwiseAbs();
    double max = 0; int idx = 0;
    for (int j = 0; j < err.size(); ++j) {
      if (err(j) > max) {
        max = err(j); idx = j;
      }
    }

    // check termination criteria
    if (max < tol*fmax) {
      tRes = te;
      return;
    }

    // (iii) add this node to our set of nodes
    t.push_back(sp(idx));
    y.push_back(fvals(idx));
  } 
  std::cout << "Desired accuracy could not be reached.\n";
  tRes = sp; // return all sampling points
}





