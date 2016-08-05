# include <iostream>
# include <cmath>
# include <Eigen/Dense>
# include "intpolyval.hpp"

template <class Function>
void adaptivepolyintp(const Function& f, const double a, const double b, 
                      const double tol, const unsigned N, 
                      Eigen::VectorXd& tRes, Eigen::VectorXd& errRes) {
  // get sampling points and evaluate f there
  Eigen::VectorXd sp = Eigen::VectorXd::LinSpaced(N, a, b),
                  fvals = sp.unaryExpr(f);
  double fmax = fvals.cwiseAbs().maxCoeff(); // approximate max |f(x)|
  std::vector<double> t { (a+b)/2. },     // set of interpolation nodes
                      y { f((a+b)/2.) },  // values in interpolation nodes
                      errors;               // save errors

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
      // if terminate save results in correct variables
      tRes = te; 
      Eigen::Map<Eigen::VectorXd> tmp(errors.data(), errors.size());
      errRes = tmp;  
      return;
    }

    // (iii) add this node to our set of nodes and save error
    errors.push_back(max);
    t.push_back(sp(idx));
    y.push_back(fvals(idx));
  } 
  std::cout << "Desired accuracy could not be reached.\n";
  tRes = sp; // return all sampling points
}





