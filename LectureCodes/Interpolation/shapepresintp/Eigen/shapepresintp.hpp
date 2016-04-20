# include <iostream>
# include <cmath>
# include <string>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>

using Eigen::VectorXd;
using Eigen::ArrayXd;

// wait function (implemented below) to pause after each step
void wait(); 

void shapepresintp(const VectorXd& t, const VectorXd& y) {

  const long n = t.size() - 1;

  // ======= Step 1: choice of slopes ========
  // shape-faithful slopes (c) in the nodes using harmonic mean of data slopes
  // the final interpolant will be tangents to these segments
  std::cout << "STEP 1 - Shape-faithful slopes ...\n";

  const VectorXd h = t.tail(n) - t.head(n); // n = t.size() - 1
  const VectorXd delta = (y.tail(n) - y.head(n)).cwiseQuotient(h);

  std::cout << delta.transpose() << "\n";

  VectorXd c = VectorXd::Zero(n + 1);
  for (long j = 0; j < n - 1; ++j) {
    if (delta(j)*delta(j + 1) > 0) { // no change of sign
      c(j + 1) = 2./(1./delta(j) + 1./delta(j + 1));
    }
    else { // change of sign
      // as c is already zero we have nothing to do here
    }
  }
  // take care of first and last slope
  c(0) = 2*delta(0) - c(1);
  c(n) = 2*delta(n - 1) - c(n - 1);

  std::cout << c.transpose() << "\n";

  // plot segments indicating the slopes c(i)
  mgl::Figure segs;
  std::string seg_style = "c";
  std::vector<double> x(2), s(2);
  for (long j = 1; j < n; ++j) {
    x[0] = t(j) - 0.3*h(j - 1); x[1] = t(j) + 0.3*h(j);
    s[0] = y(j) - 0.3*h(j - 1)*c(j); s[1] = y(j) + 0.3*h(j)*c(j); 
    segs.plot(x, s, seg_style);
  }
  // plot first segment
  x[0] = t(0); x[1] = t(0) + 0.3*h(0);
  s[0] = y(0); s[1] = y(0) + 0.3*h(0)*c(0);
  segs.plot(x, s, seg_style);

  // plot last segment
  x[0] = t(n) - 0.3*h(n - 1); x[1] = t(n);
  s[0] = y(n) - 0.3*h(n - 1)*c(n); s[1] = y(n);
  segs.plot(x, s, seg_style);

  // plot data points
  segs.plot(t, y, " #m^");
  
  // save plot
  segs.grid(off);
  segs.save("segments");


  
  // ==== Step 2: choice of middle points ====
  // fix point \Blue{$p_j$} in \Blue{$t_j, t_{j+1}$} depending on the slopes \Blue{$c_j, c_{j+1}$}
  
  std::cout << "Step 2 - Middle points ...\n";

  VectorXd p = (t(0) - 1)*VectorXd::Ones(n);
  for (long j = 0; j < n; ++j) {
    if (c(j) != c(j + 1)) { // avoid division by zero
      p(j) = (y(j + 1) - y(j) + t(j)*c(j) - t(j + 1)*c(j + 1)) / (c(j) - c(j + 1));
    }
    // check and repair if \Blue{$p_j$} is outside its interval
    if (p(j) < t(j) || p(j) > t(j + 1)) { 
      p(j) = 0.5*(t(j) + t(j + 1));
    }
  }


  // ==== Step 3: auxiliary line spline =====
  // build the linear spline with nodes in:
  // - \Blue{$t_j$}
  // - the middle points between \Blue{$t_j$} and \Blue{$p_j$}
  // - the middle points between \Blue{$p_j$} and \Blue{$t_{j+1}$}
  // - \Blue{$t_{j+1}$}
  // and with slopes \Blue{$c_j$} in \Blue{$t_j$}, for every \Blue{$j$}
  
  std::cout << "Step 3 - Auxiliary line spline ...\n";

  x = std::vector<double>(4);
  s = std::vector<double>(4);

  mgl::Figure lin;

  for (long j = 0; j < n; ++j) {
    x[0] = t(j); s[0] = y(j);
    x[1] = (p(j) + t(j))/2.; s[1] = y(j) + c(j)*(p(j) - t(j))/2.;
    x[2] = (p(j) + t(j+1))/2.; s[2] = y(j+1) + c(j+1)*(p(j) - t(j+1))/2.;
    x[3] = t(j+1); s[3] = y(j+1);

    // plot linear spline
    lin.plot(x, s, "r");

    // plot additional evaluation points \Blue{$\frac{p_j + t_j}{2}$}, \Blue{$\frac{p_j + t_{j+1}}{2}$}
    std::vector<double> pt = {x[1], x[2]},
                        py = {s[1], s[2]};
    lin.plot(pt, py, " r^");
  }
  // plot data points
  lin.plot(t, y, " #r^").label("Data");
  lin.addlabel("Additional points", " r^");
  lin.addlabel("Linear splines", "r");
  // lin.legend(0, 0.3);
  lin.grid(false);
  lin.save("linearsplines");

  
  // ======= Step 4: quadratic spline ========
  // final quadratic shape preserving spline
  // quadratic polynomial in the intervals \Blue{$[t_j, p_j]$} and \Blue{$[p_j, t_{j+1}}$}
  // tangent in \Blue{$t_j$} and \Blue{$p_j$} to the linear spline of step 3
  
  std::cout << "Step 4 - Quadratic spline ...\n";

  // for every interval 2 quadratic interpolations
  // a, b, ya, yb = extremes and values in subinterval
  // w = value in middle point that gives the red slope
  mgl::Figure quad;
  for (long j = 0; j < n; ++j) {
    // handling the interval \Blue{$[t_j, p_j]$}
    double a = t(j), b = p(j), ya = y(j),
           w = y(j) + 0.5*c(j)*(p(j) - t(j)),
           yb = ( (t(j+1) - p(j))*(y(j) + 0.5*c(j)*(p(j) - t(j))) +
                  (p(j) - t(j))*(y(j+1) + 0.5*c(j+1)*(p(j) - t(j+1))) ) / (t(j+1) - t(j));

    ArrayXd eval = ArrayXd::LinSpaced(100, a, b),
            pb = (ya*(b - eval).pow(2) + 2*w*(eval - a)*(b - eval) + yb*(eval - a).pow(2))/std::pow(b - a, 2);

    quad.plot(eval.matrix(), pb.matrix(), "r");

    // now the same for interval \Blue{$p_j, t_{j+1}]$}
    a = b; b = t(j+1); ya = yb; yb = y(j+1);
    w = y(j+1) + 0.5*c(j+1)*(p(j) - t(j + 1));
    eval = ArrayXd::LinSpaced(100, a, b);
    pb = (ya*(b - eval).pow(2) + 2*w*(eval - a)*(b - eval) + yb*(eval - a).pow(2))/std::pow(b - a, 2);

    quad.plot(eval.matrix(), pb.matrix(), "r");
  }
  
  // plot data points
  quad.plot(t, y, " #mo");
  quad.grid(false);
  quad.save("quadraticspline");
  

}
