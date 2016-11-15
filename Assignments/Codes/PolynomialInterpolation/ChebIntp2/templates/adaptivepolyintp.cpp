# include <cmath>
# include <figure/figure.hpp>

template <class Function>
void adaptivepolyintp(const Function& f, const double a, const double b,
                      const double tol, const unsigned N,
                      Eigen::VectorXd& tRes, Eigen::VectorXd& errRes) {
  // get sampling points and evaluate f there
  Eigen::VectorXd sp = Eigen::VectorXd::LinSpaced(N, a, b),
                  fvals = sp.unaryExpr(f);
  double fmax = fvals.cwiseAbs().maxCoeff(); // approximate max |f(x)|
  std::vector<double> t { (a+b)/2. },     // set of interpolation nodes
                      y { f((a+b)/2.) },  // values in nodes
                      errors;               // save errors

  for (int i = 0; i < N; ++i) {
    // (i) interpolate with current nodes
    // need to convert std::vector to Eigen::VectorXd to use the function interpoyval
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


int main2() {
  // declare test functions
  auto f = [](double x) { return 1./(1 + std::exp(-x*x)); };
  auto g = [](double x) { return x > 0 ? 1 : -1; };
  auto h = [](double x) { return std::cos(x)*std::sin(2*x); };

  // test intervals
  const double af = -10,   bf = 10,
               ag = -1,    bg = 1,
               ah = -M_PI, bh = M_PI;

  // get interpolation nodes and print runtimes
  const unsigned N = 1000; // no. of sampling points
  const double tol = 1e-10; // tolerance
  Eigen::VectorXd tf, tg, th, // nodes
                  ef, eg, eh; // errors
  Timer fTimer, gTimer, hTimer;

  std::cout << "== Timing & Sizes =================================\n";
  fTimer.start(); adaptivepolyintp(f, af, bf, tol, N, tf, ef); fTimer.stop();
  std::cout << "Analytic function took " << fTimer.duration() << "s\n"
            << "Interval:     [" << af << ", " << bf << "]\n"
            << "No. of nodes: " << tf.size() << " (max = " << N << ")\n\n";

  gTimer.start(); adaptivepolyintp(g, ag, bg, tol, N, tg, eg); gTimer.stop();
  std::cout << "Step function took " << gTimer.duration() << "s\n"
            << "Interval: [" << ag << ", " << bg << "]\n"
            << "No. of nodes: " << tg.size() << " (max = " << N << ")\n\n";

  hTimer.start(); adaptivepolyintp(h, ah, bh, tol, N, th, eh); hTimer.stop();
  std::cout << "Trigonometric function took " << hTimer.duration() << "s\n"
            << "Interval: [" << ah << ", " << bh << "]\n"
            << "No. of nodes: " << th.size() << " (max = " << N << ")\n\n";

  // interpolate
  const Eigen::VectorXd xf = Eigen::VectorXd::LinSpaced(500, af, bf),
                        xg = Eigen::VectorXd::LinSpaced(500, ag, bg),
                        xh = Eigen::VectorXd::LinSpaced(500, ah, bh);
  Eigen::VectorXd If, Ig, Ih;
  intpolyval(tf, tf.unaryExpr(f), xf, If);
  intpolyval(tg, tg.unaryExpr(g), xg, Ig);
  intpolyval(th, th.unaryExpr(h), xh, Ih);

  // compute error
  const Eigen::VectorXd fx = xf.unaryExpr(f), // exact function values
                        gx = xg.unaryExpr(g),
                        hx = xh.unaryExpr(h);
  std::cout << "== Errors =========================================\n"
            << "Analytic function: \n"
            << "L^2:   " << (If - fx).lpNorm<2>() << "\n"
            << "L^inf: " << (If - fx).lpNorm<Eigen::Infinity>() << "\n\n"
            << "Step function: \n"
            << "L^2:   " << (Ig - gx).lpNorm<2>() << "\n"
            << "L^inf: " << (Ig - gx).lpNorm<Eigen::Infinity>() << "\n\n"
            << "Trigonometric function: \n"
            << "L^2:   " << (Ih - hx).lpNorm<2>() << "\n"
            << "L^inf: " << (Ih - hx).lpNorm<Eigen::Infinity>() << "\n\n";

  // plot
  std::cout << "== Plots ==========================================\n";
  mgl::Figure pf, pg, ph,
              pef, peg, peh;
  pf.title("f(t) = 1/(1 + e^{-t^2})");
  pf.plot(xf, xf.unaryExpr(f), "b").label("Function");
  pf.plot(xf, If, "r|").label("Interpolation");
  pf.plot(tf, tf.unaryExpr(f), " co").label("Data");
  pf.legend(0,0);
  pf.save("analytic");

  pef.title("f(t) = 1/(1 + e^{-t^2})");
  pef.setlog(true, true);
  pef.xlabel("No. of interpolation nodes");
  pef.ylabel("max_x |f(x) - I_T(x)|");
  pef.plot(ef, " r+").label("Error");
  pef.save("analyticError");

  std::cout << "Analytic function -> analytic.eps & analyticError.eps\n";

  pg.title("Step function");
  pg.plot(xg, xg.unaryExpr(g), "b").label("Function");
  pg.plot(xg, Ig, "r|").label("Interpolation");
  pg.plot(tg, tg.unaryExpr(g), " co").label("Data");
  pf.legend(0.5, 1);
  pg.save("step");

  peg.title("f(t) = 1/(1 + e^{-t^2})");
  peg.setlog(true, true);
  peg.xlabel("No. of interpolation nodes");
  peg.ylabel("max_x |f(x) - I_T(x)|");
  peg.plot(eg, " r+").label("Error"); // will throw an error as contains -nan
  peg.save("stepError");

  std::cout << "Step function -> step.eps & stepError.eps\n";

  ph.title("f(t) = cos(t)sin(2t)");
  ph.plot(xh, xh.unaryExpr(h), "b").label("Function");
  ph.plot(xh, Ih, "r|").label("Interpolation");
  ph.plot(th, th.unaryExpr(h), " co").label("Data");
  ph.legend(0,1);
  ph.save("trigonometric");

  peg.title("f(t) = cos(t)sin(2t)");
  peg.setlog(true, true);
  peg.xlabel("No. of interpolation nodes");
  peg.ylabel("max_x |f(x) - I_T(x)|");
  peg.plot(eh, " r+").label("Error");
  peg.save("trigonometricError");

  std::cout << "Trigonometric function -> trigonometric.eps & trigonometricError.eps\n\n";

  return 0;
}


int main() {
  // declare test functions
  auto f1 = [](double t) { return std::sin(std::exp(2*t)); };
  auto f2 = [](double t) { return std::sqrt(t)/(1 + 16*t*t); };

  // test interval
  const double a = 0, b = 1;

  // get interpolation nodes and print runtimes
  const unsigned N = 1000; // no. of sampling points
  const double tol = 1e-6; // tolerance
  Eigen::VectorXd tf1, tf2, // nodes
                  ef1, ef2; // errors

  adaptivepolyintp(f1, a, b, tol, N, tf1, ef1);
  adaptivepolyintp(f2, a, b, tol, N, tf2, ef2);

  // plot
  mgl::Figure fig;
  fig.title("Errors at each step");
  fig.setlog(false, true);
  fig.xlabel("No. of interpolation nodes");
  fig.ylabel("max |f(t) - I_Tf(t)|");
  fig.plot(ef1, "r").label("f_1(t) = sin(e^{2t})");
  fig.plot(ef2, "b").label("f_2(t) = \\sqrt{t}/(1 + 16t^2)");
  fig.legend();
  fig.save("cvg");

  return 0;
}
