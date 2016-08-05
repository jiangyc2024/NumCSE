# include <cmath>
# include <figure/figure.hpp>
# include "timer.h"
# include "./adaptivepolyintp.hpp"

int main() {
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
