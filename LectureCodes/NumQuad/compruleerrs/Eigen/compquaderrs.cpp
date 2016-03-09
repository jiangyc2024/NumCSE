# include <cmath>
# include "../../simpson/Eigen/simpson.hpp"
# include "../../trapezoidal/Eigen/trapezoidal.hpp"
# include <figure/figure.hpp>

void errs() {
  const double I_ex_f = 0.2*std::atan(5),
               I_ex_g = 2./3;

  auto f = [](double x) { return 1./(1 + 25*x*x); };
  auto g = [](double x) { return std::sqrt(x); };

  const unsigned steps = 200;
  std::vector<double> err_simp_f(steps), err_simp_g(steps),
                      err_trap_f(steps), err_trap_g(steps);

  for (unsigned i = 0; i < steps; ++i) {
    err_simp_f[i] = std::abs(simpson(f, 0, 1, steps) - I_ex_f);
    err_simp_g[i] = std::abs(simpson(g, 0, 1, steps) - I_ex_g);
    err_trap_f[i] = std::abs(trapezoidal(f, 0, 1, steps) - I_ex_ f);
    err_trap_g[i] = std::abs(trapezoidal(g, 0, 1, steps) - I_ex_g);
  }

  mgl::Figure fig_f, fig_g;
  fig_f.title("NumQuad for f(t) = \\frac{1}{1 + (5t)^2}");
  fig_f.xlabel("\\# of intervals");
  fig_f.ylabel("|quadrature error|");
  fig_f.plot(err_simp_f).label("Simpson");
  fig_f.plot(err_trap_f).label("Trapezoidal");
  fig_f.legend();
  fig_f.save("compruleerr1");

  fig_g.title("NumQuad for f(t) = \\sqrt{t}");
  fig_g.xlabel("\\# of intervals");
  fig_g.ylabel("|quadrature error|");
  fig_g.plot(err_simp_g).label("Simpson");
  fig_g.plot(err_trap_g).label("Trapezoidal");
  fig_g.legend();
  fig_g.save("compruleerr2");

  return;
}

int main () {
  errs();
  return 0;
}
