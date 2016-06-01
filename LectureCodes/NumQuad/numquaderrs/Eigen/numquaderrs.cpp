# include <cmath>
# include "numquad.hpp"
# include "figure.hpp"

void errs () {
  const unsigned N = 20;
  
  // \Blue{$f(x) = 1/(1 + 25x^2)$} on \Blue{$[0, 1]$}
  // \Blue{$g(x) = sqrt(x)$} on \Blue{$[0,1]$}
  auto f = [](double x){ return 1./(1 + 25*x*x); };
  auto g = [](double x){ return std::sqrt(x); };
  const double I_ex_f = std::atan(5)/5,
               I_ex_g = 2./3;
  std::vector<double> equi_f = numquad(f, 0, 1, N, "equidistant"),
                      cheb_f = numquad(f, 0, 1, N, "chebychev"),
                      gauss_f = numquad(f, 0, 1, N, "gauss"),
                      equi_g = numquad(g, 0, 1, N, "equidistant"),
                      cheb_g = numquad(g, 0, 1, N, "chebychev"),
                      gauss_g = numquad(g, 0, 1, N, "gauss");

  // compute errors
  std::vector<double> err_equi_f(N), err_equi_g(N),
                      err_cheb_f(N), err_cheb_g(N),
                      err_gauss_f(N), err_gauss_g(N);

  for (unsigned i = 0; i < N; ++i) {
    err_equi_f[i] = std::abs(equi_f[i] - I_ex_f);
    err_cheb_f[i] = std::abs(cheb_f[i] - I_ex_f);
    err_gauss_f[i] = std::abs(gauss_f[i] - I_ex_f);
    err_equi_g[i] = std::abs(equi_g[i] - I_ex_g);
    err_cheb_g[i] = std::abs(cheb_g[i] - I_ex_g);
    err_gauss_g[i] = std::abs(gauss_g[i] - I_ex_g);
  }

  // convergence plot for f
  mgl::Figure fig_f;
  fig_f.title("f(x) = \\frac{1}{1 + \\ (5t)^2}");
  fig_f.setlog(true, true);
  fig_f.plot(err_equi_f, "+b").label("Equidistant");
  fig_f.plot(err_cheb_f, "+r").label("Chebychev");
  fig_f.plot(err_gauss_f, "+g").label("Gauss");
  fig_f.legend();
  fig_f.save("numquaderr1");


  // convergence plot for g
  mgl::Figure fig_g;
  fig_g.title("f(x) = \\sqrt{x}");
  fig_g.setlog(true, true);
  fig_g.plot(err_equi_g, "+b").label("Equidistant");
  fig_g.plot(err_cheb_g, "+r").label("Chebychev");
  fig_g.plot(err_gauss_g, "+g").label("Gauss");
  fig_g.legend();
  fig_g.save("numquaderr2");

  return;
}

int main () {
  errs();
  return 0;
}
