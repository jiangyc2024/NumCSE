# include <vector>
# include "./trigpolyvalequid.hpp"
# include <figure/figure.hpp>

int main() {

  mgl::Figure fig_nonp;
  fig_nonp.fplot("x^2", "k|").label("Exact");

  std::vector<int> N = {5, 9, 13, 17, 25, 55};
  const int neval = 300;
  for (auto n : N) {
    VectorXd y = VectorXd::LinSpaced(n, 0, 1).array().pow(2).matrix();

    VectorXd q, x = VectorXd::LinSpaced(neval, 0, 1);
    q = trigpolyvalequid(y, neval);

    fig_nonp.plot(x, q).label("N = " + std::to_string(n));
  }
 
  fig_nonp.title("f(x) = x^2 (not 1-periodic)");
  fig_nonp.legend(0, 1);
  fig_nonp.save("non-periodic");

  mgl::Figure fig_p;
  fig_p.fplot("sin(2*pi*x)*cos(4*pi*x)", "k|").label("Exact");

  for (auto n : N) {
    Eigen::VectorXd t = Eigen::ArrayXd::LinSpaced(n, 0, 1);
	Eigen::VectorXd y = (2*M_PI*t.array()).sin()*(4*M_PI*t.array()).cos();

    VectorXd q, x = VectorXd::LinSpaced(neval, 0, 1);
    q = trigpolyvalequid(y, neval);

    fig_p.plot(x, q).label("N = " + std::to_string(n));
  }
 
  fig_p.title("f(x) = sin(2\\pi x) cos(4\\pi x)");
  fig_p.legend(0, 1);
  fig_p.save("periodic");


  return 0;
}
