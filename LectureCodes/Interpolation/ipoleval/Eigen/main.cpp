# include "./ipolyeval.hpp"
# include <figure/figure.hpp>

int main() {
  /* testing ipoleval for sine for different degrees of interpolation
   * polynomials: 5, 10, 50 */
  VectorXd t0 = VectorXd::LinSpaced(5, 0, 6*M_PI),
           t1 = VectorXd::LinSpaced(10, 0, 6*M_PI),
           t2 = VectorXd::LinSpaced(50, 0, 6*M_PI),
           y0 = t0.array().sin().matrix(),
           y1 = t1.array().sin().matrix(),
           y2 = t2.array().sin().matrix(),
           x = VectorXd::LinSpaced(200, 0, 6*M_PI),
           v0, v1, v2;

  ipoleval(t0, y0, x, v0);
  ipoleval(t1, y1, x, v1);
  ipoleval(t2, y2, x, v2);

  mgl::Figure fig;
  fig.title("Polyfit + Polyval");
  fig.plot(x, v0, "r").label("deg(p) = 5");
  fig.plot(t0, y0, " ro");
  fig.plot(x, v1, "g").label("deg(p) = 10");
  fig.plot(t1, y1, " go");
  fig.plot(x, v2, "c").label("deg(p) = 50");
  fig.plot(t2, y2, " co");
  fig.fplot("sin(x)", "b|").label("Exact");
  fig.addlabel("Input data", " ko");
  fig.legend(0,0);
  fig.save("polyfit_polyval_intp");

  return 0;
}
