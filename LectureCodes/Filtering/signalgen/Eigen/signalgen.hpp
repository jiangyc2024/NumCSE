# include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::ArrayXd;

VectorXd signalgen() {
  ArrayXd t = ArrayXd::LinSpaced(64, 0, 63),
          x = (2*M_PI*t/64.).sin() + (7*2*M_PI*t/64.).sin();
  return (x + ArrayXd::Random(64)).matrix();
}
