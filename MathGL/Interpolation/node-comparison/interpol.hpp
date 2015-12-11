# ifndef INTERPOL_HPP
# define INTERPOL_HPP

# include <Eigen/Dense>

class Interpol {
  public:
    Interpol(Eigen::VectorXd& t, Eigen::VectorXd& y);
    double operator()(double x);

  private:
    Eigen::VectorXd t_;
    Eigen::VectorXd y_;
    Eigen::VectorXd lambda_;
};

# endif
