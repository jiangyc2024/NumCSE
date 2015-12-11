# ifndef NATCSI_HPP
# define NATCSI_HPP

# include <Eigen/Dense>

class NatCSI {
  public:
    NatCSI(const Eigen::VectorXd& t, const Eigen::VectorXd& y);
    double operator() (double x) const;

  private:
    Eigen::MatrixXd t_;
    Eigen::VectorXd y_;
    Eigen::VectorXd c_;
};

# endif
