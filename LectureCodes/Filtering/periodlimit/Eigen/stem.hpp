# ifndef STEM_HPP
# define STEM_HPP

# include <string>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
using Eigen::VectorXd;

struct Stem {
  std::string title, xlabel, ylabel, file;

  void plot(const VectorXd x, const VectorXd& y, const std::string& style) {
    if (x.size() != y.size()) {
      std::cout << "x and y vectors must have same sizes!\n";
      return;
    }

    mgl::Figure fig;
    fig.title(title);
    fig.xlabel(xlabel);
    fig.ylabel(ylabel);
    for (std::size_t i = 0; i < x.size(); ++i) {
      std::vector<double> xvec = {x(i), x(i)},
                          yvec = {0, y(i)};
      fig.plot(xvec, yvec, style);
    }
    fig.save(file);
  }
};

# endif // STEM_HPP
