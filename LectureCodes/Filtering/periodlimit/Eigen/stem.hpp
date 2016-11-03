# ifndef STEM_HPP
# define STEM_HPP

# include <string>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
using Eigen::VectorXd;

// stem plot, produces the same output as matlab's version:
// https://ch.mathworks.com/help/matlab/ref/stem.html
struct Stem {
  // title, xlabel, ylabel of plot and file to save it in
  std::string title, xlabel, ylabel, file;

  void plot(const VectorXd x, const VectorXd& y, const std::string& style) {
    if (x.size() != y.size()) {
      std::cout << "x and y vectors must have same sizes!\n";
      return;
    }

    mgl::Figure fig;
    fig.title(title); // set title
    fig.xlabel(xlabel); // set xlabel
    fig.ylabel(ylabel); // set ylabel
    for (std::size_t i = 0; i < x.size(); ++i) {
      // plot line between two y values for same x value
      std::vector<double> xvec = {x(i), x(i)},
                          yvec = {0, y(i)};
      fig.plot(xvec, yvec, style);
    }
    fig.save(file); // save plot
  }
};

# endif // STEM_HPP
