# include <iostream>
# include <Eigen/Dense>
# include <Eigen/Eigenvalues>
# include <figure/figure.hpp>
# include "circul.hpp"
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXcd; // complex for eigenvectors

// @brief Plot real and imaginary part of eigenvectors of n x n 
//        circulant matrix. 'ind' is just for the plot title
void circeig(const int n, const int ind) {
  const std::string title = "Circ matrix " + std::to_string(ind),
                    saveas = "circeig" + std::to_string(ind) + "ev";
  // get random circulant matrix and get eigenvectors
  MatrixXd C; circul(C, VectorXd::Random(n));
  Eigen::EigenSolver<MatrixXd> eig(C);
  MatrixXcd V = eig.eigenvectors();

  for (int j = 0; j < n; ++j) {
    MatrixXd data(n, 2); // prepare bar data
    // try: data << V.col(j).real(), V.col(j).imag();
    data.col(0) = V.col(j).real(); data.col(1) = V.col(j).imag();
    mgl::Figure fig;  fig.title(title);
    fig.bar(data, "br"); 
    fig.xlabel("vector component index");
    fig.ylabel("vector component value");
    fig.addlabel("real part", "b");
    fig.addlabel("imaginary part", "r");
    fig.legend(0,0);
    fig.save(saveas + std::to_string(j));
  }
}

int main() {
  circeig(8, 1);
  circeig(8, 2);
  return 0;
}
