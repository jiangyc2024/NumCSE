#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iomanip>
#include <iostream>

#include "spai.hpp"

int main() {
  std::cout
      << "Type 0 if you want to test spai(), 1 if you want to see the table of "
         "needed iterations for CG w/ or w/out preconditioner.\n";
  unsigned int input;
  std::cin >> input;
  switch (input) {
    case 0: {
      std::cout << "Test spai():\n";
      const unsigned int n = 40;
      Eigen::SparseMatrix<double> M(5, 5);
      M.coeffRef(3, 4) = 1;
      M.coeffRef(4, 4) = 2;
      M.coeffRef(1, 4) = 3;
      M.coeffRef(3, 3) = 4;
      M.coeffRef(3, 2) = 4;
      M.coeffRef(2, 3) = 4;
      M.coeffRef(2, 2) = 5;
      M.coeffRef(3, 1) = 6;
      M.coeffRef(0, 0) = 9;

      Eigen::SparseMatrix<double> N = spai(M);
      Eigen::SparseMatrix<double> I(5, 5);
      I.setIdentity();

      std::cout << "Error (n = 5): " << (I - M * N).norm() << std::endl;

      Eigen::SparseMatrix<double> M2(n * n, n * n);
      Eigen::SparseMatrix<double> I2(n, n);
      I2.setIdentity();

      Eigen::MatrixXd R = Eigen::MatrixXd::Random(n, n);
      M2 = kroneckerProduct(R, I2);

      Eigen::SparseMatrix<double> N2 = spai(M2);

      Eigen::SparseMatrix<double> Ibig(n * n, n * n);
      Ibig.setIdentity();

      std::cout << "Error (n = " << n * n << "): " << (Ibig - M2 * N2).norm()
                << std::endl;
      break;
    }
    case 1: {
      std::vector<std::pair<unsigned int, unsigned int>> tv = testSPAIPrecCG(6);
      std::cout << "Table for the number of iterations:\n";
      std::cout << std::setw(30) << "N" << std::setw(30)
                << "#it.'s w/out preconditioner" << std::setw(30)
                << "#it.'s w/ preconditioner\n";
      int i = 1;
      for (auto& pair : tv) {
        i *= 4;
        std::cout << std::setw(30) << i << std::setw(30) << std::get<0>(pair)
                  << std::setw(30) << std::get<1>(pair) << std::endl;
      }
      break;
    }
    default:
      break;
  }
}
