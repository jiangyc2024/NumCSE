#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "gradedmeshes.hpp"

using namespace Eigen;



int main() {
  std::cout << "\nEnter \"0\" to test all functions.\n"
            << "Enter \"1\" to only test pwlintpMaxError().\n"
            << "Enter \"2\" to only test cvplotEquidistantMesh().\n"
            << "Enter \"3\" to only test cvgrateEquidistantMesh().\n"
            << "Enter \"4\" to only test testcvgEquidistantMesh().\n"
            << "Enter \"5\" to only test cvgrateGradedMesh().\n"
            << "Enter \"6\" to only test testcvgGradedMesh().\n"
            << "\n";
  
  int ans=0;
  std::cin >> ans;
  
  if (ans==1 || ans==0) {
    VectorXd Mesh(11);
    Mesh << 0.0, 1e-6, 1e-4, 0.01, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0;
    std::cout << "Testing pwlintpMaxError() on a non-uniform mesh M = { "
              << Mesh.transpose() << " }\n";
    double err_sqrt = pwlintpMaxError( [](double x) { return std::sqrt(x); }, Mesh);
    double err_sin = pwlintpMaxError( [](double x) { return std::sin(10*x); }, Mesh);
    std::cout << "Error for f(x) = sqrt(x): " << err_sqrt << "\n";
    std::cout << "Error for f(x) = sin(10*x): " << err_sin << "\n\n";
  }
  std::cout << std::setprecision(3);
  if (ans==2 || ans==0) {
    VectorXd alpha(3);
    alpha << 0.3, 0.5, 0.7;
    std::cout << "Running cvgplotEquidistantMesh() for alpha = "
              << alpha.transpose() << " \n\n";
    cvgplotEquidistantMesh(alpha);
  }
  if (ans==3 || ans==0) {
    int n_alphas = 6;
    VectorXd alpha = VectorXd::LinSpaced(n_alphas,0.25,2.75);
    std::cout << "Testing cvgrateEquidistantMesh() for alpha = "
              << alpha.transpose() << " \n";
    VectorXd ConvRates = cvgrateEquidistantMesh(alpha);
    std::cout << "Convergence rates: "
              << ConvRates.transpose() << "\n\n";
  }
  if (ans==4 || ans==0) {
    std::cout << "Running testcvgEquidistantMesh()...\n\n";
    testcvgEquidistantMesh();
  }
  if (ans==5 || ans==0) {
    int n_alphas = 5;
    VectorXd alpha = VectorXd::LinSpaced(n_alphas,0.25,2.25);
    int n_betas = 4;
    VectorXd beta = VectorXd::LinSpaced(n_betas,0.5,2.0);
    std::cout << "Testing cvgrateGradedMesh() for\n\talpha = "
              << alpha.transpose()
              << ", and\n\tbeta = "
              << beta.transpose()
              << " \n";
    MatrixXd ConvRatesGraded = cvgrateGradedMesh(alpha,beta);
    std::cout << "Convergence rates:\n"
              << ConvRatesGraded << "\n\n";
  }
  if (ans==6 || ans==0) {
    std::cout << "Running testcvgGradedMesh()...\n\n";
    testcvgGradedMesh();
  }

  return 0;
    
}