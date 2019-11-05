
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

#include "matplotlibcpp.h"


using namespace Eigen;
namespace plt = matplotlibcpp;


/* @brief Least squares polynomial fit to the data (x,y).
 * @param[in] x, y vectors of same size
 * @param[in] order, the degree of the fitted polynomial
 * @param[out] coeff, vector of size order+1 containing the monomial coefficients of
 *        the fitted polynomial p(x) = coeff[0] + coeff[1]*x + ... + coeff[order]*x^order.
 */
 /* SAM_LISTING_BEGIN_9 */
VectorXd polyfit(const VectorXd &x, const VectorXd &y, size_t order)
{
  Eigen::MatrixXd A(x.size(), order+1);

  assert(x.size() == y.size());
  assert(x.size() >= order+1);

  // Create matrix
  for (size_t i=0; i<x.size(); ++i) {
    for (size_t j=0; j<order+1; ++j) {
      A(i, j) = pow(x(i), j);
    }
  }

  // Solve for linear least squares fit
  VectorXd coeff = A.householderQr().solve(y);
  coeff.conservativeResize(order + 1);

  return coeff;
}
/* SAM_LISTING_END_9 */


/* @brief Compute the $L^\infty$ error of piecewise linear interpolation.
 * @param[in] f function with evaluation operator double operator()(double)const.
 * @param[in] t sorted mesh with t[0] = 0 and t[t.size()-1] = 1.
 * @param[out] $\|f - I_{\mathcal{M}}f\|_{L^\infty}$
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FUNCTION>
double pwlintpMaxError(FUNCTION &&f, const Eigen::VectorXd &t) {
  int n = t.size() - 1;
  int N = 1e5; // Sampling resolution for approximating L-infinty norm.
  double maxerr = 0.0; // L-infinity error so far.
  
  // TO DO (7-1.c): Approximate the maximum norm of (f - If) on [0,1],
  // where If is the piecewise linear interpolation of f on the mesh t.
  // START
  // We iterate over the sampling points x = 0, 1/N, 2/N, ... 1.
  // For each x we need to interpolate f in the correct mesh interval.
  // We are currently in the interval [tj, tj1] = [0,t(1)].
  int j = 0;
  double tj = t[j], tj1 = t[j+1];
  double fj = f(tj), fj1 = f(tj1);
  for(int i=0; i<=N; i++) {
    // We are currently in the mesh interval [tj, tj1].
    double x = (1.0*i)/N;
    // If x is outside the current mesh interval,
    // then find the interval containing x.
    while( x > tj1 ) { // x >= tj always satisfied.
      tj = tj1;
      fj = fj1;
      j++;
      tj1 = t[j+1];
      fj1 = f(tj1);
    }
    // Now we have tj <= x <= tj1.
    
    // Evaluate linear interpolant [tj,tj1] at x.
    double Ifx = fj + (fj1-fj)/(tj1-tj)*(x-tj);
    // Update L-infinity error on [0,x]
    double err = std::abs(f(x) - Ifx);
    if(err > maxerr){
      maxerr = err;
    }
  }
  // END
  return maxerr;
}
/* SAM_LISTING_END_0 */

/* @brief Estimates convergence rates of piecewise linear interpolation of $f(t)=t^\alpha$ using equidistant meshes on [0,1].
 * @param[in] alpha, a vector of positive values not equal to 0 or 1.
 * @param[out] Rates, a vector with Rates[j] = convergence rate of piecewise linear interpolation of t^alpha[j] using equidistant meshes on [0,1].
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd cvgrateEquidistantMesh(const Eigen::VectorXd &alpha) {
  int n_alphas = alpha.size();
  VectorXd Rates(n_alphas);
  
  // TO DO (7-1.d): Fill in the entries of Rates.
  // Hint: For each alpha(j), use polyfit() to estimate convergence
  // rates of the maximum norm (obtained by pwlintMaxError()) with
  // number of mesh intervals 32, 64, ..., 4096.
  // START
  int kmin=5, kmax = 12;
  int n_meshes = kmax-kmin+1;
  // n=2^kmin, 2^(kmin+1), ..., 2^kmax.
  ArrayXd N = Eigen::pow(2,ArrayXd::LinSpaced(n_meshes,kmin,kmax));
  
  // Errors(i,j) = L-infinity error of interpolation of t^alpha(j)
  // using a mesh M with #M = N(i)+1.
  MatrixXd Errors(n_meshes,n_alphas);
  for(int i=0; i<n_meshes; i++) {
    int n = N(i); // The mesh has #M = (n+1) = (N(i)+1) nodes.
    VectorXd mesh = VectorXd::LinSpaced(n+1,0.0,1.0);
    for(int j=0; j<n_alphas; j++) {
      // Compute L-infinity error for f_j(t) = t^alpha(j).
      double alphaj = alpha(j);
      auto fj = [alphaj](double t) { return std::pow(t,alphaj); };
      Errors(i,j) = pwlintpMaxError<std::function<double (double)>>(fj,mesh);
    }
  }
  
  for(int j=0; j<n_alphas; j++) {
    // Estimate convergence rate for f_j(t) = t^alpha(j) by the slope
    // of line fitted through the errors on a log-log scale.
    VectorXd coeff = polyfit(N.array().log(), Errors.col(j).array().log(), 1);
    Rates(j) = -coeff[1];
  }
  // END
  return Rates;
}
/* SAM_LISTING_END_1 */

/* @brief Tabulates and plots convergence rate of piecewise linear interpolation of t^alpha on equidistant meshes.
 */
/* SAM_LISTING_BEGIN_2 */
void testcvgEquidistantMesh(void) {
  plt::figure();
  // TO DO (7-1.e): Plot and tabulate the convergence rates of pw.
  // lin. intp. of &t^\alpha$ using equidistant meshes on [0,1],
  // for alpha = 0.1, 0.3, 0.5 ..., 2.9.
  // START
  int n_alphas = 15;
  VectorXd alpha = VectorXd::LinSpaced(n_alphas,0.1,2.9);
  
  VectorXd ConvRates = cvgrateEquidistantMesh(alpha);
  std::cout << std::setw(10) << "alpha" << std::setw(20) << "Convergence rate" << std::endl;
  for(int j=0; j<n_alphas; j++) {
    std::cout << std::setw(10) << alpha(j) << std::setw(20) << ConvRates(j) << std::endl;
  }
  
  plt::plot(alpha,ConvRates);
  plt::xlabel("alpha");
  plt::ylabel("conv. rate");
  plt::title("Pw. lin. intp. of f(t)=t^alpha on uniform meshes");
  // END
  plt::savefig("./cx_out/cvgrateEquidistant.png");
}
/* SAM_LISTING_BEGIN_2 */


/* @brief Estimates convergence rates of piecewise linear interpolation of $f(t)=t^\alpha$ on [0,1].
 * @param[in] alpha, a vector of positive values not equal to 0 or 1.
 * @param[in] beta, a vector of positive values.
 * @param[out] Rates, a matrix with Rates(k,j) = convergence rate of piecewise linear interpolation of t^alpha[j] using a beta[k]-graded mesh on [0,1].
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd cvgrateGradedMesh(const Eigen::VectorXd &alpha, const Eigen::VectorXd &beta) {
  int n_alphas = alpha.size();
  int n_betas = beta.size();
  MatrixXd Rates(n_betas,n_alphas);
  
  // TO DO (7-1.h): Fill in the entries of Rates.
  // START
  int kmin=5, kmax = 12;
  int n_meshes = kmax-kmin+1;
  // n=2^kmin, 2^(kmin+1), ..., 2^kmax.
  ArrayXd N = Eigen::pow(2,ArrayXd::LinSpaced(n_meshes,kmin,kmax));
  
  // Errors(i,k*n_alphas+j) = L-infinity error of interpolation of t^alpha(j)
  // using a mesh M with #M = N(i)+1 and grading parameter beta(k).
  MatrixXd Errors(n_meshes,n_alphas*n_betas);
  for(int i=0; i<n_meshes; i++) {
    int n = N(i); // The mesh has #M = (n+1) = (N(i)+1) nodes.
    VectorXd mesh = VectorXd::LinSpaced(n+1,0.0,1.0);
    for(int k=0; k<n_betas; k++) {
      // Transform the mesh using the grading parameter beta.
      VectorXd gradedMesh = pow(mesh.array(),beta(k)).matrix();
      for(int j=0; j<n_alphas; j++) {
        // Compute L-infinity error for f_j(t) = t^alpha(j).
        double alphaj = alpha(j);
        auto fj = [alphaj](double t) { return std::pow(t,alphaj); };
        int ind = k*n_alphas + j;
        Errors(i,ind) = pwlintpMaxError<std::function<double (double)>>(fj,gradedMesh);
      }
    }
  }
  
  for(int k=0; k<n_betas; k++) {
    for(int j=0; j<n_alphas; j++) {
      // Estimate convergence rate for f_j(t) = t^alpha(j) by the slope
      // of line fitted through the errors on a log-log scale.
      int ind = k*n_alphas + j;
      VectorXd coeff = polyfit(N.array().log(), Errors.col(ind).array().log(), 1);
      Rates(k,j) = -coeff[1];
    }
  }
  // END
  return Rates;
}
/* SAM_LISTING_END_3 */


/* @brief Prints convergence rates of piecewise linear interpolation of t^alpha on graded meshes.
 */
/* SAM_LISTING_BEGIN_4 */
void testcvgGradedMesh(void) {
  // TO DO (7-1.i): Print the table of convergence rates from
  // cvgrateGradedMesh() using alpha = 0.5, 0.15, ..., 1.95
  // and beta = 0.1, 0.2, ..., 2.0.
  // START
  int n_alphas = 20;
  VectorXd alpha = VectorXd::LinSpaced(n_alphas,0.05,1.95);
  int n_betas = 20;
  VectorXd beta = VectorXd::LinSpaced(n_betas,0.1,2.0);
  std::cout << "alpha = " << alpha.transpose()
            << "\nbeta = " << beta.transpose() << "\n\n";
            
  MatrixXd ConvRates = cvgrateGradedMesh(alpha,beta);
  std::cout << ConvRates << std::endl;
  
  plt::figure();
  // plt::contour(ConvRates); // Not implemented
  VectorXd Y = ConvRates.row(4);
  plt::plot(alpha,Y);
  Y = ConvRates.row(9);
  plt::plot(alpha,Y);
  Y = ConvRates.row(14);
  plt::plot(alpha,Y);
  Y = ConvRates.row(19);
  plt::plot(alpha,Y);
  plt::savefig("./cx_out/alphabeta.png");
  // END
}
/* SAM_LISTING_BEGIN_4 */
