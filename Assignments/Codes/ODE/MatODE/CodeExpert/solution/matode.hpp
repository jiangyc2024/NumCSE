#ifndef MATODE_HPP
#define MATODE_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/QR>


using namespace Eigen;


//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
/* SAM_LISTING_BEGIN_3 */
MatrixXd eeulstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
  // TO DO: Implement ONE step of explicit euler applied to Y0, 
  // for the ODE Y' = A*Y
  // START
  return Y0 + h*A*Y0;
  // END
}
/* SAM_LISTING_END_3 */


//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
/* SAM_LISTING_BEGIN_4 */
MatrixXd ieulstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
  // TO DO: Implement ONE step of implicit euler applied to Y0, 
  // for the ODE Y' = A*Y
  // START
  return (MatrixXd::Identity(3,3) - h*A).partialPivLu().solve(Y0);
  // END
}
/* SAM_LISTING_END_4 */


//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
/* SAM_LISTING_BEGIN_5 */
MatrixXd impstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
  // TO DO: Implement ONE step of implicit midpoint rule applied to Y0, 
  // for the ODE Y' = A*Y
  // START
  return (MatrixXd::Identity(3,3) - h*0.5*A)
          .partialPivLu()
          .solve(Y0+h*0.5*A*Y0);
  // END
}
/* SAM_LISTING_END_5 */


/* SAM_LISTING_BEGIN_6 */
std::tuple<double,double,double> checkOrthogonality(void) {
  // TO DO (12-1.c): compute and tabulate the Frobenius norms of Y_k'*Y_k - I
  // for 20 steps of eeulstep, ieulstep and impstep. 
  // Return the values corresponding to the last step.
  // START
  unsigned int n = 3;
  double h = 0.01;
  MatrixXd M(n,n);
  M << 8, 1, 6, 3, 5, 7, 9, 9, 2;

  // Build Q
  HouseholderQR<MatrixXd> qr(n, n);
  qr.compute(M);
  MatrixXd Q = qr.householderQ();

  // Build A
  MatrixXd A(n,n);
  A << 0, 1, 1, -1, 0, 1, -1, -1, 0;
  MatrixXd I = MatrixXd::Identity(n,n);

  MatrixXd Meeul = Q, Mieul = Q, Mimp = Q;
  
  std::vector<int> sep = {5,15};
  std::cout << std::setw(sep[0]) << "step"
            << std::setw(sep[1]) << "exp. Eul"
            << std::setw(sep[1]) << "imp. Eul"
            << std::setw(sep[1]) << "Mid-Pt"
            << std::endl;
  // Norm of Y'Y-I for initial value
  std::cout << std::setw(sep[0]) << "-1"
            << std::setw(sep[1]) << (Meeul.transpose()*Meeul - I).norm() 
            << std::setw(sep[1]) << (Mieul.transpose()*Mieul - I).norm()
            << std::setw(sep[1]) << (Mimp.transpose()*Mimp - I).norm() 
            << std::endl;
            
  std::vector<double> norms(3);
  // Norm of Y'Y-I for 20 steps
  for(unsigned int j = 0; j < 20; ++j) {
    Meeul = eeulstep(A, Meeul, h);
    Mieul = ieulstep(A, Mieul, h);
    Mimp = impstep(A, Mimp, h);
    
    norms[0] = (Meeul.transpose()*Meeul - I).norm();
    norms[1] = (Mieul.transpose()*Mieul - I).norm();
    norms[2] = (Mimp.transpose()*Mimp - I).norm(); 
    
    std::cout  << std::setw(sep[0]) << j
              << std::setw(sep[1]) << norms[0]
              << std::setw(sep[1]) << norms[1]
              << std::setw(sep[1]) << norms[2]
              << std::endl;
  }
  return std::make_tuple( norms[0], norms[1], norms[2] );
  // END
}
/* SAM_LISTING_END_6 */

#endif