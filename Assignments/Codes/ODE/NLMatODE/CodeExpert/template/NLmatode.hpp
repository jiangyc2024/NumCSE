#ifndef NLMATODE_HPP
#define NLMATODE_HPP

#include "ode45.hpp"

using namespace Eigen;


//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_1 */
MatrixXd matode(const MatrixXd &Y0, double T) {
  // TO DO (12-5.a): use the ode45 class to find an approximation
  // of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
  
  // START
  // Define the RHS
  auto F = [] (const MatrixXd & M) {
      return -(M  - M.transpose())*M;
  };
  ode45<MatrixXd> O(F);

  // Set tolerances
  O.options.atol = 10e-10;
  O.options.rtol = 10e-8;

  // Return only matrix at $T$, (solution is vector
  // of pairs $(y(t_k), t_k)$ for each step k
  MatrixXd YT = O.solve(Y0, T).back().first;
  // END
  return YT;
}
/* SAM_LISTING_END_1 */


//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const MatrixXd &M, double T) {
  // TO DO (12-5.c): check if $Y'*Y$ is preserved at the time $T$ by matode.
  // START
  MatrixXd N = matode(M, T);

  if( (N.transpose()*N-M.transpose()*M).norm() <
      10 * std::numeric_limits<double>::epsilon()* M.norm()) {
      return true;
  } else {
      return false;
  }
  // END
}
/* SAM_LISTING_END_2 */

#endif