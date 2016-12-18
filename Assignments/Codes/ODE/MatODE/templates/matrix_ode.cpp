#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/QR>

#include "ode45.hpp"

using namespace Eigen;

//! \file matrix_ode.cpp Solution for MatODE, involving ode45 and matrix ODEs

//! \brief Solve matrix IVP Y' = -(Y-Y')*Y using ode45 up to time T
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return Matrix of solution of IVP at t = T
/* SAM_LISTING_BEGIN_1 */
MatrixXd matode(const MatrixXd & Y0, double T) {
    // TODO: solve matrix ODE with ode45 class
    return Y0;
}
/* SAM_LISTING_END_1 */

//! \brief Find if invariant is preserved after evolution with matode
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return true if invariant was preserved (up to round-off),
//! i.e. if norm was less than 10*eps
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const MatrixXd & M, double T) {
    // TODO: test wether matrix satisfy invariant.
    return false;
}
/* SAM_LISTING_END_2 */

//! \brief Implement ONE step of explicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_3 */
MatrixXd expeulstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
    // TODO: implement one explicit Euler step
    return Y0;
}
/* SAM_LISTING_END_3 */

//! \brief Implement ONE step of implicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_4 */
MatrixXd impeulstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
    // TODO: implement one implicit Euler step
    return Y0;
}
/* SAM_LISTING_END_4 */

//! \brief Implement ONE step of implicit midpoint ruler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_5 */
MatrixXd impstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
    // TODO: implement one implicit midpoint step.
    return Y0;
}
/* SAM_LISTING_END_5 */

int main() {
    /* SAM_LISTING_BEGIN_6 */
    double T = 1;
    unsigned int n = 3;

    MatrixXd M(n,n);
    M << 8,1,6,3,5,7,4,9,2;

    std::cout << "SUBTASK c)" << std::endl;
    // Test preservation of orthogonality

    // Build Q
    HouseholderQR<MatrixXd> qr(M.rows(), M.cols());
    qr.compute(M);
    MatrixXd Q = qr.householderQ();

    // Build A
    MatrixXd A(n,n);
    A << 0, 1, 1, -1, 0, 1, -1, -1, 0;
    MatrixXd I = MatrixXd::Identity(n,n);

    // TODO
    /* SAM_LISTING_END_6 */

    /* SAM_LISTING_BEGIN_7 */
    std::cout << "Test implementation of ode45" << std::endl;
    // TODO
    /* SAM_LISTING_END_7 */

    /* SAM_LISTING_BEGIN_8 */
    std::cout << "Test whether invariant was preserved or not" << std::endl;
    // TODO
    /* SAM_LISTING_END_8 */
}
