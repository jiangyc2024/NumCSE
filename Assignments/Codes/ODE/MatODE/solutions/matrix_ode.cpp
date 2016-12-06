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
    auto F = [] (const MatrixXd & M) {
        return -(M  - M.transpose())*M;
    };
    ode45<MatrixXd> O(F);

    // Set tolerances
    O.options.atol = 10e-10;
    O.options.rtol = 10e-8;

    // Return only matrix at $T$, (solution is vector
    // of pairs $(y(t_k), t_k)$ for each step k
    return O.solve(Y0, T).back().first;
}
/* SAM_LISTING_END_1 */

//! \brief Find if invariant is preserved after evolution with matode
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return true if invariant was preserved (up to round-off),
//! i.e. if norm was less than 10*eps
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const MatrixXd & M, double T) {
    MatrixXd N(3,3);

    N = matode(M, T);

    if( (N.transpose()*N-M.transpose()*M).norm() <
        10 * std::numeric_limits<double>::epsilon()* M.norm()) {
        return true;
    } else {
        return false;
    }
}
/* SAM_LISTING_END_2 */

//! \brief Implement ONE step of explicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_3 */
MatrixXd expeulstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
    return Y0 + h*A*Y0;
}
/* SAM_LISTING_END_3 */

//! \brief Implement ONE step of implicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_4 */
MatrixXd impeulstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
    return (MatrixXd::Identity(3,3) - h*A).partialPivLu().solve(Y0);
}
/* SAM_LISTING_END_4 */

//! \brief Implement ONE step of implicit midpoint ruler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_5 */
MatrixXd impstep(const MatrixXd & A, const MatrixXd & Y0, double h) {
    return (MatrixXd::Identity(3,3) - h*0.5*A)
            .partialPivLu()
            .solve(Y0+h*0.5*A*Y0);
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

    // Norm of Y'Y-I for 20 steps
    MatrixXd Mexpeul = Q, Mimpeul = Q, Mimp = Q;
    double h = 0.01;
    std::vector<int> sep = {8,15,15,15};
    std::cout << "Evolution of norm(Y_k'*Y_k - I) for three methods:" << std::endl;
    std::cout   << std::setw(sep[0]) << "step"
                << std::setw(sep[1]) << "exp. Eul"
                << std::setw(sep[2]) << "imp. Eul"
                << std::setw(sep[3]) << "IMP"
                << std::endl;
    std::cout   << std::setw(sep[0]) << "-1"
                << std::setw(sep[1]) << (Mexpeul.transpose()*Mexpeul - I).norm()
                << std::setw(sep[2]) << (Mimpeul.transpose()*Mimpeul - I).norm()
                << std::setw(sep[3]) << (Mimp.transpose()*Mimp - I).norm()
                << std::endl;
    for(unsigned int j = 0; j < 20; ++j) {
        Mexpeul = expeulstep(A, Mexpeul, h);
        Mimpeul = impeulstep(A, Mimpeul, h);
        Mimp = impstep(A, Mimp, h);

        std::cout   << std::setw(sep[0]) << j
                    << std::setw(sep[1]) << (Mexpeul.transpose()*Mexpeul - I).norm()
                    << std::setw(sep[2]) << (Mimpeul.transpose()*Mimpeul - I).norm()
                    << std::setw(sep[3]) << (Mimp.transpose()*Mimp - I).norm()
                    << std::endl;
    }
    /* SAM_LISTING_END_6 */

    std::cout << "SUBTASK d)" << std::endl;
    /* SAM_LISTING_BEGIN_7 */
    // Test implementation of ode45

    std::cout << "M = " << std::endl
              << M << std::endl;
    MatrixXd  N = matode(M, T);
    std::cout << "N = " << std::endl
              << N << std::endl;
    /* SAM_LISTING_END_7 */

    std::cout << "SUBTASK g)" << std::endl;
    /* SAM_LISTING_BEGIN_8 */
    // Test whether invariant was preserved or not

    bool is_invariant = checkinvariant(N, T);

    if( is_invariant ) {
        std::cout << "Invariant was preserved."
                  << std::endl;
    } else {
        std::cout << "Invariant was NOT preserved."
                  << std::endl;
    }
    /* SAM_LISTING_END_8 */

    return 0;
}
