#include <iostream>
#include <list>
#include <ctime>
#include <random>
#include <cmath>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

/*!
 * \brief shape_ident_matrix Build matrix $\mathbf{B}$.
 * Build the overdetermined system matrix arising from the
 * point $\mathbf{x}^i$.
 * \param[in] X A $2 \times n$ matrix with model points.
 * \return The system matrix $\mathbf{B}$.
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXd shape_ident_matrix(const MatrixXd & X) {
    assert(X.rows() == 2 && "X must have 2 rows!");
    unsigned n = X.cols();

    MatrixXd B = MatrixXd::Zero(2*n, 4);

    for(unsigned int row = 0; row < n; ++row) {
        // Odd row, first two columns
        B(2*row,  0) = X(0,row);
        B(2*row,  1) = X(1,row);
        // Even row, last two columns
        B(2*row+1,2) = X(0,row);
        B(2*row+1,3) = X(1,row);
    }

    return B;
}
/* SAM_LISTING_END_0 */

/*!
 * \brief solve_lsq Solve least square problem.
 * Build the overdetermined system, find best solution in least square sense.
 * Then return norm of residual and store in $\mathbf{A}$ the $2 \times 2$ linear
 * transformation that is the LSQ solution of the system.
 * \param[in] X Model points, a $2 \times n$ Matrix.
 * \param[in] P Real points, a $2 \times n$ Matrix.
 * \param[out] A Best $2 \times 2$ linear trtansfomration (in LSQ sense).
 * \return Norm of residual.
 */
/* SAM_LISTING_BEGIN_1 */
double solve_lsq(const MatrixXd & X,
                 const MatrixXd & P,
                 MatrixXd & A) {
    assert(X.rows() == 2 && "X must have 2 rows!");
    assert(P.rows() == 2 && "P must have 2 rows!");
    assert(X.cols() == P.cols() && "P and X must have same size!");
    unsigned n = X.cols();
    // Build system matrix
    MatrixXd B = shape_ident_matrix(X);

    // Solve LSQ system using normal equation
    // We need to do some reshaping in order to properly set up the system
    A = Map<MatrixXd>(
                MatrixXd(
                    // Solve LSQ problem using normal equation
                    (B.transpose() * B).ldlt()
                                       .solve(B.transpose() *
                                              // Ned to vectorize matrix
                                              Map<const MatrixXd>(P.data(), 2*n, 1)
                                              )
                    // Pass an array to "Map"
                    ).data(),
                // Need to reshape vector to 2x2 matrix
                2, 2).transpose();

    // Residual: must reshape P
    return (Map<const MatrixXd>(P.data(), 2, n) - A*X).norm();
}
/* SAM_LISTING_END_1 */

//! Enum used for classification.
enum Shape { Stop, Priority };

/*!
 * \brief identify Choose if points represent a stop or priority sign.
 * Use LSQ to find best linear transformation for both models, then
 * decide which one fits the best.
 * \param[in] Xstop "Model" points for stop sign.
 * \param[in] Xpriority "Model" points for priority sign.
 * \param[in] P Real points (e.g. photographed).
 * \param[in] A Return the $2 \times 2$ linear transformation matrix.
 * \return The decided shape.
 */
/* SAM_LISTING_BEGIN_2 */
Shape identify(const MatrixXd Xstop,
               const MatrixXd Xpriority,
               const MatrixXd & P,
               MatrixXd & A) {
    // Best linear transformation for both "model" shape
    MatrixXd Astop, Apriority;

    // Residuals and linear transform for both "models"
    double res_stop = solve_lsq(Xstop, P, Astop);
    double res_priority = solve_lsq(Xpriority, P, Apriority);

    std::cout << "'Stop' residual norm: " << res_stop << std::endl
              << "'Prioriy' residual norm: " << res_priority << std::endl;

    // If residual with stop Model is bigger than the one with priority sign,
    // probably it it is priority sign, otherwise is a stop sign.
    if(res_stop <= res_priority) {
        std::cout << "Points appear to define a stop sign!" << std::endl;
        A = Astop;
        return Stop;
    }
    else {
        std::cout << "Points appear to define a priority road sign!" << std::endl;
        A = Apriority;
        return Priority;
    }
}
/* SAM_LISTING_END_2 */


int main(int argc, char **argv) {
    const unsigned int n = 8;
    // X_{stop}
    MatrixXd Xstop(2, n);
    Xstop << // x-coords
             1, 3, 3, 1, -1, -3, -3, -1,
            // y-coords
            -3, -1, 1, 3, 3, 1, -1, -3;

    // X_{priority}
    MatrixXd Xpriority(2, n);
    Xpriority << // x-coords
                 0, 3, 0, -3, 0, 2.5, 0, -2.5,
            // y-coords
            -3, 0, 3, 0, -2.5, 0, 2.5, 0;




    /* SAM_LISTING_BEGIN_3 */
    {
        MatrixXd P1(2, n);
        P1 << 0.23657, 1.35369, -0.13624, -1.33702,
                0.0989619, 0.993235, -0.0735973, -1.11657,
                -2.76114, -2.60103, 2.90403, 2.66831,
                -2.44302, -2.04656, 2.31922, 2.20296;

        std::cout << "****************** Set 1 ******************"
                  << std::endl;
        MatrixXd A;
        Shape s = identify(Xstop, Xpriority, P1, A);
    }
    {
        MatrixXd P2(2, n);
        P2 << -1.12783, -1.75868, -1.40935, -0.0664574,
                1.09654, 1.75873, 1.5195, -0.0607661,
                1.72169, 0.344036, -0.889686, -1.87847,
                -1.57535, -0.41511, 0.834371, 1.88514;

        std::cout << "****************** Set 2 ******************"
                  << std::endl;
        MatrixXd A;
        Shape s = identify(Xstop, Xpriority, P2, A);
    }
    {
        MatrixXd P3(2, n);
        P3 << -1.23988, -0.731643, 0.00492048, 1.08039,
                1.34128, 0.670982, -0.101797, -1.02859,
                1.53076, 2.02881, 1.36163, -0.340912,
                -1.47697, -1.99975, -1.47947, 0.374859;

        std::cout << "****************** Set 23 ******************"
                  << std::endl;
        MatrixXd A;
        Shape s = identify(Xstop, Xpriority, P3, A);
    }
    /* SAM_LISTING_END_3 */

    return 0;
}
