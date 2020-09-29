#ifndef SHAPE_IDENT_HPP
#define SHAPE_IDENT_HPP

#include <iostream>
#include <list>
#include <ctime>
#include <random>
#include <cmath>

#include <Eigen/Dense>

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
    // TODO: Build system matrix $\mathnbf{B}$.
    // START
    MatrixXd B = MatrixXd::Zero(2*n, 4);
    return B;
    // END
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
    // TODO: solve LSQ problem, return best linear approximation and residual
    // START
    return 0.;
    // END
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
    // TODO: Use residual do decide wether shape defines stop or priority road sign.
    // START
    return Stop;
    // END
}
/* SAM_LISTING_END_2 */

#endif
