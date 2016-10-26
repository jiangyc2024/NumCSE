//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
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
MatrixXd shape_ident_matrix(const MatrixXd & X) {
    assert(X.rows() == 2 && "X must have 2 rows!");
    unsigned n = X.cols();

    MatrixXd B = MatrixXd::Zero(2*n, 4);

    // TODO: Build system matrix $\mathnbf{B}$.

    return B;
}

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
double solve_lsq(const MatrixXd & X,
                 const MatrixXd & P,
                 MatrixXd & A) {
    assert(X.rows() == 2 && "X must have 2 rows!");
    assert(P.rows() == 2 && "P must have 2 rows!");
    assert(X.cols() == P.cols() && "P and X must have same size!");
    unsigned n = X.cols();
    // TODO: solve LSQ problem, return best linear approximation and residual
    return 0.;
}

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
Shape identify(const MatrixXd Xstop,
               const MatrixXd Xpriority,
               const MatrixXd & P,
               MatrixXd & A) {
    // TODO: Use residual do decide wether shape defines stop or priority road sign.
    return Stop;
}


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

        std::cout << "****************** Set 2 ******************"
                  << std::endl;
        MatrixXd A;
        Shape s = identify(Xstop, Xpriority, P3, A);
    }

    return 0;
}
