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
            B(2*row,  0) = X(0,row);
            B(2*row,  1) = X(1,row);
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
                MatrixXd((B.transpose() * B).ldlt().solve(B.transpose() *
                                                 Map<const MatrixXd>(P.data(), 2*n, 1)
                                                 )).data(),
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
    MatrixXd Xstop(2, n);
    Xstop << // x-coords
             1, 3, 3, 1, -1, -3, -3, -1,
             // y-coords
             -3, -1, 1, 3, 3, 1, -1, -3;

    MatrixXd Xpriority(2, n);
    Xpriority << // x-coords
                 0, 3, 0, -3, 0, 2.5, 0, -2.5,
                 // y-coords
                 -3, 0, 3, 0, -2.5, 0, 2.5, 0;




    {
        MatrixXd P(2, n);
        P << -0.370407, 0.737504, 1.85296, 2.61666, 3.36683,
            2.61022, 3.74738, 0.370185, 1.11543, -1.10877,
            -0.745235, -2.97745, -2.23456, -1.86617, -1.12343,
            0.02186, -1.1121, -1.17265, -1.58382, -1.61935,
            -1.20599, -1.26846, -0.0132163, -0.0618339, 0.784893,
            0.765241, 1.60338, 0.825417, 0.808752, 0.0582148;

       std::cout << "****************** Set 1 ******************"
                  << std::endl;
       MatrixXd A;
       Shape s = identify(Xstop, Xpriority, P, A);
    }
    {
        MatrixXd P(2, n);
        P << 0.00384823, -0.594978, -0.453248, -0.537803, -0.253519,
                -0.137096, 0.0778092, 0.0626776, 0.196791, 0.272765,
                0.330656, 0.356591, 0.104327, 0.0293895, -0.187963,
                0.251845, -0.218366, -0.674227, -0.8617, -1.33682,
                -1.10144, -1.54197, -0.166215, -0.526587, 0.309162,
                0.262935, 1.17761, 0.857399, 0.719816, 0.562366;

       std::cout << "****************** Set 2 ******************"
                  << std::endl;
       MatrixXd A;
       Shape s = identify(Xstop, Xpriority, P, A);
    }
    {
        MatrixXd P(2, n);
        P << -0.867156, -0.11333, 0.252925, 0.509669, 1.1114,
                0.526225, 0.42909, -0.110557, -0.846383, -0.495474,
                -0.598265, 0.184648, 0.694457, 0.259439, -0.558467,
                -0.65055, -0.499286, -0.751277, -0.201265, 0.292899,
                0.346837, 0.807498, 0.253455, 0.240331, -0.121647,
                -0.380944, -0.477591, 0.144312, 0.554164, 0.135438;

        std::cout << "****************** Set 2 ******************"
                  << std::endl;
        MatrixXd A;
        Shape s = identify(Xstop, Xpriority, P, A);
    }


    return 0;
}
