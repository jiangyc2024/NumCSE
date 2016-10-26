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

#if SOLUTION
    for(unsigned int row = 0; row < n; ++row) {
        B(2*row,  0) = X(0,row);
        B(2*row,  1) = X(1,row);
        B(2*row+1,2) = X(0,row);
        B(2*row+1,3) = X(1,row);
    }
#else // TEMPLATE
    // TODO: Build system matrix $\mathnbf{B}$.
#endif // TEMPLATE

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
#if SOLUTION
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
#else // TEMPLATE
    // TODO: solve LSQ problem, return best linear approximation and residual
    return 0.;
#endif // TEMPLATE
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
#if SOLUTION
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
#else // TEMPALTE
    // TODO: Use residual do decide wether shape defines stop or priority road sign.
    return Stop;
#endif
}
/* SAM_LISTING_END_2 */

#if INTERNAL
// Internal: take points and "noisify" them.
MatrixXd generate_noisy_shape(const MatrixXd & X, double eps) {
    MatrixXd Xnoisy;
    Xnoisy.resizeLike(X);

    MatrixXd A = MatrixXd::Random(2,2);

    for(unsigned int i = 0; i < X.cols(); ++i) {
        Xnoisy.col(i) = A*X.col(i) + eps*VectorXd::Random(2);
    }

    return Xnoisy;

}

/*!
 * \brief transform Linearly transform points defined in $\mathbf{X}$.
 * Used to pull-back or push nforward transformed points.
 * \param X Points to transform.
 * \param A Linear transformation, $2 \ŧimes 2$ matrix.
 * \return Transformed points.
 */
MatrixXd transform(const MatrixXd & X,
                   const MatrixXd & A) {
    MatrixXd Xret;
    Xret.resizeLike(X);

    for(unsigned int i = 0; i < X.cols(); ++i) {
        Xret.col(i) = A*X.col(i);
    }

    return Xret;
}

/*!
 * \brief plot Save figure with points.
 * \param Xshape Base shape (if zero size, do not plot).
 * \param X Transformed shape.
 * \param title Title of the plot.
 * \param name Name of the file.
 * \param A Transform point with this matrix, only if specified.
 */
void plot(const MatrixXd & Xshape,
          const MatrixXd & X,
          const std::string & title,
          const std::string & name,
          const MatrixXd & A = MatrixXd::Identity(2,2)) {
    MatrixXd Xtransform = transform(X, A);

    mgl::Figure fig;
    fig.xlabel("x");
    fig.ylabel("y");
    fig.plot(Xtransform.row(0), Xtransform.row(1), " *b").label("Noisy points");
    if(Xshape.rows() == 2)
        fig.plot(Xshape.row(0), Xshape.row(1), "-r").label("Training shape");
    fig.title(title);
    fig.legend();
    fig.save(name + ".eps");
}
#endif // INTERNAL

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


#if INTERNAL
    MatrixXd _Xstop(2,n+1);
    _Xstop << Xstop, Xstop.leftCols<1>();
    MatrixXd _Xpriority(2, n+2);
    _Xpriority << Xpriority.leftCols<4>(), Xpriority.leftCols<1>(),
            Xpriority.rightCols<4>(), Xpriority.col(4);

    auto points = [&_Xstop, &_Xpriority] (Shape s) {
        return s == Stop ? _Xstop : _Xpriority;
    };
    {
        mgl::Figure fig;
        fig.xlabel("x");
        fig.ylabel("y");
        fig.plot(_Xstop.row(0), _Xstop.row(1));
        fig.plot(_Xstop.row(0), _Xstop.row(1), " b*");
        fig.title("Stop sign");
        fig.legend();
        fig.save("stop-pts.eps");
    }{
        mgl::Figure fig;
        fig.xlabel("x");
        fig.ylabel("y");
        MatrixXd Xpriority1(2, n/2+1);
        Xpriority1 << Xpriority.leftCols<4>(), Xpriority.leftCols<1>();
        MatrixXd Xpriority2(2, n/2+1);
        Xpriority2 <<  Xpriority.rightCols<1>(), Xpriority.rightCols<4>();
        fig.plot(Xpriority1.row(0), Xpriority1.row(1), "b");
        fig.plot(Xpriority2.row(0), Xpriority2.row(1), "b");
        fig.plot(Xpriority1.row(0), Xpriority1.row(1), " b*");
        fig.plot(Xpriority2.row(0), Xpriority2.row(1), " b*");
        fig.title("Priority road sign");
        fig.legend();
        fig.save("priority-pts.eps");
    }
#endif // INTERNAL

#if INTERNAL
    if(argc > 1) {
        srand(time(nullptr));

        MatrixXd X;
        double eps = 10e-3;
        if(argc > 2) {
            eps = std::stof(argv[2]);
        }

        if( std::stoi(argv[1]) == 1)
            X = generate_noisy_shape(Xstop, eps);
        else
            X = generate_noisy_shape(Xpriority, eps);

        for(unsigned int i = 0; i < n; ++i) {
            std::cout << X(0, i) << ", ";
        }
        for(unsigned int i = 0; i < n; ++i) {
            std::cout << X(1, i) << ", ";
        }

        return 0;
    }
#endif // INTERNAL

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
#if INTERNAL
        {
            mgl::Figure fig;
            fig.xlabel("x");
            fig.ylabel("y");
            fig.plot(P1.row(0), P1.row(1), " b*");
            fig.title("Example of points $\\mathbf{p}^i$");
            fig.legend();
            fig.save("example-p.eps");
        }

        plot(transform(points(s), A), P1,
             "Set 1: original points",
             "set-1-original");
        plot(points(s), P1,
             "Set 1: transformed back",
             "set-1-transformed", A.inverse());
#endif // INTERNAL
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
#if INTERNAL
        plot(transform(points(s), A), P2,
             "Set 2: original points",
             "set-2-original");
        plot(points(s), P2,
             "Set 2: transformed back",
             "set-2-transformed", A.inverse());
#endif // INTERNAL
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
#if INTERNAL
        plot(transform(points(s), A), P3,
             "Set 3: original points",
             "set-3-original");
        plot(points(s), P3,
             "Set 3: transformed back",
             "set-3-transformed", A.inverse());
#endif // INTERNAL
    }
    /* SAM_LISTING_END_3 */

    return 0;
}
