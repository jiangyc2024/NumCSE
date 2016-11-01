//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <fstream>

#include <mgl2/mgl.h>

#include <figure/figure.hpp>

#include <Eigen/Dense>

// Contains PGMObject
#include "pgm.hpp"

// Contains definition of "set_focus"
#include "autofocus.hpp"

// Contains FFT utilities
#include "FFT/fft.hpp"

using namespace Eigen;

/*!
 * \brief save_image
 * \param focus
 */
void save_image(double focus) {
    // Create empty object
    PGMObject q;

    // Set data using function "set\_data"
    // Data obtained from "set\_focus"
    q.set_data(set_focus(focus));

    // Create and save file
    std::stringstream ss;
    ss << "image_focus="
       << (int) focus
       << ".pgm";
    std::ofstream file(ss.str());
    file << q;
}

/*!
 * \brief plot_freq
 * \param focus
 */
void plot_freq(double focus) {
    int a = 0;
    int b = 8000;
    auto clamp = [a,b] (double x) {
        return x < a ? a : x > b ? b : x;
    };

    MatrixXd D = fft2r(set_focus(focus))
            .cwiseAbs()
            .unaryExpr(clamp);

    // Plot values of $\mathbf{X}$.
    mglData Xd(D.cols(), D.rows(), D.data());

    mglGraph gr;
    gr.SetRange('c', 0, b);
    gr.Colorbar("bcwyr");
    std::stringstream ss;
    ss << "Specturm with f = "
        << focus
        << ".";
    gr.Title(ss.str().c_str());
    gr.Axis();
    gr.Tile(Xd, "bcwyr");
    std::stringstream ss2;
    ss2 << "spectrum_focus="
        << focus
//        << ".eps";
        << ".png";
//    gr.WriteEPS(ss2.str().c\_str());
    gr.WritePNG(ss2.str().c_str());

}

/*!
 * \brief high_frequency_content
 * \param M
 * \return
 */
double high_frequency_content(const MatrixXd & M) {

    int n = M.rows();
    int m = M.cols();

    double V = 0;
        for(unsigned int j = 0; j < M.cols(); ++j) {
            for(unsigned int i = 0; i < M.rows(); ++i) {
            double a = n/2 - std::abs(i - n/2);
            double b = m/2 - std::abs(j - m/2);
            V += (a*a + b*b) * M(i,j);
        }
    }

    return V;
}

/*!
 * \brief plotV
 */
void plotV() {

    unsigned int N = 100;

    VectorXd x(N), y(N);

    for(unsigned int i = 0; i < N; ++i) {
        double V = high_frequency_content(
                    fft2r(
                        set_focus(5. / (N-1) * i)
                        )
                    .cwiseAbs()
                    );
        x(i) = 5. / (N-1) * i;
        y(i) = V;
    }

    mgl::Figure fig;
    fig.title("High frequency content.");
//    fig.ranges(2, 9000, 1e-8, 1e3);
    fig.plot(x, y, " r+").label("V(\mathbf{B}(f))");
    fig.xlabel("f");
    fig.ylabel("V(\mathbf{B}(f))");
    fig.legend(0, 1);
    fig.save("focus_plot.eps");
    fig.save("focus_plot.png");
}

/*!
 * \brief autofocus
 * \return
 */
double autofocus() {

    // Max number of iteration
    unsigned int Niter = 6;

    // Maximum focus
    unsigned int max_focus = 5;
    // Starting guess
    double f0 = max_focus / 2.;
    // Finite differences increment
    double df = max_focus / 1e2;
    // Starting step
    double step = max_focus / 2.;
    // Returns $V(\mathbf{B}(f))$
    auto computeV = [] (double focus) {
        return high_frequency_content(
                    fft2r(
                        set_focus(focus)
                        ).cwiseAbs()
                    );
    };

    // Bisection method
    for(unsigned int i = 0; i < Niter; ++i) {
        double dV = computeV(f0+df) - computeV(f0);

        step = step / 2.;
        f0 = f0 + (dV > 0 ? 1 : -1) * step;
    }

    return f0;
}

// Comment to disable compilation of subproblem
#define SUBPROBLEM1
#define SUBPROBLEM2
#define SUBPROBLEM3
#define SUBPROBLEM4

int main() {

    //// SUBPROBLEM 1: save differently blurred images
#ifdef SUBPROBLEM1
    for(unsigned int i = 0; i <= 3; ++i) {
        save_image(i);
    }
#endif

    //// SUBPROBLEM 2: plot spectrum for different $f$
#ifdef SUBPROBLEM2
    for(unsigned int i = 0; i <= 3; ++i) {
        plot_freq(i);
    }
#endif

    //// SUBPROBLEM 3: plot V(\mathbf{B}(f))
#ifdef SUBPROBLEM3
    plotV();
#endif

    //// SUBPROBLEM 4: find most focused image
#ifdef SUBPROBLEM4
    std::cout << "Autofocus returns:"
              << autofocus()
              << std::endl;
#endif
}
