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

    // TODO: read matrix of image generated
    // by "set_focus" and same as an image in format ".pgm"
}

/*!
 * \brief plot_freq
 * \param focus
 */
void plot_freq(double focus) {
    int a = 0;
    int b = 8000;
    // TODO: compute D containing the
    // spectrum of set_focus(focus)
    // "clamp" the data between 0 and 8000
    MatrixXd D;

    // Plot values of $\mathbf{X}$.
    mglData Xd(D.cols(), D.rows(), D.data());
    mglGraph gr;
    gr.Colorbar("bcwyr");
    std::stringstream ss;
    ss << "Specturm with f = " << focus << ".";
    gr.Title(ss.str().c_str());
    gr.Axis(); gr.Tile(Xd, "bcwyr");
    std::stringstream ss2;
    ss2 << "spectrum_focus="  << focus << ".png";
    gr.WritePNG(ss2.str().c_str());
}

/*!
 * \brief high_frequency_content
 * \param M
 * \return
 */
double high_frequency_content(const MatrixXd & M) {
  int n = M.rows(),m = M.cols();
  double V = 0;
    // TODO: compute $V(\mathbf{M}).
  return V;
}

/*!
 * \brief plotV
 */
void plotV() {

    unsigned int N = 100;

    VectorXd x(N), y(N);

    // TODO: plot $V(\mathbf{B}(f))$

    mgl::Figure fig;
    fig.title("High frequency content.");
    fig.plot(x, y, "r+").label("$V(\\mathbf{B}(f))$");
    fig.xlabel("$f$");
    fig.ylabel("$V(\\mathbf{B}(f))$");
    fig.legend(0, 1);
    fig.save("focus_plot.eps");
    fig.save("focus_plot.png");
}

/*!
 * \brief autofocus
 * \return
 */
double autofocus() {
    // Minimum focus
    const double  min_focus = 0;
    // Maximum focus
    const double max_focus = 5;
    // Min step
    const double min_step = 0.05;
    // Starting guess
    double f0 = (max_focus - min_focus) / 2.;
    // Finite differences increment
    double df = min_step;
    // Starting step
    double step = max_focus / 2.;
    // Max number of iteration
    unsigned int Niter = std::log2(
                (max_focus - min_focus) / min_step
                );
    // TODO: use bisection method to find best focus

    return f0;
}

// Comment to disable compilation of subproblem
#define SUBPROBLEMa
#define SUBPROBLEMb
#define SUBPROBLEMc
#define SUBPROBLEMd

int main() {

    //// SUBPROBLEM a: save differently blurred images
#ifdef SUBPROBLEMa
    std::cout << "*** Subproblem a ***"
              << std::endl;
    for(unsigned int i = 0; i <= 3; ++i) {
        std::cout << "Saving image..."
                  << std::endl;
        save_image(i);
    }
#endif

    //// SUBPROBLEM b: plot spectrum for different $f$
#ifdef SUBPROBLEMb
    std::cout << "*** Subproblem b ***"
              << std::endl;
    for(unsigned int i = 0; i <= 3; ++i) {
        std::cout << "Saving plot..."
                  << std::endl;
        plot_freq(i);
    }
#endif

    //// SUBPROBLEM c: plot $V(\mathbf{B}(f))$
#ifdef SUBPROBLEMc
    std::cout << "*** Subproblem c ***"
              << std::endl;
    plotV();
#endif

    //// SUBPROBLEM d: find most focused image
#ifdef SUBPROBLEMd
    std::cout << "*** Subproblem d ***"
              << std::endl;
    std::cout << "Autofocus returns:"
              << autofocus()
              << std::endl;
#endif
}
