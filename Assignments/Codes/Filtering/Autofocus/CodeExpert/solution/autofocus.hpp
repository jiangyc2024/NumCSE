//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <fstream>
#include <Eigen/Dense>

#include <mgl2/mgl.h> // TODO: deprecate
#include <figure/figure.hpp>

// Contains PGMObject
#include "pgm.hpp"

// Contains definition of "set_focus"
#include "autofocus.hpp"

// Contains FFT utilities
#include "FFT/fft.hpp"

using namespace Eigen;

/*!
 * \brief Save differently blurred images.
 * \param focus
 */
void save_image(double focus) {
    // Create empty object
    PGMObject q;

    // TO DO: (a) Read matrix of image generated
    // by "set_focus" and save as an image in format ".pgm".
    // START
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
    // END
}

/*!
 * \brief Plot spectrum for different $f$.
 * \param focus
 */
void plot_freq(double focus) {
    int a = 0;
    int b = 8000;
    auto clamp = [a,b] (double x) {
        return x < a ? a : x > b ? b : x;
    };

    // TO DO: (b) compute D containing the spectrum of set_focus(focus).
    // "clamp" the data between 0 and 8000.
    // START
    MatrixXd D = fft2r(set_focus(focus))
            .cwiseAbs()
            .unaryExpr(clamp) / b;
    // END

    // Plot values of $\mathbf{X}$.
    mglData Xd(D.cols(), D.rows(), D.data());
    mglGraph gr;
    gr.Colorbar("bcwyr");
    std::stringstream ss;
    ss << "Spectrumm with f = "
       << focus
       << ".";
    gr.Title(ss.str().c_str());
    gr.Axis(); gr.Tile(Xd, "bcwyr");
    std::stringstream ss2;
    ss2 << "spectrum_focus="
        << focus
        << ".png";
    gr.WritePNG(ss2.str().c_str());
}

/*!
 * \brief Compute $V(\mathbf{B}(f))$.
 * \param M
 * \return V
 */
double high_frequency_content(const MatrixXd & M) {
    int n = M.rows(),m = M.cols();
    double V = 0;
    // TO DO: compute $V(\mathbf{M}).
    // START
    for(unsigned int i = 0; i < M.rows(); ++i) {
        for(unsigned int j = 0; j < M.cols(); ++j) {
            double a = n/2. - std::abs(i - n/2.);
            double b = m/2. - std::abs(j - m/2.);
            V += (a*a + b*b) * M(i,j) * M(i,j);
        }
    }
    // END
    return V;
}

/*!
 * \brief Plot $V(\mathbf{B}(f))$.
 */
void plotV() {

    unsigned int N = 100;

    VectorXd x(N), y(N);

    // TO DO: (c) Plot $V(\mathbf{B}(f))$.
    // START
    for(unsigned int i = 0; i < N; ++i) {
        double V = high_frequency_content(
                    // Find 2D spectrum of matrix $\mathbf{B}(t)$
                    fft2r(
                        // Evaluate set\_focus at equidistant points
                        set_focus(5. / (N-1) * i)
                        )
                    .cwiseAbs()
                    );
        x(i) = 5. / (N-1) * i;
        y(i) = V;
    }
    
    // END
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
 * \brief Find most focused image.
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
    unsigned int Niter = std::ceil(std::log2(
                (max_focus - min_focus) / min_step
                ));
    // TO DO: (d) Use bisection method to find best focus.
    // START
    // Returns $V(B(f))$
    auto computeV = [] (double focus) {
        return high_frequency_content(
                    fft2r(
                        set_focus(focus)
                        ).cwiseAbs()
                    );
    };

    // Bisection method
    for(unsigned int i = 0; i < Niter; ++i) {
        double dV = computeV(f0+df) - computeV(f0-df);

        step = step / 2.;
        f0 = f0 + (dV > 0 ? 1 : -1) * step;
    }
    // END
    return f0;
}