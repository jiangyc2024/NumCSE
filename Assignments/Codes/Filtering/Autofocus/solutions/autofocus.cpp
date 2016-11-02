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
/* SAM_LISTING_BEGIN_1 */
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
/* SAM_LISTING_END_1 */

/*!
 * \brief plot_freq
 * \param focus
 */
/* SAM_LISTING_BEGIN_0 */
void plot_freq(double focus) {
    int a = 0;
    int b = 8000;
    auto clamp = [a,b] (double x) {
        return x < a ? a : x > b ? b : x;
    };

    MatrixXd D = fft2r(set_focus(focus))
            .cwiseAbs()
            .unaryExpr(clamp) / 8000;

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
        << ".png";
    gr.WritePNG(ss2.str().c_str());

}
/* SAM_LISTING_END_0 */

/*!
 * \brief high_frequency_content
 * \param M
 * \return
 */
/* SAM_LISTING_BEGIN_2 */
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
/* SAM_LISTING_END_2 */

/*!
 * \brief plotV
 */
/* SAM_LISTING_BEGIN_3 */
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
    fig.plot(x, y, "r+").label("$V(\mathbf{B}(f))$");
    fig.xlabel("$f$");
    fig.ylabel("$V(\mathbf{B}(f))$");
    fig.legend(0, 1);
    fig.save("focus_plot.eps");
    fig.save("focus_plot.png");
}
/* SAM_LISTING_END_3 */

/*!
 * \brief autofocus
 * \return
 */
/* SAM_LISTING_BEGIN_4 */
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
        double dV = computeV(f0+df) - computeV(f0);

        step = step / 2.;
        f0 = f0 + (dV > 0 ? 1 : -1) * step;
    }

    return f0;
}
/* SAM_LISTING_END_4 */

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
