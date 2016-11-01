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

    // TODO: read matrix of image generated
    // by "set_focus" and same as an image in format ".pgm"
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
    // TODO: compute D containing the
    // spectrum of set_focus(focus)
    // "clamp" the data between 0 and 8000
    MatrixXd D;

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
//    gr.WriteEPS(ss2.str().c_str());
    gr.WritePNG(ss2.str().c_str());

}
/* SAM_LISTING_END_1 */

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
    // TODO: compute $V(\mathbf{M}).

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

    // TODO: plot $V(\mathbf{B}(f))$

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
    // TODO: use bisection method to find best focus

    return f0;
}
/* SAM_LISTING_END_4 */

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
