#ifndef AUTOFOCUS_HPP
#define AUTOFOCUS_HPP

////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <fstream>
#include <vector>

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

// Contains FFT utilities
#include "fft2.hpp"

/**
 * \brief set_focus Given a double, retuns an image "taken" with
 * the double as "focusing parameter".
 * The focus parameter "f0" simulate the focal length chosen by
 * a digital camera. The function returns a $n \times m$ Matrix of doubles,
 * whose
 * values represent a grayscale image $n \times m$, with values
 * between 0 and 255, and where 0 represents black, and 255 represents
 * white.
 *
 * \param f0 The focal length parameter.
 * \return A MatrixXd, which contains the grey scale values of the image.
 */
Eigen::MatrixXd set_focus(double f0);

#include "_SUPER_SECRET_FILE.hpp"

/**
 * \brief Save differently blurred images.
 *
 * \param focus
 */
/* SAM_LISTING_BEGIN_1 */
void save_image(double focus) {
  // Create empty object
  PGMObject q;

  // TODO: (4-1.a) Read matrix of image generated
  // by "set\_focus" and save as an image in format ".pgm".
  // START

  // END
}
/* SAM_LISTING_END_1 */

/**
 * \brief Plot spectrum for different $f$.
 *
 * \param focus
 */
/* SAM_LISTING_BEGIN_0 */
void plot_freq(double focus) {
  constexpr unsigned int a = 0;
  constexpr unsigned int b = 8000;
  auto clamp = [](double x) { return x < a ? a : x > b ? b : x; };

  // TODO: (4-1.b) compute D containing the spectrum of set\_focus(focus).
  // "clamp" the data between 0 and 8000.
  // START
  Eigen::MatrixXd D;
  // END

  // Plot values of $\mathbf{X}$.
  // Axis labels
  std::vector<double> xticks(5);
  for (int i = 0; i < 5; ++i) {
    xticks[i] = i * (D.cols() - 1) / 4;
  }
  std::vector<double> yticks(5);
  for (int i = 0; i < 5; ++i) {
    yticks[i] = i * (D.rows() - 1) / 4;
  }
  std::vector<std::string> labels{"-1", "-0.5", "0", "0.5", "1"};

  plt::figure();
  plt::imshow(D, {{"cmap", "viridis"}, {"origin", "lower"}});
  plt::colorbar();
  plt::xticks(xticks, labels);
  plt::yticks(yticks, labels);
  std::stringstream ss;
  ss << "Spectrum with f = " << focus << ".";
  plt::title(ss.str().c_str());
  std::stringstream ss2;
  ss2 << "./cx_out/spectrum_focus=" << focus << ".png";
  plt::savefig(ss2.str().c_str());
}
/* SAM_LISTING_END_0 */

/**
 * \brief Compute $V(\mathbf{B}(f))$.
 *
 * \param M
 * \return V
 */
/* SAM_LISTING_BEGIN_2 */
double high_frequency_content(const Eigen::MatrixXd& M) {
  const unsigned int n = M.rows(), m = M.cols();
  double V = 0;
  // TODO: (4-1.d) compute $V(\mathbf{M})$.
  // START

  // END
  return V;
}
/* SAM_LISTING_END_2 */

/**
 * \brief Plot $V(\mathbf{B}(f))$.
 *
 */
/* SAM_LISTING_BEGIN_3 */
void plotV() {
  constexpr unsigned int N = 20;

  Eigen::VectorXd x(N), y(N);

  // TODO: (4-1.d) Plot $V(\mathbf{B}(f))$.
  // START

  // END
  plt::figure();
  plt::title("High frequency content.");
  plt::plot(x, y, "r+", {{"label", "$V(\\mathbf{B}(f))$"}});
  plt::xlabel("$f$");
  plt::ylabel("$V(\\mathbf{B}(f))$");
  plt::savefig("./cx_out/focus_plot.png");
}
/* SAM_LISTING_END_3 */

/**
 * \brief Find most focused image.
 *
 * \return
 */
/* SAM_LISTING_BEGIN_4 */
double autofocus() {
  // Minimum focus
  constexpr double min_focus = 0;
  // Maximum focus
  constexpr double max_focus = 5;
  // Min step
  constexpr double min_step = 0.05;
  // Starting guess
  double f0 = (max_focus - min_focus) / 2.;
  // Finite differences increment
  double df = min_step;
  // Starting step
  double step = max_focus / 2.;
  // Max number of iteration
  const unsigned int Niter =
      std::ceil(std::log2((max_focus - min_focus) / min_step));
  // TODO: (4-1.e) Use bisection method to find best focus.
  // START

  // END
  return f0;
}
/* SAM_LISTING_END_4 */

#endif
