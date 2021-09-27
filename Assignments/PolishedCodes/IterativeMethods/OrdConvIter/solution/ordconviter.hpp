#ifndef ORDCONVITER_HPP
#define ORDCONVITER_HPP

#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <vector>

#include "matplotlibcpp.h"
#include "meshgrid.hpp"

namespace plt = matplotlibcpp;

void kminplot() {
  constexpr double p = 1.5;
  constexpr double C = 2.;
  const double eps_max = std::pow(C, 1. / (1. - p));
  constexpr unsigned int ngp = 100;  // number of grid points

  plt::figure();

  // TODO: (9-5.c) Plot kmin using matplotlibcpp. To that end, start by creating
  // a mesh using the provided meshgrid function in meshgrid.hpp. START
  Eigen::VectorXd eps_lin = Eigen::VectorXd::LinSpaced(ngp, 0, eps_max);
  Eigen::VectorXd tau_lin = eps_lin;
  // we need an open interval, i.e. remove first and last element
  Eigen::VectorXd eps_lin_seg = eps_lin.segment(1, eps_lin.size() - 2);
  Eigen::VectorXd tau_lin_seg = tau_lin.segment(1, tau_lin.size() - 2);
  Eigen::MatrixXd eps_msh, tau_msh;
  meshgrid(eps_lin_seg, tau_lin_seg, eps_msh, tau_msh);

  // apply formula to mesh
  Eigen::MatrixXd kmin;
  Eigen::MatrixXd num = tau_msh.array().log();
  num += (1. / (p - 1.)) * std::log(C) *
         Eigen::MatrixXd::Ones(tau_msh.rows(), tau_msh.cols());
  Eigen::MatrixXd den = (eps_msh * std::pow(C, 1. / (p - 1.))).array().log();
  kmin = (num.cwiseQuotient(den)).array().log() / std::log(p);

  // Consider only gridpoints where eps is larger than tau
  for (unsigned int i = 0; i < kmin.rows(); ++i) {
    for (unsigned int j = 0; j < kmin.cols(); ++j) {
      kmin(i, j) = std::ceil(kmin(i, j));
      if (i > j) {
        kmin(i, j) = 0;
      }
    }
  }

  std::vector<double> ticks(5);
  std::vector<std::string> xlabels(5), ylabels(5);
  for (unsigned int i = 0; i < 5; ++i) {
    ticks[i] = i * (ngp - 3) / 4;
    xlabels[i] = std::to_string(eps_lin_seg(ticks[i]));
    ylabels[i] = std::to_string(tau_lin_seg(ticks[i]));
  }

  plt::imshow(kmin, {{"cmap", "viridis"}, {"origin", "lower"}});
  plt::colorbar();
  plt::title("Minimal number of iterations for error < tau");
  plt::xlabel("epsilon_0");
  plt::ylabel("tau");
  plt::xticks(ticks, xlabels);
  plt::yticks(ticks, ylabels);

  // END
  plt::savefig("./cx_out/k_min.png");
}

#endif
