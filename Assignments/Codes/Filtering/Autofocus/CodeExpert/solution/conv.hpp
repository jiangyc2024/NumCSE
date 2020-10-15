#pragma once

#include <Eigen/Dense>
#include <iostream>

enum class ConvolutionType { Full, Same, Valid, Periodic };

Eigen::MatrixXd conv2(const Eigen::MatrixXd& P, const Eigen::MatrixXd& S,
                      ConvolutionType type = ConvolutionType::Same) {
  const int m = P.rows(), n = P.cols(), M = S.rows(), N = S.cols();

  switch (type) {
    default:
    case ConvolutionType::Same: {
      Eigen::MatrixXd C(m, n);

      int _tj = N / 2 + 1;
      int _ti = M / 2 + 1;
      for (int j = 0; j < n; ++j) {
        int _lj = std::max(0, _tj + j - n);
        int _rj = std::min(N, _tj + j);
        for (int i = 0; i < m; ++i) {
          double s = 0;
          int _li = std::max(0, i + _ti - m);
          int _ri = std::min(M, i + _ti);
          for (int k = _li; k < _ri; ++k) {
            for (int l = _lj; l < _rj; ++l) {
              int _i = i - k + M / 2;
              int _j = j - l + N / 2;

              s += P(_i, _j) * S(k, l);
            }
          }

          C(i, j) = s;
        }
      }

      //        int _tj = N/2;
      //        int _ti = M/2;
      //        for (int j = 0; j < n; ++j) {
      //            int _lj = std::max(0, _tj - j);
      //            int _rj = std::min(N, _tj - j + n);
      //            for (int i = 0; i < m; ++i) {

      //                int _li = std::max(0, i - _ti);
      //                int _ri = std::min(M, i - _ti + m);

      //                C(i, j) = (
      //                            P.block(i + _li - M/2, i + _ri - M/2,
      //                                    j + _lj - N/2, j + _rj - N/2
      //                                    ).array() *
      //                            S.block(_li, _ri,
      //                                    _lj, _rj
      //                                    ).array()
      //                            ).sum();
      //            }
      //        }

      //        for (int j = 0; j < n; ++j) {
      //            for (int i = 0; i < m; ++i) {

      //                double s = 0;

      //                for (int l = 0; l < N; ++l) {
      //                    for (int k = 0; k < M; ++k) {
      //                        int _i = i - k + M/2;
      //                        int _j = j - l + N/2;
      //                        if(_i < 0 || _i >= m || _j < 0 || _j >= n) {
      //                            continue;
      //                        }

      //                        s += P(_i, _j) * S(k, l);
      //                    }
      //                }

      //                C(i, j) = s;
      //            }
      //        }

      return C;
    }
    case ConvolutionType::Full:
    case ConvolutionType::Valid:
    case ConvolutionType::Periodic:
      std::cerr << "Error: mode not available!" << std::endl;
      exit(-1);
  }
}
