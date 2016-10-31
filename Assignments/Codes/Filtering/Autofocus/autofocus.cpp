#include <fstream>

#include <mgl2/mgl.h>
#include <Eigen/Dense>

#include "pgm.hpp"

#include "autofocus.hpp"

#include "FFT/fft2.hpp"


using namespace Eigen;

//! FFT for matrices
//  This just implements FFT for each column of X
MatrixXcd fftr(const MatrixXd& X) {
  const long m = X.rows(), n = X.cols();
  MatrixXcd Y(m, n);
  Eigen::FFT<double> fft;
  for (long j = 0; j < n; ++j) {
    VectorXd Xj = X.col(j);
    Y.col(j) = fft.fwd(Xj);
  }
  return Y;
}

//! Inverse FFT for matrices
//  This just implements inverse FFT for each column of X
MatrixXd ifftr(const MatrixXcd& X) {
  const long m = X.rows(), n = X.cols();
  MatrixXd Y(m, n);
  Eigen::FFT<double> fft;
  for (long j = 0; j < n; ++j) {
    VectorXcd Xj = X.col(j);
    Y.col(j) = fft.inv(Xj);
  }
  return Y;
}

//! 2-dimensional FFT
//  Implementation based on: https://ch.mathworks.com/help/matlab/ref/fft2.html
MatrixXcd fft2r(const MatrixXd& X) {
  return fft(fftr(X).transpose()).transpose();
}

//! 2-dimensional inverse FFT
MatrixXd ifft2r(const MatrixXcd& X) {
  return ifftr(ifft(X).transpose()).transpose();
}

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


double autofocus() {

    unsigned int Niter = 6;

    unsigned int max_focus = 5;
    double df = max_focus / 1e2;
    double f0 = max_focus / 2.;
    double step = max_focus / 2.;

    auto computeV = [] (double focus) {
        return high_frequency_content(
                    fft2r(
                        set_focus(focus)
                        ).cwiseAbs()
                    );
    };

    for(unsigned int i = 0; i < Niter; ++i) {
        double dV = computeV(f0+df) - computeV(f0);

        step = step / 2.;
        f0 = f0 + (dV > 0 ? 1 : -1) * step;
    }

    return f0;
}

int main() {
//    std::ifstream file("image.pgm");
//    PGMObject p;
//    file >> p;

//    MatrixXd M;
//    p.get_data(M);

//    p.set_data(M.transpose() / 2);

//    file_out << p;

    //    PGMObject q;
//    p.set_data(set_focus(4));
//    std::ofstream file_blur("image_blur.pgm");
//    file_blur << p;

    for(unsigned int i = 1; i <= 4; ++i) {
        std::stringstream ss;
        ss << "image_focus"
           << i
           << ".pgm";
        std::ofstream file_out(ss.str());

        MatrixXd D = fft2r(set_focus(i)).cwiseAbs();

//        D

        int a = 0;
        int b = 8000;
        auto clamp = [a,b] (double x) {
            return x < a ? a : x > b ? b : x;
        };
//        .unaryExpr(clamp)
//        D.maxCoeff()

        std::cout << D.maxCoeff();
//        PGMObject p;
//        p.set_data(D.unaryExpr(clamp) / b * 254);

        D = D.unaryExpr(clamp);
//        file_out << p;


        // Plot values of $\mathbf{X}$.
        mglData Xd(D.cols(), D.rows(), D.data());

        mglGraph gr;
//        gr.SubPlot(1,1,0,"<_");
//        gr.SetRanges(0,D.rows(),0,D.cols());
        gr.SetRange('c', 0, b);
        gr.Colorbar("bcwyr");
        gr.Title("Visualization of $X$");
        gr.Axis();
        gr.Tile(Xd, "bcwyr");
        std::stringstream ss2;
        ss2 << "image_focus"
           << i
//           << ".eps";
        << ".png";
//        gr.WriteEPS(ss2.str().c_str());
        gr.WritePNG(ss2.str().c_str());

    }

    unsigned int N = 100;
    for(unsigned int i = 0; i < N; ++i) {
        double V = high_frequency_content(fft2r(set_focus(5. / (N-1) * i)).cwiseAbs());
        std::cout << 5. / (N-1) * i << "\t" << V << std::endl;
    }

    std::cout << "Autofocus returns:"
              << autofocus()
              << std::endl;
}
