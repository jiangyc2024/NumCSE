#include "imread.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <libgen.h>
#include <span>

// tests the bmp read functions
int main(int argc, char **argv) {
  const std::span args(argv, argc);
  const std::string path = std::string(dirname(args[0])) + "/test.bmp"; //NOLINT(concurrency-mt-unsafe)
  const auto img = imread::readBMP(path);

  std::cout << "Size: " << img.rows() << " x " << img.cols() << std::endl
            << std::endl;

  for (int i = 0; i < 3; ++i) {
    std::cout << "color channel " << i << std::endl;
    const Eigen::MatrixXd G = imread::getcolor(img, i);

    std::cout << G << std::endl;
  }

  std::cout << "greyscale " << std::endl;
  const Eigen::MatrixXd G = imread::greyscale(img);

  std::cout << G << std::endl;
}
