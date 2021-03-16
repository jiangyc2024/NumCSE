#include "imread.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <libgen.h>

// tests the bmp read functions
int main(int arc, char **argv) {
  std::string path = std::string(dirname(argv[0])) + "/test.bmp";
  auto img = readBMP(path);

  std::cout << "Size: " << img.rows() << " x " << img.cols() << std::endl
            << std::endl;

  for (int i = 0; i < 3; ++i) {
    std::cout << "color channel " << i << std::endl;
    Eigen::MatrixXd G = getcolor(img, i);

    std::cout << G << std::endl;
  }

  std::cout << "greyscale " << std::endl;
  Eigen::MatrixXd G = greyscale(img);

  std::cout << G << std::endl;
}
