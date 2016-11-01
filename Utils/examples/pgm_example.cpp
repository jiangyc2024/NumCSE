#include <fstream>

#include <Eigen/Dense>

#include "pgm.hpp"

using namespace Eigen;

int main() {

    //// EXAMPLE USAGE OF PGM OBJECT

    std::ifstream file("image.pgm");
    std::ofstream file_out("image_edited.pgm");

    PGMObject p;
    file >> p;

    MatrixXd M;
    p.get_data(M);

    p.set_data(M.transpose() / 2);

    file_out << p;
}
