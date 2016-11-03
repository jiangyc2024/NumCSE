#include <fstream>

#include <Eigen/Dense>

#include "pgm.hpp"

using namespace Eigen;

int main() {

    //// EXAMPLE USAGE OF PGM OBJECT

    // This "input" stream object will be used to read an image
    std::ifstream file("image.pgm");
    // This "output" stream object will be used to save the image
    std::ofstream file_out("image_edited.pgm");

    // We create an empty object
    PGMObject p;
    // We "pass" the input stream to the object:
    // the object now contains the content of "file"
    file >> p;

    // PGMObject has a method to obtain is data
    // as a matrix.
    MatrixXd M;
    p.get_data(M);

    // We now set the object to contain a new image
    // we rotate the image and reduce the luminosity
    p.set_data(M.transpose() / 2);

    // We now save the file to disk
    file_out << p;
}
