#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

int main() {
    // TODO: Compute approximation of the derivative of sin(x)
    // Print the error of each computation

    // Plot
    mgl::Figure fig;
    fig.setlog(true, true);
    fig.legend();
    fig.title("Error of approximation of f'(x_0)");
    fig.xlabel("h");
    fig.ylabel("| f'(x_0) - g_i(x_0, h) |");
    // TODO: plot errors and save figure
}
