#pragma once

#include "autofocus.hpp"

// TODO: copy this file to ../solution once finished editing
int main() {
    //// SUBPROBLEM a: save differently blurred images
    std::cout << "*** Subproblem a ***"
              << std::endl;
    for(unsigned int i = 0; i <= 3; ++i) {
        std::cout << "Saving image..."
                  << std::endl;
        save_image(i);
    }

    //// SUBPROBLEM b: plot spectrum for different $f$
    std::cout << "*** Subproblem b ***"
              << std::endl;
    for(unsigned int i = 0; i <= 3; ++i) {
        std::cout << "Saving plot..."
                  << std::endl;
        plot_freq(i);
    }

    //// SUBPROBLEM c: plot $V(\mathbf{B}(f))$
    std::cout << "*** Subproblem c ***"
              << std::endl;
    plotV();

    //// SUBPROBLEM d: find most focused image
    std::cout << "*** Subproblem d ***"
              << std::endl;
    std::cout << "Autofocus returns:"
              << autofocus()
              << std::endl;
}
