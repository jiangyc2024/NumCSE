#pragma once

#include <iostream>

#include "autofocus.hpp"

int main() {
    std::cout << "\nEnter \"0\" to test all functions.\n"
              << "Enter \"1\" to only test save_image().\n"
              << "Enter \"2\" to only test plot_freq().\n"
              << "Enter \"3\" to only test plotV().\n"
              << "Enter \"4\" to only test autofocus().\n";
  
    int ans = 0;
    std::cin >> ans;
    switch (ans) {
        case 0: std::cout << "Using focus = 0, 1, 2, 3" << std::endl;
                for (unsigned int i = 0; i <= 3; ++i) {
                    std::cout << "Saving image..." << std::endl;
                    save_image(i);
                }

                std::cout << "Using focus = 0, 1, 2, 3" << std::endl;
                for (unsigned int i = 0; i <= 3; ++i) {
                    std::cout << "Saving plot..." << std::endl;
                    plot_freq(i);
                }

                plotV();

                std::cout << "Autofocus returns:"
                << autofocus()
                << std::endl;
                break;
        case 1: std::cout << "Using focus = 0, 1, 2, 3" << std::endl;
                for (unsigned int i = 0; i <= 3; ++i) {
                    std::cout << "Saving image..." << std::endl;
                    save_image(i);
                }
                break;
        case 2: std::cout << "Using focus = 0, 1, 2, 3" << std::endl;
                for (unsigned int i = 0; i <= 3; ++i) {
                    std::cout << "Saving plot..." << std::endl;
                    plot_freq(i);
                }
                break;
        case 3: plotV();
        case 4: std::cout << "Autofocus returns:"
                << autofocus()
                << std::endl;
                break;
    }
}