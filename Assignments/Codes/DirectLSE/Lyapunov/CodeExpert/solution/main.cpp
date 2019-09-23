
#include "lyapunov.hpp"

int main() {
    bool ans = testLyapunov();
    std::cout << "Error less than tol: " << ans << std::endl;
}