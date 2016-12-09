#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>
#include <iomanip>

#include "errors.hpp"

using namespace Eigen;

/* SAM_LISTING_BEGIN_1 */
template <class Function, class State>
void rk4step(const Function &odefun, double h,
             const State & y0, State & y1)
{
    // TODO: implement a single step of the classical Runge-Kutta method of order 4
}
/* SAM_LISTING_END_1 */

int main() {
/* SAM_LISTING_BEGIN_0 */
    // TODO: compute $f$, $T$, $y_0$, $\VA$, $\Vb$ to run "errors(f, T, y0, A, b);"
/* SAM_LISTING_END_0 */
}
