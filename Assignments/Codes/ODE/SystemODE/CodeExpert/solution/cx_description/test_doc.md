

# Tests of system.hpp

Test `rk4step()` using $f(y) = (y_2y_3, y_1y_2, 3y_3)$, $h=\frac{1}{2}$ and $y^{(0)}=(-1,1,2)$.
The correct output is $y^{(1)} = (0.931516 0.879332  8.79688)$.

Test of `testcvgRK4()` simply checks if `std::round(testcvgRK4())` is equal to $4$.