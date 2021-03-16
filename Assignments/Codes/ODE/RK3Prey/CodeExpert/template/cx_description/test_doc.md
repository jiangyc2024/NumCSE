

# Tests of `rkintegrator.hpp`

Test `RKIntegrator` the Butcher tableau 
$
\begin{matrix}
0 & 0\\
1 & 0\\
\hline
\frac{1}{2} & \frac{1}{2}
\end{matrix}
$
and solving
$y'=f(y) = (-y_1/2, y_1y_2)$, $y^{(0)}=(-1,1)$ up to time $T=2$ in $N=10$ steps.
The correct output is $y^{(N)} = (-0.368541  0.284085)$.

Test of `RK3prey()` simply checks if `std::round(std::abs(RK3prey()))` is equal to $3$.
Taking absolute value in order to avoid unnecessary confusion.
