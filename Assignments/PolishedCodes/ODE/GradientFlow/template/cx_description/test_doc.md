# Tests of gradientflow.hpp

***
> We test using $f(\bm{x}) = \begin{bmatrix} ||\bm{x}||^2 \\ (\bm{x})_0  (\bm{x})_1\end{bmatrix}$, $h = 0.1$, $n = 30$ and $y_0 = \begin{bmatrix}10\\1 \end{bmatrix}$.
***
> The stages from `computeStages` are compared against the solution.
***
> For `discEvolSDIRK`, we check the final state against the solution.
***
> `solveGradientFlow` is tested with random $\bm{d}$ and $\bm{y}$, $\lambda=0.5$, $T=10$, $M=50$. Again, the last state is compared against the solution.
***
