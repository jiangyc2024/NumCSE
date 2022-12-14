# Tests of mirk.hpp

***
> `Newton2Steps` is tested on $f(x, y) = \begin{bmatrix} e^x \\ e^y \end{bmatrix}$ for initial guess $x = \begin{bmatrix} -1 \\ 1 \end{bmatrix}$ against the correct solution.
***
> `MIRKStep` is tested on the IVP $y' = \lambda y$ with $\lambda = -50$, $h=0.001$, $y_0 = 1.$ against the correct solution.
***
> `MIRKSolve` is tested on the same IVP as before with 100 steps.
***
