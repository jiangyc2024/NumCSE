# Tests of laserquad.hpp

> Please note that tests in this problem are independent, i.e. you can have an incorrect implementation of a function from a subproblem, call it in another function that depends on it then the tests will pass if your implementation of this particular subproblem is correct.

> Test `evalquad(-2,3,f,Q)` for $f(x) = x^3 \mathrm{log}(|x|+1)$ and 
```
    Q.nodes_ << -1., -std::sqrt(3./7.), 0, std::sqrt(3./7.), 1.;
    Q.weights_ << 0.1, 49./90., 32./45., 49./90., 0.1;
```
The output should be equal to $20.9099$.

***

> Test `gaussquadtriangle(g,5)` for $g(x,y) = x(y-1)\mathrm{sin}(6xy)$.
The output should be $-0.0585825$.


***

>Run the program to test  `convtest2DQuad()`
