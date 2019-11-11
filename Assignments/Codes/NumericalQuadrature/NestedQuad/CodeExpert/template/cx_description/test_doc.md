
# Tests of laserquad.hpp

> Test `evalquad(-2,3,f,Q)` for $f(x) = x^3 \mathrm{log}(|x|+1)$ and 
```
    Q.nodes << -1., -std::sqrt(3./7.), 0, std::sqrt(3./7.), 1.;
    Q.weights << 0.1, 49./90., 32./45., 49./90., 0.1;
```
The output should be equal to $20.9099$.

***

> Test `gaussquadtriangle(g,5)` for $g(x,y) = x(y-1)\mathrm{sin}(6xy)$.
The output should be $-0.0585825$.


***

> To test `convtest2DQuad()`, we only check if the correct error is printed for N=2. I.e. if the number $0.00129915$ is on the screen.