# Tests of quadinf.hpp

***
> `double quadinf(const int n, const Function& f)`: We test this function on the function $f(t) =  \frac{1}{1+t^2}$ and $n = 10$. Its exact integral over $\mathbb{R}$ equals $\pi$. 

> This calls the function 
`double quad(const Function& f, const Eigen::VectorXd& w, const Eigen::VectorXd& x,double a, double b)`, which is not tested since it is optional. 

***

> `void cvgQuadInf(void)`: produces the convergence plot. No test.