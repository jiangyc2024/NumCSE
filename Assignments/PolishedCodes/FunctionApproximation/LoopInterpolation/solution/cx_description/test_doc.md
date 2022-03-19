

# Tests of adaptiveintp.hpp

***
> `std::vector<Vector2d> closedPolygonalInterpolant(std::vector<Vector2d> Sigma, const VectorXd &x)`: We test this function with `Sigma` being a point evaluation of the function given in (6.10.1) at $50$ points, linearly spaced  over $\[0,1\]$ and `x` contains $33$ values, linearly spaced over $\[0,1\]$.

***

> `std::vector<Vector2d> closedHermiteInterpolant(std::vector<Vector2d> Sigma, const VectorXd &x)`: The inputs for this test are the same as in the previous test.

***

> `std::pair<std::vector<Vector2d>, std::vector<Vector2d>> adaptedHermiteInterpolant(CurveFunctor &&c, unsigned int nmin, const VectorXd &x, double tol = 1.0e-3)`: Here in the test we use for `c` the function fiven in (6.10.1), $5$ for `nmin` and `x` being the same as in the previous testcases.

***
> `void plotKite()`: This function is not tested.
