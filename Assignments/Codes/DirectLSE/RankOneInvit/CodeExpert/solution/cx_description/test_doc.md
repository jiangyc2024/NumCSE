

# Tests of rankoneinvit.hpp

> In this exercise, the student implements two versions of an algorithm.
We run both implementations in the following way:
```
std::srand(21);
double tol = 1e-3;
VectorXd d = VectorXd::Random(n);
double lmin = rankoneinvit(d, tol);
double lmin_fast = rankoneinvit_fast(d, tol);
```

1. The first test checks if `lmin` has the correct value: `1.0118`

2. The second test checks if `lmin_fast - lmin` is close to zero (tolerance $10^{-6}$).

In other words, the first test is passed if `rankoneinvit` returns the correct answer. The second test is passed if `rankoneinvit_fast` and `rankoneinvit` return the same value.

***
> There is no test for `rankoneinvit_runtime()`.