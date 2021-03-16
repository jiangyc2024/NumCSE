

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
The table printed by the student should look similar to the following table. Note that we do not include a line for $n=2^9$ because we exceed the time limits of Code Expert.
```
              n           Slow           Fast
              4      5.190e-04      6.390e-05
              8      1.769e-03      9.975e-05
             16      4.552e-03      1.081e-04
             32      1.233e-02      1.576e-04
             64      5.138e-02      2.541e-04
            128      2.738e-01      4.427e-04
            256      1.587e+00      8.242e-04
```