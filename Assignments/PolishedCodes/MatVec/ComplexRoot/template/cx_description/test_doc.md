
# Tests of compelexroot.hpp

> Three numbers are tested for the square root
```
std::complex<double> w(1e20, 5);
w = std::complex<double>(-5, 1e20);
w = std::complex<double>(1e-8, 0);
```

> The expected output is
```
The square root of (1e+20,5) is (1e+10,2.5e-10)
The square root of (-5,1e+20) is (7.07107e+09,7.07107e+09)
The square root of (1e-08,0) is (0.0001,0)
```