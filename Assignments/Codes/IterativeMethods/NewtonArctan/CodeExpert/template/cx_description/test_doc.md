

# Tests of newtarctan.hpp


> `double newton_arctan(double x0_ = 2)`: We test this function on the default initial guess $2$. The expected output is
```
x0 = 1.39175
x1 = -1.39175
x2 = 1.39175
```
That is, the sign alternates and the absolute value stays unchanged.