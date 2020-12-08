
# Tests of quadraticconvergence.hpp

> The test button only tests `steffensen`, for the other functions run the program.

> main.cpp calls the functions in the following way:
***
> `double steffensen(Function &&f, double x0)` is tested by calling `testSteffensen()`: the student should input the function $f(x) = xe^{x} - 1$  
and a suitable initial guess (e.g. $x_0 = 1$). The result should be
```
x = 0.567143
```


***

> `class Logger`: In this class the student should fill little gaps. The overloaded `operator (T val)` should call a `push_back` on the member `std::vector<T> info`.
The test is passed if 
```
Logger<double> test_logger; 
test_logger(2.7);
test_logger(2.0);
test_logger(2.5);
```
results in the vector  $[2.7, 2, 2.5]^T$

***

> `double steffensen_log(Function &&f, double x0, Logger<double> *logger_p = nullptr)`: this is tested by calling the function `orderSteffensen ()`. The test is passed if the sequence of approximate rates
```
1.542345498
1.650553642
1.770024324
1.883754996
1.964598249
1.995899954
1.999927866
0.741551601
```
is printed in a table. The student can format the table freely: the test does not check the presence of other columns.

