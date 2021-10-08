
# Tests of adaptedlinreg.hpp

> To test `linReg()` we define 
```
  n = 25;
  t = VectorXd::LinSpaced(n,0.0,1.0);
  std::srand(41);
  noise = VectorXd::Random(n);
```
and use `linReg()` to fit a line to `(t,y)` where
```
y = 12*t + VectorXd::Constant(n,-154) + noise;
```
The correct output is:
```
  (alpha, beta) = -153.796  11.6682
```

***
> To test `expFit()` we define `t` and `noise` as before.
Now we fit an exponential curve to `(t,y)` where
```
y = 17*exp( -3*t.array() ).matrix() + 0.5*noise;
```
The correct output is
```
  (alpha, beta) =  17.8476 -3.14235
```