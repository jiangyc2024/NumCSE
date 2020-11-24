# Tests of chebyApprox.hpp

***
> The function `lengthEllipse` form _(7-5.a)_ is tested with random $\rho>1$ and fixed $n=100$ against
> a series expansion approximation for the arc length of an ellipse

***
> The function `bestBound` from _(7-5.d)_ is tested against the reference solution.
> and the output should be 
```
10.0753467546, 1.9129663501
```
***
> The function `compareErrorAndBound` from _(7-5.f)_ is test against the reference solution with `n_max = 20`,
> and the output should be
```
Order Error          Upperbound     
4     0.0122079      0.409276       
5     0.00263694     0.209879       
6     0.00202698     0.108634       
7     0.000423312    0.0564976      
8     0.000330751    0.0294552      
9     6.79868e-05    0.0153762      
10    5.3577e-05     0.00803203     
11    1.09199e-05    0.00419713     
12    8.64874e-06    0.00219361     
13    1.75395e-06    0.00114659     
14    1.3936e-06     0.000599344    
15    2.8172e-07     0.000313297    
16    2.24322e-07    0.000163773    
17    4.52499e-08    8.56114e-05    
18    3.60853e-08    4.47531e-05    
19    7.26804e-09    2.33945e-05    
20    5.80245e-09    1.22294e-05    
```