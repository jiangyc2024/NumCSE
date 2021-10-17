
# Tests of polydiv.hpp

> We use $u = 1 + 2x + 3x^2 + 4x^3$ and $v = 10 + 20x + 30x^3$ as our two testing polynomials.

> We test the correctness of the coef vector of the product $uv = 10 + 40x + 100x^2 + 160x^3 + 170x^4 + 120x^5$ 
> computed from function `polyMult_naive` and `polyMult_fast`

> We also test the correctness of the coef vector of the division $v = uv / u$ 
> computed from function `polyDiv`
