# Tests of tridiagqr.hpp

> Note that you can use *compGivensRotation* even if you did not implement it correctly. The tests will automatically call the solution of this function. However, your implementation of this function will also be tested.

> The function *compGivensRotation* from _(4-4.b)_ is tested using three random input vectors.

> The constructor of *TriDiagonalQR* from _(4-4.c)_ is tested by comparing the Q and R factors from a "linspaced" and a random tridiagonal matrix.

> The member function *applyQT* from _(4-4.d)_ is tested using the same tridiagonal matrices as in the constructor test. The vector is a random vector. Please note that this test will only pass if you already implemented the constructor in a correct manner.

> The member function *solve* from _(4-4.e)_ is tested with a "linspaced" tridiagonal matrix and a linspaced vector and with a random tridiagonal matrix and a random vector. It is also tested whether your function throws an exception in case of a singular matrix.

> The specialized function *invit* from _(4-4.g)_ is tested with a random vector.
