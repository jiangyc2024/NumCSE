

# Tests of strassen.hpp

> Please note that the tests are independent. This means you can call functions from previous subproblems even if you did not solve them or failed to solve them. 

> There are three functions in strassen.hpp and main.cpp. The first two are tested in tests. The third function is a timing exercise, which should not be tested automatically.

***
> `strassenMatMult(const MatrixXd& A, const MatrixXd& B)`: We test this function by calculating norm(Strassen(A,Id) - A) and also norm(Strassen(Id,A) - A) for a random 128 $\times$ 128 matrix A. The tolerance is 10$^{-9}$.

***
> `test_strassen()`: This function has to be completed. It should create two random matrices A and B, and return err = norm(Strassen(A,B) - A$\cdot$B).
The test is passed if err < 10$^{-9}$.
