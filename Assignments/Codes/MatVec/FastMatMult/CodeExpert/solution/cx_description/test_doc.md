

# Tests of gramschmidt.hpp

> There are three functions in strassen.hpp and main.cpp has tests for the first two. The third function is a timing exercise, which should not be tested automatically.

***
> `strassenMatMult(const MatrixXd& A, const MatrixXd& B)`: We test this function by calculating norm(Strassen(A,Id) - A) and also norm(Strassen(Id,A) - A) for a random 128*128 matrix A. The tolerance is 10^{-9}.

***
> `test_strassen()`: This function is a test written by the student. It should create two random matrices A and B, and return err = norm(Strassen(A,B) - A*B).
The test is passed if err < 10^{-9}.