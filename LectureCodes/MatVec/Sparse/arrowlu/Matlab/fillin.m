% Simple example for dense inverse despite sparse LU-factors
A = [diag(1:10),ones(10,1);ones(1,10),2];
[L,U] = lu(A); spy(A); spy(L); spy(U); spy(inv(A));
