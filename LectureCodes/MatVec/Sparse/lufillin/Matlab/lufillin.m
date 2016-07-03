% Demonstration of fill-in for LU-factorization of sparse matrices
n = 100; 
A = [gallery('tridiag',n,-1,3,-1), speye(n); speye(n) , gallery('tridiag',n,-1,3,-1)]; 
[L,U,P] = lu(A); // LU-decomposition
figure; spy(A); title('Sparse matrix'); 
print -depsc2 '../PICTURES/sparseA.eps';
figure; spy(L); title('Sparse matrix: L factor'); 
print -depsc2 '../PICTURES/sparseL.eps';
figure; spy(U); title('Sparse matrix: U factor'); 
print -depsc2 '../PICTURES/sparseU.eps';




