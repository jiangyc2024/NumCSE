% Study of fill-in with LU-factorization due to pivoting
n = 10; D = diag(1./(1:n));
A = [ D , 2*ones(n,1); 2*ones(1,n), 2];
[L,U,P] = lu(A);
figure; spy(A,'r'); title('{\bf arrow matrix A}');
print -depsc2 '../PICTURES/fillinpivotA.eps';
figure; spy(L,'r'); title('{\bf L factor}');
print -depsc2 '../PICTURES/fillinpivotL.eps';
figure; spy(U,'r'); title('{\bf U factor}');
print -depsc2 '../PICTURES/fillinpivotU.eps';
        