function [L,U] = luenv(A)
% envelope aware recursive LU-factorization
% of structurally symmetric matrix
n = size(A,1); 
if (size(A,2) ~= n), 
  error('A must be square'); end
if (n == 1), L = eye(1); U = A; 
else
  mr = rowbandwidth(A); 
  [L1,U1] = luenv(A(1:n-1,1:n-1));
  u = substenv(L1,A(1:n-1,n),mr); 
  l = substenv(U1',A(n,1:n-1)',mr); 
  if (mr(n) > 0)
    gamma = A(n,n) - l(n-mr(n):n-1)'*u(n-mr(n):n-1);
  else gamma = A(n,n); end
  L = [L1,zeros(n-1,1); l' , 1];
  U = [U1,u;zeros(1,n-1) , gamma];
end
