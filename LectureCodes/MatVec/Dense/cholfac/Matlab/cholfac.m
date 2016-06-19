function R = cholfac(A)
% simple Cholesky factorization
n = size(A,1);
for k = 1:n
  for j=k+1:n
    A(j,j:n) = A(j,j:n) - A(k,j:n)*A(k,j)/A(k,k);
  end
  A(k,k:n) = A(k,k:n)/sqrt(A(k,k));
end
R = triu(A);
