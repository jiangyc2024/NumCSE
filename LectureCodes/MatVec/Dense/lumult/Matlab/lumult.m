function A = lumult(L,U)
% Multiplication of normalized 
% lower/upper triangular matrices
n = size(L,1); A = zeros(n,n);
if ((size(L,2) ~= n) || (size(U,1) ~= n) || (size(U,2) ~= n))
    error('size mismatch'); end
for k=n:-1:1
    for j=k:n
        A(k,j) = U(k,j) + L(k,1:k-1)*U(1:k-1,j);
    end
    for i=k+1:n
        A(i,k) = L(i,1:k-1)*U(1:k-1,k) + L(i,k)*U(k,k);
    end
end
    