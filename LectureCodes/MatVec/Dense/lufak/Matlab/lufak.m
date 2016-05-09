function [L,U] = lufak(A)
% Algorithm of Crout:
% LU-factorization of \Blue{$\VA\in\bbK^{n,n}$}
n = size(A,1); if (size(A,2) ~= n), error('n ~= m'); end
L = eye(n); U = zeros(n,n); 
for k=1:n
  % Compute row of \Blue{$\VU$}
    for j=k:n
        U(k,j) = A(k,j) - L(k,1:k-1)*U(1:k-1,j);
    end
  % Compute column of \Blue{$\VL$}  
    for i=k+1:n
        L(i,k) = (A(i,k) - L(i,1:k-1)*U(1:k-1,k)) /U(k,k);
    end
end
        