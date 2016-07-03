function d = eigqr(A,tol)
n = size(A,1);
while (norm(tril(A,-1)) > tol*norm(A))
% shift by eigenvalue of lower right 2$\times$2 block closest to $(\VA)_{n,n}$    
  sc = eig(A(n-1:n,n-1:n));      
  [dummy,si] = min(abs(sc-A(n,n)));
  shift = sc(si);
  [Q,R] = qr( A - shift * eye(n));
  A = Q'*A*Q;
end
d = diag(A);
