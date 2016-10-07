function y = lsqsvd(A,b)
[U,S,V] = svd(A,0);
sv = diag(S);
r = max(find(sv/sv(1) > eps));
y = V(:,1:r)*(diag(1./sv(1:r))*...
          (U(:,1:r)'*b));