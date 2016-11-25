function [V,alph,bet] = lanczos(A,k,z0)
% Note: this implementation of the Lanczos process also records the orthonormal CG residuals in the columns of the matrix \texttt{V}, which is not needed when only eigenvalue approximations are desired.
V = z0/norm(z0); 
% Vectors storing entries of tridiagonal matrix \eqref{eq:tdcg}
alph=zeros(k,1); bet = zeros(k,1);
for j=1:k
  q = A*V(:,j); alph(j) = dot(q,V(:,j));
  w = q - alph(j)*V(:,j);
  if (j > 1), w = w - bet(j-1)*V(:,j-1); end;
  bet(j) = norm(w); V = [V,w/bet(j)];
end
bet = bet(1:end-1);
