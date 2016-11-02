function [ev,V] = sspowitrp(A,k,m,tol,maxit)
% Power iteration with Ritz projection for matrix \Blue{$\VA=\VA^T\in\bbR^{n,n}$}:
% Subspace of dimension \Blue{$m\leq n$} is used to compute the \Blue{$k\leq m$} largest
% eigenvalues of \Blue{$\VA$} and associated eigenvectors.
n = size(A,1); V = eye(n,m); d = zeros(m,1); % (Arbitrary) initial eigenvectors
% The approximate eigenvectors are stored in the columns of \Blue{$\VV\in\bbR^{n,m}$}
for i=1:maxit
  [Q,R] = qr(A*V,0); % Power iteration and orthonormalization
  [U,D] = eig(Q'*A*Q); % Small \Blue{$m\times m$} eigenvalue problem for Ritz projection
  [ev,idx] = sort(diag(D)); % eigenvalue approximations in diagonal of \Blue{$\VD$}
  V = Q*U; % 2nd part of Ritz projection
  if (abs(ev-d) < tol*max(abs(ev))), break; end
  d = ev;
end
ev = ev(m-k+1:end);
V = V(:,idx(m-k+1:end));
