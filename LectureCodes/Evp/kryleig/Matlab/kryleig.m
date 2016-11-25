function [V,D] = kryleig(A,m)
% Ritz projection onto Krylov subspace. An orthonormal basis of \Blue{$\Kryl{m}{\VA}{\mathbf{1}}$} is assembled into the columns of \Blue{$\VV$}. 
n = size(A,1); V = (1:n)'; V = V/norm(V); 
for l=1:m-1
  V = [V,A*V(:,end)]; [Q,R] = qr(V,0); 
  [W,D] = eig(Q'*A*Q); V = Q*W;
end
