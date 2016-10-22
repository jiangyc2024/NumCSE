function [v,w] = sspowitstep(A,v,w)
v = A*v; w = A*w; 
% v = v/norm(v); w = w - dot(v,w)*v; w = w/norm(w);
[Q,R] = qr([v,w],0); 
[V,D] = eig(Q'*A*Q);
w = Q*V(:,1); v = Q*V(:,2); 