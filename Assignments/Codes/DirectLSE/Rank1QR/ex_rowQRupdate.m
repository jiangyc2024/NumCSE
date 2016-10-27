function ex_rowQRupdate(m,n,j)
% parameters:
if nargin==0;  m=200; n=100; j=25;  end;
if n>m;   n=m;   end;
if j>m;   j=m;   end;

A = rand(m,n);
v = [1:n]';

% A tilde:
At = [A(1:j,:); v'; A(j+1:end,:)];
% Permutation matrix
P = [zeros(1,j), 1 ,zeros(1,m-j);
     eye(j), zeros(j, m+1-j);
     zeros(m-j,j+1), eye(m-j) ];
% A hat:
Ah = P*At;

% economical decomposition of A (see help qr):
[Q,R] = qr(A, 0);
% decomposition of Ah:
[Qh, Rh] = toprowQRupdate(Q,R,v);
% decomposition of At:
Qt = P'*Qh;
Rt = Rh;

% test the code:
RelativeError = norm(At-Qt*Rt) / norm(At)