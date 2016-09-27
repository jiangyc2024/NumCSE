% Solve Lyapunov equation   A*X+X*A'=I
function X = solveLyapunov(A)
n = length(A);
% use sparse identity, it provides sparse B and b
I = speye(n,n);
b = reshape(I,[],1);
B = kron(A,I) + kron(I,A);
x = B\b;
X = reshape(x,n,n);

% check with:
% n=50; A=randn(n); X=solveLyapunov(A); norm(A*X+X*A.'-eye(n))