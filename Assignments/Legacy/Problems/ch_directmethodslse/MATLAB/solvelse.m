% Input:
% nxn     upper triangular matrix R,
% nx1     vector v,
% nx1     vector u,
% (n+1)x1 vector b:
% Solve x=[R, v; u',0]\b

function x = solvelse(R,v,u,b)
n = length(b)-1;
% solve L y = b through forward insertion:
y = zeros(n+1,1);
% The first n entries are trivial because of the identity
% in the first n lines of L:
y(1:n) = b(1:n);
% y(n+1) needs the scalar u'*inv(R)*y(1:n):
y(n+1) = b(n+1) - u'*solveR(R,y(1:n));

% Now solve 'U x = y' through backward insertion:
x = zeros(n+1,1);
% the last entry needs the scalar -u'*inv(R)*v:
x(n+1) = y(n+1) / ( -u'*solveR(R,v) );
% The other entries are x(1:n) = inv(R) * (y(1:n) - v*x(n+1)):
x(1:n) = solveR(R, y(1:n) - v * x(n+1));
end

% ============================
% This function determines the solution of 'Rx = b'
% through backward insertion.
% R must be a nxn upper triangular matrix and b an nx1 vector.
function x = solveR(R,b)
n = length(b);
x = zeros(n,1);
% use the triangular structure of R:
x(n) = b(n)/R(n,n);
for i = 1 : (n-1)
    x(n-i) = b(n-i);
    for j = 0 : (i-1)
        x(n-i) = x(n-i)-R(n-i,n-j)*x(n-j);
    end
    x(n-i) = x(n-i)/R(n-i,n-i);
end
end

%  Check your program with:
% n=10; R=triu(rand(n)); v=rand(n,1); u=rand(n,1); b=rand(n+1,1);
% norm(solvelse(R,v,u,b)-[R,v;u',0]\b)