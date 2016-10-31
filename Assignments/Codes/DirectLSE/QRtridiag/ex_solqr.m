% Solve a 3-banded LSE with Givens rotations 
% the transformation is given by ex_qrtridiag.m
function ex_solqr(n)

if nargin==0; n=100; end
% initialize matrix and rhs
c = 2*rand(n-1,1)-1;
e = 2*rand(n-1,1)-1;
%d= 2*ones(n  ,1);       % enforce non-singularity
d = 2*rand(n  ,1)-1;
b = 2*rand(n  ,1)-1;     % right hand side
bb= b;

[R,rho] = ex_qrtridiag(c,d,e);

% apply all the rotations to b
for k=1:n-1
    if rho(k) == 1
        gamma = 0;
        sigma = 1; 
    elseif abs(rho(k))<1
        sigma = 2*rho(k); 
        gamma = sqrt(1-sigma^2);
    else
        gamma = 2/rho(k);
        sigma = sqrt(1-gamma^2);
    end 
    G=[gamma sigma; -sigma gamma];
    b(k:k+1) = G * b(k:k+1);
end
% back substitution:
x = R\b;

% compute the residual ||Ax-b|| to check the code:
resid = norm(d.*x+[c.*x(2:n); 0]+[0; e.*x(1:n-1)] - bb) / norm(bb); 
sprintf('Relative Residual = %g \n', resid)