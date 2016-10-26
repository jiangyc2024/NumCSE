% Solve a 3-banded LSE with Givens rotations 
% 2nd version: enforcing G(1,1)>=0, different sign convention for rho
function ex_solqr_pos(n)
if nargin==0;    n=100;    end
c = 2*rand(n-1,1)-1;   e = 2*rand(n-1,1)-1;
d = 2*rand(n  ,1)-1;   b = 2*rand(n  ,1)-1;   bb= b;

[R,rho] = ex_qrtridiag_pos(c,d,e);
for k=1:n-1
    if rho(k) == 1
        gamma = 0;        sigma = 1; 
    elseif abs(rho(k))<1
        sigma = 2*rho(k);        %get sign from rho
        gamma = sqrt(1-sigma^2); %positive
    else
        gamma = 2/abs(rho(k));   %add abs
        sigma = sign(rho(k))*sqrt(1-gamma^2); %add sign
    end 
    G=[gamma sigma; -sigma gamma];
    b(k:k+1) = G * b(k:k+1);
end
x = R\b;
resid = norm(d.*x+[c.*x(2:n); 0]+[0; e.*x(1:n-1)] - bb) / norm(bb); 
sprintf('Relative Residual = %g \n', resid)

% === transform the matrix and store the rotations in rho-format ===
function [R,rho] = ex_qrtridiag_pos(c,d,e)
n =length(d);
R = spdiags([[e;0], d, [0;c]], [-1 0 1],  n,n);
rho = zeros(n-1,1);
for k=1:n-1 
    % if the term under the diagonal is 0 no rotation is needed
    if R(k+1,k)==0 
        rho(k) = 0;
    else
        G = planerot([R(k,k); R(k+1,k)]);
        G = sign(G(1,1))*G;   % make sure G(1,1) is always positive
        R([k,k+1],k:min(end,k+2)) = G*R([k,k+1],k:min(end,k+2)); 
        gamma = G(1,1);        sigma = G(1,2);
        if gamma==0
           rho(k) = 1;
        elseif abs(sigma) < abs(gamma)
           rho(k) = 0.5 * sigma; % remove sign(gamma), always 1
        else
           rho(k) = 2   * sign(sigma) / gamma;
        end
    end
end