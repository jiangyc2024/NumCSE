n = 10; u = (1:n)'/3; v = (1./u).*(-1).^((1:n)'); 
x = ones(10,1); nx = norm(x,'inf');

result = [];
for epsilon = 10.^(-5:-0.5:-14)
    A = u*v' + epsilon*eye(n); 
    b = A*x; nb = norm(b,'inf');
    xt = A\b; % Gaussian elimination
    r = b - A*xt; % residual
    result = [result; epsilon, norm(x-xt,'inf')/nx, norm(r,'inf')/nb, cond(A,'inf')];
end

