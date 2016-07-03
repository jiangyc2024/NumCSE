n = 10; u = (1:n)'/3; v = (1./u).*(-1).^u; 
x = ones(10,1); nx = norm(x,'inf');

result = [];
for epsilon = 10.^(-5:-0.5:-14)
    A = u*v' + epsilon*rand(n,n); 
    b = A*x; nb = norm(b,'inf');
    xt = A\b;      % Gaussian elimination
    r = b - A*xt;  % residualB
    B = inv(A); xi = B*b; 
    ri = b - A*xi; % residual
    R = eye(n) - A*B; % residual
    result = [result; epsilon, norm(r,'inf')/nb, norm(ri,'inf')/nb, norm(R,'inf')/norm(B,'inf') ];
end

