function [m, r] = circ_alg_fit(x, y)
% algebraic L.S. fitting for circles
% data coordinates (x,y) --> center (m) and radius (r)

A = [-2*x(:), -2*y(:), -ones(size(x(:)))];
b = -(x(:).^2 + y(:).^2);
z = A \ b;
m = z(1:2);
r = sqrt(z(3) + m(1)^2 + m(2)^2);