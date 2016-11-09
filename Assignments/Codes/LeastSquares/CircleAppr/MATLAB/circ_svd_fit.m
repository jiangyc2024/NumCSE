function [m,r] = circ_svd_fit(x,y)
A = [x(:).^2+y(:).^2, x(:), y(:), ones(size(x(:)))];
[U,S,V] = svd(A);
v = V(:,end);
m = -v(2:3)/(2*v(1));
r = sqrt((v(2)^2 + v(3)^2)/(4*v(1)^2) - v(4)/v(1));