function [x,y] = levinson(u,b)
k = length(u)-1;
if (k == 0), x=b(1); y = u(1); return; end
[xk,yk] = levinson(u(1:k),b(1:k));
sigma = 1-dot(u(1:k),yk);
t = (b(k+1)-dot(u(k:-1:1),xk))/sigma;
x = [ xk-t*yk(k:-1:1);t];
s = (u(k+1)-dot(u(k:-1:1),yk))/sigma;
y = [yk-s*yk(k:-1:1); s];

