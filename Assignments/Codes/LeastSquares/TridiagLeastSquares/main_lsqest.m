function main_lsqest
z = (1:10)';
c = (10:-1:1)';
[alpha, beta] = lsqest(z,c)

% z and c are column vectors
function [alpha, beta] = lsqest(z,c)
A = [z, [z(2:end);0] + [0;z(1:end-1)]];
v = A \ c;              % least squares
alpha = v(1);
beta = v(2);