% This function calculates eigenvalues with the rayleight-quotient iteration,
% which is a improved power method.
%
% Originaly, there were the input parameter d and maxit,
% but in the skript, d = [1:10]' and maxit = 10
% See also rqi.m


function [res,cvg] = rquid()

% Inputparameters
d = [1:10]'; maxit = 10;

% This is for producing matrix A from vector d
n = length(d);
Z = diag(sqrt(1:n),0) + ones(n,n);
[Q,R] = qr(Z);
A = Q*diag(d,0)*Q';

% Subroutine rqi() makes the rayleight-quotient iteration
res = rqi(A,maxit);

% plot and output:
semilogy(res(:,1),res(:,4),'r-*',res(:,1),res(:,3),'m-o');
print -depsc2 rqui.eps


cvg = [res(2:end,1), res(2:end,3)./res(1:end-1,3),  res(2:end,4)./res(1:end-1,4)];
cvg = [cvg,...
	[0;log(cvg(2:end,2))./log(cvg(1:end-1,2))],...
	[0;log(cvg(2:end,3))./log(cvg(1:end-1,3))]];
