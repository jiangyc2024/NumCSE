function [x,prot] = simplenewt(x,F,DF,N,sol)
% MATLAB template for simplified Newton method

J = DF(x);
prot = [];
for k = 1:N
  s = J\F(x); x = x-s;
  prot = [prot; k,norm(s),norm(F(x)),norm(x-sol)];
  fprintf('Simplified Newton-Iteration %d: |s| = %f, |F(x)| = %f \n',...
	  k,norm(s),norm(F(x)));
end
