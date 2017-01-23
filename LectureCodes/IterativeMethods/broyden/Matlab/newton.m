function [x,prot] = newton(x,F,DF,N,sol)
% MATLAB template for Newton method

prot = [];
for k = 1:N
  s = DF(x)\F(x); x = x-s;
  prot = [prot; k,norm(s),norm(F(x)),norm(x-sol)];
  fprintf('Newton-Iteration %d: |s| = %e, |F(x)| = %e \n',k,norm(s),norm(F(x)));
end
