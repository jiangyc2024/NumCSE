function [x,rn,xk] = pcgbase(evalA,b,tol,maxit,invB,x)
% \texttt{evalA} must pass a handle to a function implementing \texttt{A*x}
% \texttt{invB} is to be a handle to a function providing the action of the
% preconditioner on a vector. The other arguments like for \matlab's \texttt{pcg}.
r = b - evalA(x); rho = 1; rn = []; 
if (nargout > 2), xk = x; end
for i = 1 : maxit
  y = invB(r);
  rho_old = rho; rho = r' * y; rn = [rn,rho];
  if (i == 1), p = y; rho0 = rho; 
  elseif (rho < rho0*tol), return;
  else beta = rho/rho_old; p = y+beta*p; end
  q = evalA(p); alpha = rho /(p' * q);
  x = x + alpha * p; 
  r = r - alpha * q;
  if (nargout > 2), xk = [xk,x]; end
end
