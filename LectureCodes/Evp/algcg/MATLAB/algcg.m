function x = cg(evalA,b,x,tol,maxit)
% \texttt{x} supplies initial guess, \texttt{maxit} maximal number of CG steps
% \texttt{evalA} must pass a \Red{handle} to a MATLAB function realizing \texttt{A*x}
r = b - evalA(x); rho = 1; n0 = norm(r);
for i = 1 : maxit
   rho1 = rho; rho = r' * r;
   if (i == 1), p = r; 
   else  beta = rho/rho1; p = r + beta * p; end
   q = evalA(p); alpha = rho/(p' * q);
   x = x + alpha * p;    % update of approximate solution
   if (norm(b-A*x) <= tol*n0) return; end % termination, see Rem.~\ref{rem:cgterm}
   r = r - alpha * q;      % update of residual
end  