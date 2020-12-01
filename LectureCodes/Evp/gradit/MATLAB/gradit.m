function x = gradit(A,b,x,rtol,atol,maxit)
r = b-A*x; % residual $\to$ Def.~\ref{def:residual}
for k=1:maxit
  p = A*r;
  ts = (r'*r)/(r'*p); % \emph{cf.} \eqref{eq:tmin}
  x = x + ts*r;
  cn = (abs(ts)*norm(r); % norm of correction
  if ((cn < rtol*norm(x)) || (cn < atol))
    return;
  end
  r = r - ts*p; % \label{gi:10}
end
