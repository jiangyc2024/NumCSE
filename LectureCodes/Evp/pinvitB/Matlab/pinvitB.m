function [lmin,z,res] = pinvit(evalA,n,invB,tol,maxit)
% \texttt{invB} $\hat{=}$ handle to function implementing preconditioner \Blue{$\VB^{-1}$}
z = (1:n)'; z = z/norm(z); % initial guess
res = []; rho = 0;
for i=1:maxit
  v = evalA(z); rhon = dot(v,z); % Rayleigh quotient
  r = v - rhon*z;    % residual 
  z = z - invB(r);   % iteration according to \eqref{eq:pcinvit}
  z = z/norm(z);     % normalization
  res = [res; rhon]; % tracking iteration
  if (abs(rho-rhon) < tol*abs(rhon)), break;
  else rho = rhon; end
end
lmin =  dot(evalA(z),z); res = [res; lmin],

