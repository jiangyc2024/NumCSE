function [z,lmin] = rqui(A,tol,maxit)
alpha = 0; n = size(A,1);
z = rand(size(A,1),1); z = z/norm(z); % \Blue{$\Vz^{(0)}$}
for i=1:maxit
  z = (A-alpha*speye(n))\z; % \Blue{$\Vz^{(k+1)} = (\VA-\rho_{\VA}(\Vz^{(k)})\VI)^{-1}\Vx^{(k)}$}\label{rqui:1}
  z = z/norm(z); lmin=dot(A*z,z); % Computation of \Blue{$\rho_{\VA}(\Vz^{(k+1)})$}
  if (abs(alpha-lmin) < tol*lmin) % Desired relative accuracy reached ?
    break; end; 
  alpha = lmin;
end
