function [lmin,y] = invit(A,tol)
[L,U] = lu(A); % \Magenta{single} intial LU-factorization, see Rem.~\ref{rem:seqsolvelse}
n = size(A,1); x = rand(n,1); x = x/norm(x); % random initial guess
y = U\(L\x); lmin = 1/norm(y); y = y*lmin; lold = 0;
while(abs(lmin-lold) > tol*lmin) % termination, if small \emph{relative change}
  lold = lmin; x = y; 
  y = U\(L\x); % core iteration: \Blue{$\Vy = \VA^{-1}\Vx$}, 
  lmin = 1/norm(y); % new  approxmation of \Blue{$\lambda_{\min}(\VA)$}
  y = y*lmin;  % normalization \Blue{$\Vy := \frac{\Vy}{\N{\Vy}_2}$}
end

