function [x,res] = qrlsqsolve(A,b)
% Solution of linear least squares problem \eqref{eq:LSQ1} by means of QR-decomposition
% Note: \Blue{$\VA\in\bbR^{m,n}$} with \Blue{$m>n$}, \Blue{$\operatorname{rank}(\VA) = n$} is assumed
[m,n] = size(A); 
R = triu(qr([A,b],0)), % economical QR-decomposition of extended matrix \label{qrl:1}
x = R(1:n,1:n)\R(1:n,n+1); % \Blue{$\wh{\Vx} = (\VR)_{1:n,1:n}^{-1}(\VQ^T\Vb)_{1:n}$}
res = R(n+1,n+1); % \Blue{$= \N{\VA\wh{\Vx}-\Vb}_2$} (why ?) \label{qrl:2}
