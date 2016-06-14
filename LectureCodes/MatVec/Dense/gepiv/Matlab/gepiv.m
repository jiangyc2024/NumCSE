function x = gepiv(A,b)
% Solving an LSE \Blue{$\VA\Vx=\Vb$} by Gaussian elimination with partial pivoting
n = size(A,1); A = [A,b]; % \label{lup:0}
% Forward elimination by rank-1 modification, see Rem.~\ref{rem:blockgs}
for k=1:n-1
   [p,j] = max(abs(A(k:n,k))./max(abs(A(k:n,k:n))')') %\label{lup:1}
   if (p < eps*norm(A(k:n,k:n),1)), %\label{lup:2} 
     disp('A nearly singular'); end 
   A([k,j+k-1],k:n+1) = A([j+k-1,k],k:n+1);                     % \label{lup:3} 
   A(k+1:n,k+1:n+1) = A(k+1:n,k+1:n+1)-(A(k+1:n,k)*A(k,k+1:n+1))/A(k,k); % \label{lup:4} 
end
% \Hyperlink{RUECKSUBST}{Back substitution} (same as in Code~\ref{mc:gausselimsolve})
A(n,n+1) = A(n,n+1) /A(n,n);
for i=n-1:-1:1
  A(i,n+1) = (A(i,n+1) - A(i,i+1:n)*A(i+1:n,n+1))/A(i,i);
end
x = A(:,n+1); % \label{lup:last}

