function x = gausselimsolve(A,b)
% Gauss elimination without pivoting, \texttt{x = A\symbol{92}b}
% \texttt{A} must be an \Blue{$n\times n$}-matrix, \texttt{b} an \Blue{$n$}-vector
n = size(A,1); A = [A,b]; % \label{gse:1}
% Forward elimination (\textit{cf.} step \ding{192} in Ex.~\ref{ex:GE})
for i=1:n-1, pivot = A(i,i);
  for k=i+1:n, fac = A(k,i)/pivot;
      A(k,i+1:n+1) = A(k,i+1:n+1) - fac*A(i,i+1:n+1); % \label{gse:vec}
  end
end
% \Hyperlink{RUECKSUBST}{Back substitution} (\textit{cf.} step \ding{193} in Ex.~\ref{ex:GE})
A(n,n+1) = A(n,n+1) /A(n,n);
for i=n-1:-1:1
    for l=i+1:n
        A(i,n+1) = A(i,n+1) - A(l,n+1)*A(i,l);
    end
    A(i,n+1) = A(i,n+1)/A(i,i);
end
x = A(:,n+1); % \label{gse:last}



