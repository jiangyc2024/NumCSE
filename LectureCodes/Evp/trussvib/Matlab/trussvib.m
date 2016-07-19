function [lambda,V] = trussvib(pos,top)
% Computes vibration modes of a truss structure, see Ex.~\ref{ex:truss}. Mass point
% positions passed in the $n\times 2$-matrix \texttt{poss} and the connectivity encoded in
% the sparse symmetric matrix \texttt{top}. In addition \texttt{top(i,j)} also stores the
% Young's moduli \Blue{$\alpha_{ij}$}.
% The \Blue{$2n$} resonant frequencies are returned in the vector \texttt{lambda}, the
% eigenmodes in the column of \texttt{V}, where entries at odd positions contain the
% \Blue{$x_{1}$}-coordinates, entries at even positions the \Blue{$x_{2}$}-coordinates
n = size(pos,1); % no. of point masses
% Assembly of stiffness matrix according to \eqref{eq:stiffmat}
A = zeros(2*n,2*n);
[Iidx,Jidx] = find(top); idx = [Iidx,Jidx]; % Find connected masses
for ij = idx'
  i = ij(1); j = ij(2);
  dp = [pos(j,1);pos(j,2)] - [pos(i,1);pos(i,2)]; % \Blue{$\Delta\Vp^{ji}$}
  lij = norm(dp);                                 % \Blue{$l_{ij}$} 
  A(2*i-1:2*i,2*j-1:2*j) = -(dp*dp')/(lij^3);
end
% Set Young's moduli \Blue{$\alpha_{ij}$} (stored in \texttt{top} matrix)
A = A.*full(kron(top,[1 1;1 1]));
% Set $2\times2$ diagonal blocks
for i=1:n
  A(2*i-1:2*i,2*i-1) = -sum(A(2*i-1:2*i,1:2:end)')';
  A(2*i-1:2*i,2*i)   = -sum(A(2*i-1:2*i,2:2:end)')';
end
% Compute eigenvalues and eigenmodes
[V,D] = eig(A); lambda = diag(D);

  


