function V = sspowitsteprp(A,V)
V = A*V;           % power iteration applied to columns of \Blue{$\VV$}
[Q,R] = qr(V,0);   % orthonormalization, see \cref{sec:evporth}
[U,D] = eig(Q'*A*Q); % Solve Ritz projected \Blue{$m\times m$} eigenvalue problem
V = Q*U;           % recover approximate eigenvectors
ev = diag(D);      % approximate eigenvalues