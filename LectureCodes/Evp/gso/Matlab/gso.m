function Q = gso(V)
% \Hyperlink{GSO}{Gram-Schmidt orthonormalization} of the columns of \Blue{$\VV\in\bbR^{n,m}$}, see
% \eqref{eq:GSO}. The vectors \Blue{$\Vq_1,\ldots,\Vq_m$} are returned as the columns of
% the \Magenta{\emph{orthogonal}} matrix \Blue{$\VQ$}.
m = size(V,2);
Q = V(:,1)/norm(V(:,1)); % normalization
for l=2:m
  q = V(:,l);
  % orthogonalization
  for k=1:l-1
    q = q - dot(Q(:,k),V(:,l))*Q(:,k);
  end
  Q = [Q,q/norm(q)]; % normalization \label{gso:orth}
end
