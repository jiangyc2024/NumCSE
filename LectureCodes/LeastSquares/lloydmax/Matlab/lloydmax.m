function [C,idx,cds] = lloydmax(X,C,tol)
% Lloyd-Max iterative vector quantization algorithm for discrete point sets; the
% columns of \texttt{X} contain the points \Blue{$\Vx_{i}$}, the columns of 
% \texttt{C} initial approximations for the centers of the clusters. The final centers
% are returned in \texttt{C}, the index vector \texttt{idx} specifies the association
% of points with centers.
k = size(X,1); % dimension of space
N = size(X,2); % no. of points
n = size(C,2); % no. of clusters
if (k ~= size(C,1)), error('dimension mismatch'); end
if (nargin < 3), tol = 0.0001; end

sd_old = realmax; 
[sd,idx] = distcomp(X,C),
% Terminate, if sum of squared minimal distances has not changed much
while ((sd_old-sd)/sd > tol)
  % Compute new centers of gravity according to \eqref{eq:cgrav}
  for j=1:n
    idj = find(idx == j); 
    nj = length(idj);
    if (nj > 0), C(:,j) = sum(X(:,idj)')'/nj; end
  end
  sd_old = sd;
  [sd,idx,cds] = distcomp(X,C), 
end

end

function [sumd,idx,cds] = distcomp(X,C)
% Compute squared distances 
d = [];
for j=1:size(C,2)
  Dv = X - repmat(C(:,j),1,size(X,2));
  d = [d; sum(Dv.*Dv)];
end
% Compute minimum distance point association and sum of minimal squared distances
[mx,idx] = min(d);
sumd = sum(mx);
% Computer sum of squared distances within each cluster
for j=1:size(C,2)
  cds(j) = sum(mx(find(idx == j)));
end
end