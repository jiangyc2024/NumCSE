function [i1,i2] = princaxissep(X)
% Separation of a set of points whose coordinates are stored in the columns of
% X according to their location w.r.t. the principal axis

N = size(X,2);         % no. of. points
g = sum(X')'/N;        % Compute center of gravity, \emph{cf.} \eqref{eq:cgrav}
Y = X - repmat(g,1,N); % Normalize point coordinates.

% Compute \Red{principal axes}, \emph{cf.} \eqref{eq:pax} and \eqref{lsq:maxconst}. Note
% that the SVD of a symmetric matix is available through an orthonormal basis of
% eigenvectors.
[V,D] = eig(Y*Y'); 
a = V(:,end);      % Major principal axis
c = a'*Y;    % Coordinates of points w.r.t. to major principal axis
% Split point set according to locations of projections on principal axis
i1 = find(c < 0); i2 = find(c >= 0);
