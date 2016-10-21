function [a,Q] = givenscoltrf(a)
% Orthogonal transformation of a (column) vector into a multiple of the 
% first unit vector by successive Givens transformations 
n=length(a);
Q=eye(n); % Assemble rotations in matrix, alternative see Rem.~\ref{rem:orthstore}
for j=2:n
  G=planerot([a(1);a(j)]); % see Code~\ref{planerot}
  a([1,j])=G*a([1,j]);
  Q(:,[1,j])=Q(:,[1,j])*G';
end

