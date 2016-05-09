function Q = gramschmidt(A)
% Gram-Schmidt orthogonalization of column vectors
% Arguments: Matrix A passes vectors in its columns
% Return values: orthornormal system in columns of matrix Q
[n,k] = size(A); % Get number k of vectors and dimension n of space
Q = A(:,1)/norm(A(:,1)); % First basis vector
for j=2:k
  q = A(:,j) - Q*(Q'*A(:,j)); % Orthogonal projection; loop-free implementation \Label{gs:op}
  nq = norm(q);               % Check premature termination 
  if (nq < (1E-9)*norm(A(:,j))), break; end % Safe check for == 0 \Label{gs:1}
  Q = [Q,q/nq]; % Add new basis vector as another column of Q
end
