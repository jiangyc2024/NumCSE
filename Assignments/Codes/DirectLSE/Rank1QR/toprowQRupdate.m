function [Qh, Rh] = toprowQRupdate(Q,R,v)
[m,n] = size(Q);

% initialize the two matrices
Qh = [1, zeros(1,n); zeros(m,1), Q];
Rh = [v'; R];
% Rh is (n+1)x(n+1) upper Hessenberg: we need n Givens rotations
for k=1:n
    % compute Givens rotations
    [G, Rh(k:k+1,k)] = planerot(Rh(k:k+1,k));
    % update Rh
    Rh(k:k+1,k+1:end) = G * Rh(k:k+1,k+1:end);
    % update Qh multiplying with G' from right
    Qh(:,k:k+1) = Qh(:,k:k+1) * G';
end