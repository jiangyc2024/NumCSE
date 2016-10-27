% Transform a 3-banded matrix into an upper triangular matrix R
% using Givens rotations, stored in the parameter vector rho
function [R,rho] = ex_qrtridiag(c,d,e)

n =length(d);
R = spdiags([[e;0], d, [0;c]], [-1 0 1],  n,n);
rho = zeros(n-1,1);

for k=1:n-1 
    % if the term under the diagonal is 0 no rotation is needed
    if R(k+1,k)==0 
        rho(k) = 0;
    else
% the 2x2 matrix is  G= +/-[ R(k,k), R(k+1,k); -R(k+1,k), R(k,k)]/r
% (gamma and sigma are R(k,k)/r and R(k+1,k)/r)
% new matrix:   R_k+1  =  [I_{k-1} 0 0; 0 G 0; 0 0 I] * R_k
        r = norm(R(k:k+1,k));
                
        % 3 different cases:
        if R(k,k)==0
            rho(k) = 1;
        elseif abs(R(k+1,k)) < abs(R(k,k))
            % save the rotation parameter rho
            rho(k) = 0.5 * sign(R(k,k)) * R(k+1,k)/r;
            % transform the matrix with the rotation
            % (store in R the sign of the rotation matrix G)
            R(k:k+1,k:min(k+2,n)) = sign(R(k,k))* ...
            [R(k,k), R(k+1,k); -R(k+1,k), R(k,k)] * R(k:k+1,k:min(k+2,n))/r;
        else
            rho(k) = 2   * sign(R(k+1,k)) / (R(k,k)/r);
            R(k:k+1,k:min(k+2,n)) = sign(R(k+1,k))* ...
            [R(k,k), R(k+1,k); -R(k+1,k), R(k,k)] * R(k:k+1,k:min(k+2,n))/r;
        end
    end
end
