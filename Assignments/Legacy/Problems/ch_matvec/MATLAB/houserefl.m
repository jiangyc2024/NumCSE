function Z = houserefl(v)
% Porting of houserefl.cpp to Matlab code
% v is a column vector
    % Size of v
    n = size(v,1);
    
    w = v/norm(v);
    u = w + [1;zeros(n-1,1)];
    q = u/norm(u);
    X = eye(n) - 2*q*q';
    
    % Remove first column X(:,1) \in span(v)
    Z = X(:,2:end);
end
