
% Problem 5c)
% INPUT
% X - points 
% OUTPUT
% B - matrix of the overdetermined system

function B = shapeidentmat(X)

% Preset B to its size
B = zeros(30, 4);

% Assign B's entries aptly
a1 = [1 0 0 0; 0 0 1 0];
a2 = [0 1 0 0; 0 0 0 1];
for i = 1:15
    B(2*i-1:2*i,:) = X(1,i)*a1+X(2,i)*a2;
end

