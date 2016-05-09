function y = kronmultv(A,B,x)
n = size(A,2); k = size(B,2);
if (numel(x) ~= n*k), error('size mismatch'); end
% Chop up vector into n pieces of length k
X = reshape(x,k,n); 
% Multiply each piece with matrix B
T = B*X; % \Label[line]{kv:1}
% Form weighted sums of products
% weights supplied by rows of A
T = T*transpose(A); % \Label[line]{kv:2}
% Put matrix columns on top of each other
y = T(:);