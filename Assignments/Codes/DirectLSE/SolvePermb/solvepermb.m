function X = solvepermb(A,b)
[n,m] = size(A);
if ((n ~= m)), error('A not squared'); end
if ((n ~= numel(b)) || (m ~= numel(b))), error('Size mismatch'); end
X = [];
for l=1:n
   X = [X,A\b];
   b = [b(end);b(1:end-1)];
end
