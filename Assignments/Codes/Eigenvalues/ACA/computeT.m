function T = computeT (B,L)
T = B(:,1:L)*B(:,1:L)';
