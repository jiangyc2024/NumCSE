function y = clenshaw(a,x)
% Clenshaw algorithm for evaluating \Blue{$p = \sum_{j=1}^{n+1}a_j T_{j-1}$} at points passed in vector \texttt{x}
n = length(a)-1; % degree of polynomial
d = repmat(reshape(a,n+1,1),1,length(x));
for j=n:-1:2
  d(j,:) = d(j,:) + 2*x.*d(j+1,:); % see \eqref{eq:cstr}
  d(j-1,:) = d(j-1,:) - d(j+1,:);
end
y = d(1,:) + x.*d(2,:);
