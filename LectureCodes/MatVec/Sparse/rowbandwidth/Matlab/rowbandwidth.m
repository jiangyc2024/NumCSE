function mr = rowbandwidth(A)
% computes row bandwidth numbers \Blue{$m_i^R(\VA)$} of \Blue{$\VA$}
n = size(A,1); mr = zeros(n,1);
for i=1:n, mr(i) = max(0,i-min(find(A(i,:) ~= 0))); end
