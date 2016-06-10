function chebexptest(y)
% Testing computation and evaluation of Chebychev expansions
n = length(y) - 1; % degree of polynomial
a = chebexp(y); % FFT based algorithm
ar = chebexpdir(y); % Direct evluation
fprintf('|a-ar| = %f\n',norm(a-ar));
t = cos(((0:n)+0.5)/(n+1)*pi); % Chebychev nodes on \Blue{$[-1,1]$}, see \eqref{eq:CHEBNODES}
yr = clenshaw(a,t);
fprintf('|y-yr| = %f\n',norm(y-yr));
