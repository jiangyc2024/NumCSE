function [a,b] = trigipequid(y)
% Efficient computation of coefficients in expansion \eqref{eq:trigpreal} for a trigonometric
% interpolation polynomial in equidistant points \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$}
% \texttt{y} has to be a row vector of odd length, return values are column vectors
N = length(y); if (mod(N,2)~=1), error('#pts odd!'); end; 
n = (N-1)/2; 
c = fft(exp(2*pi*i*(n/N)*(0:2*n)).*y)/N; % see \eqref{tip:FM}
% From \eqref{eq:trigpcomp}: \Blue{$\alpha_j = \frac{1}{2}(\gamma_{n-j}+\gamma_{n+j})$} and \Blue{$\beta_j = \frac{1}{2i}(\gamma_{n-j}-\gamma_{n+j})$}, \Blue{$j=1,\ldots,n$}, \Blue{$\alpha_0 = \gamma_n$}
a = transpose([c(n+1),c(n:-1:1)+c(n+2:N)]);
b = transpose(-i*[c(n:-1:1)-c(n+2:N)]);
