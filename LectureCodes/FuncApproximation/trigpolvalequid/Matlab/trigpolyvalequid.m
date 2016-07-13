function q = trigpolyvalequid(y,N)
% Evaluation of trigonometric interpolation polynomial through \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$}
% in equidistant points \Blue{$\frac{k}{M}$}, \Blue{$k=0,M-1$}
% \texttt{y} has to be a row vector of odd length, values returned as rwo vector, too.
N = length(y); if (mod(N,2)~=1), error('#pts odd!'); end; 
n = (N-1)/2; 
% Compute coefficients \Blue{$\gamma_j$} in \eqref{eq:trigpcomp}, see \eqref{tip:FM}
c = fft(exp(2*pi*i*(n/N)*(0:2*n)).*y)/N;
ch = [gamma,zeros(1,M-(2*n+1))]; % zero padding
v = conj(fft(conj(ch))); % Evaluate multiplication with conjugate Fourier matrix
q = exp(-2*pi*i*n*(0:M-1)/M).*v; % undo rescaling
