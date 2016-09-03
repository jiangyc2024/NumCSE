function q = trigipequidcomp(a,b,N)
% Efficient evaluation of trigonometric polynomial at equidistant points
% column vectors \texttt{a} and \texttt{b} pass coefficients \Blue{$\alpha_j$}, \Blue{$\beta_j$} in
% representation \eqref{eq:trigpreal} 
n = length(a)-1; if (N < (2*n-1)), error('N too small'); end;
gamma = transpose(0.5*[a(end:-1:2)+i*b(end:-1:1);...
                    2*a(1);a(2:end)-i*b(1:end)]);
ch = [gamma,zeros(1,N-(2*n+1))]; % zero padding
v = conj(fft(conj(ch))); % Multiplication with conjugate Fourier matrix
q = exp(-2*pi*i*n*(0:N-1)/N).*v; % undo rescaling
