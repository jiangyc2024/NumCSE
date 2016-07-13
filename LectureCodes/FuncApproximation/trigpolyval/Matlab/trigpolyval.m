function q = trigpolyval(t,y,x)
% Evaluation of trigonometric interpolant at numerous points
% \texttt{t}: row vector of nodes \Blue{$t_{0},\ldots,t_n\in[0,1[$}
% \texttt{y}: row vector of data \Blue{$y_{0},\ldots,y_{n}$}
% \texttt{x}: row vector of evaluation points \Blue{$x_{1},\ldots,x_{N}$}
N = length(y); if (mod(N,2)~=1), error('#pts odd required'); end
n = (N-1)/2; 
tc = exp(2*pi*i*t);     % Interpolation nodes on unit circle
z = exp(2*pi*i*n*t).*y; % Rescaled values, according to \Blue{$q(t) = e^{-2\pi int}\cdot p(e^{2\pi it})$}
% Evaluation of interpolating polynomial on unit circle, see Code~\ref{barycentricformula}
p = intpolyval(tc,z,exp(2*pi*i*x)); 
q = exp(-2*pi*i*n*x).*p; % Undo the scaling, see \eqref{eq:scale}
