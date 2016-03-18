function y = fourcoeffcomp(c,m,ovsmpl)
% Compute the Fourier coefficients \Blue{$y_{-m},\ldots,y_{m}$} of the function
% \Blue{$c:[0,1[\mapsto \bbC$} using an \Magenta{oversampling factor} \texttt{ovsmpl}.
% \texttt{c} must be a handle to a function \texttt{\symbol{64}(t)}, which accepts row
% vector arguments
if (nargin < 3), ovsmpl = 2; else ovsmpl = ceil(ovsmpl); end
N = (2*m+1)*ovsmpl; h = 1/N; % Number of quadrature points
% (Inverse) discrete Fourier transform
y = ifft(c(0:h:1-h));
% Undo oversampling and wrapping of Fourier coefficient array
y = [y(N-m+1:N),y(1:m+1)];
