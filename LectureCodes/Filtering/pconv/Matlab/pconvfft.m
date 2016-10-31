function z=pconvfft(u,x)
% Implementation of \eqref{eq:circdiag}, \emph{cf.} Lemma~\ref{lem:zirkdiag}
z = ifft(fft(u).*fft(x));

