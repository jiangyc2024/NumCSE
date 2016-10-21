function c = foursum(signal,M,N)
% Approximate evaluation of Fourier series, \texttt{signal} is a handle to a function of
% type \texttt{\symbol{64}(k)} providing the \Blue{$y_k$}, \texttt{M} specifies truncation
% of series according to \eqref{eq:fourtrunc}, \texttt{N} is the number of equidistant
% evaluation points for \Blue{$c$} in \Blue{$[0,1[$}. 

y = signal(-M:M); % Sample signal from \Blue{$-M$} to \Blue{$M$}
m = 2*M+1;        % Length of signal
% Ensure that there are more sampling points than terms in series
if (m > N), l = ceil(m/N); N = l*N; else l = 1; end

% \Hyperlink{ZEROPAD}{Zero padding} and wrapping of signal, see Code~\ref{mc:freqfilter}
y_ext = zeros(1,N); y_ext(1:M+1) = y(M+1:end); y_ext(N-M+1:N) = y(1:M);

% Perform DFT and decimate output vector
c = fft(y_ext); c = c(1:l:end);
  
