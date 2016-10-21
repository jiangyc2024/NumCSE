function c = foursumnaive(signal,M,N)
% Approximate direct evaluation of Fourier sum according to the defining
% formula \eqref{eq:fourtrunc}, \texttt{signal} is a handle to a function of
% type \texttt{\symbol{64}(k)} providing the \Blue{$y_k$}, \texttt{M} specifies truncation
% of series according to \eqref{eq:fourtrunc}, \texttt{N} is the number of equidistant
% evaluation points for \Blue{$c$} in \Blue{$[0,1[$}. 

t = 0:1/N:1-1/N;  % Evaluation points for Fourier sum \Blue{$c$}
c = signal(0)*ones(1,N);
omega = exp(-2*pi*i*t);
omp = omega; omm = 1./omega;

% Inefficient direct summation of Fourier series
for k=1:M
  c = c+signal(k)*omp;
  c = c+signal(-k)*omm;
  omp = omp.*omega;
  omm = omm./omega;
end
