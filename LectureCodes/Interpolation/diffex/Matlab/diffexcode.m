function d = diffex(f,x,h0,atol,rtol)
% \texttt{f}: handle of to a function defined in a neighborhood of \Blue{$x\in\mathbb{R}$},
% \texttt{x}: point at which approximate derivative is desired,
% \texttt{h0}: initial distance from \texttt{x},
% \texttt{tol}: relative target tolerance
h = h0;
% Aitken-Neville scheme, see Code~\ref{AitkenNeville} (\Blue{$x=0$}!)
y(1) = (f(x+h0)-f(x-h0))/(2*h0);
for i=2:10
  h(i) = h(i-1)/2;
  y(i) = f(x+h(i))-f(x-h(i)))/h(i-1);
  for k=i-1:-1:1
    y(k) = y(k+1)-(y(k+1)-y(k))*h(i)/(h(i)-h(k));
  end
  % termination of extrapolation, when desired tolerance is achieved
  errest = abs(y(2)-y(1)); % error indicator
  if ((errest < rtol*abs(y(1))) || (errest < atol)), break; end % \label{de:1}
end
d = y(1);
