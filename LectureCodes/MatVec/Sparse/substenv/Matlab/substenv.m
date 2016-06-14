function y = substenv(L,y,mc)
% evelope aware forward substitution for \Blue{$\VL\Vx=\Vy$}
% (\Blue{$\VL$} = lower triangular matrix)
% argument \texttt{mc}: column bandwidth vector
n = size(L,1); y(1) = y(1)/L(1,1);
for i=2:n
  if (mr(i) > 0)
    zeta = L(i,i-mr(i):i-1)*y(i-mr(i):i-1);
    y(i) = (y(i) - zeta)/L(i,i);
  else y(i) = y(i)/L(i,i); end
end

