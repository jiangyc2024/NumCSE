function res = numquad(f,a,b,N,mode)
% Numerical quadrature on $[a,b]$ by polynomial quadrature formula
% f -> function to be integrated (handle), must support vector arguments
% a,b -> integration interval $[a,b]$ (endpoints included)
% N -> Maximal degree of polynomial
% mode (equidistant, Chebychev = Clenshaw-Curtis, Gauss) selects quadrature rule
if (nargin < 5), mode = 'equidistant'; end
res = [];

if strcmp(mode,'Gauss')
  for deg=1:N
    [gx,w] = gaussquad(deg);
    % Gauss points for \Blue{$[a,b]$}
    x = 0.5*(b-a)*gx+0.5*(a+b);
    y = feval(f,x);
    res = [res; deg, 0.5*(b-a)*dot(w,y)];
  end
else
  p = (N+1:-1:1);
  w = (b.^p - a.^p)./p;
  w
  sum(w)
  for deg=1:N
    if strcmp(mode,'Chebychev')
      % Chebychev nodes on \Blue{$[a,b]$}, see \eqref{eq:CHEBNODES}
      x = 0.5*(b-a)*cos((2*(0:deg)+1)/(2*deg+2)*pi)+0.5*(a+b);
      x
    else
      x = (a:(b-a)/deg:b);
      x
    end
    % ``Quick and dirty'' implementation through polynomial interpolation
    y = feval(f,x);
    poly = polyfit(x,y,deg);
    poly
    res = [res; deg, dot(w(N+1-deg:N+1),poly)];
  end
end
