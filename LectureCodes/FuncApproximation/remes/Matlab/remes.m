function c = remes(f,f1,a,b,d,tol)
% \texttt{f} is a handle to the function, \texttt{f1} to its derivative
% d = polynomial degree (positive integer) 
% a,b = interval boundaries
% returns coefficients of polynomial in monomial basis 
% (\matlab{} convention, see \cref{rem:polyMatlab}).

n = 8*d;                  % n = number of sampling points
xtab=[a:(b-a)/(n-1):b]';  % Points of sampling grid
ftab = feval(f,xtab);     % Function values at sampling points
fsupn = max(abs(ftab));   % Approximate supremum norm of \Blue{$f$}
f1tab = feval(f1,xtab);   % Derivative values at sampling points
% The vector xe stores the current guess for the alternants; initial
% guess is Chebychev alternants \eqref{remez:chebalt}.
h=pi/(d+1); xe=(a+b)/2 + (a-b)/2*cos(h*[0:d+1]'); 
fxe=feval(f,xe);

maxit = 10;
% Main iteration loop of Remez algorithm
for k=1:maxit 
 % Interpolation at \Blue{$d+2$} points xe with deviations \Blue{$\pm\delta$}
 % Algorithm uses monomial basis, which is not optimal
 V=vander(xe); A=[V(:,2:d+2), (-1).^[0:d+1]']; % \com{LSE}
 c=A\fxe;                    % Solve for coefficients of polynomial \Blue{$q$}
 c1=[d:-1:1]'.*c(1:d);       % Monomial doefficients of derivative \Blue{$q'$}
 % Find initial guesses for the inner extremes by sampling; track sign
 % changes of the derivative of the approximation error
 deltab = (polyval(c1,xtab) - f1tab);
 s=[deltab(1:n-1)].*[deltab(2:n)]; 
 ind=find(s<0); xx0=xtab(ind); % approximate zeros of e'
 nx = length(ind);             % number of approximate zeros
 % Too few extrema; bail out
 if (nx < d), error('Too few extrema'); end

 % \com{Secant method} to determine zeros of derivative 
 % of approximation error
 F0 = polyval(c1,xx0) - feval(f1,xx0);
 % Initial guess from shifted sampling points
 xx1=xx0+(b-a)/(2*n);
 F1 = polyval(c1,xx1) - feval(f1,xx1);
 % Main loop of secant method 
 while min(abs(F1)) > 1e-12,
   xx2=xx1-F1./(F1-F0).*(xx1-xx0); 
   xx0=xx1; xx1=xx2; F0=F1;
   F1=polyval(c1,xx1) - feval(f1,xx1);
 end;
 
 % Determine new approximation for alternants; store in xe
 % If too many zeros of the derivative \Blue{$(f-p)'$}
 % have been found, select those, where the deviation is maximal.
 if (nx == d)
   xe=[a;xx0;b];
 elseif (nx == d+1)
   xmin = min(xx0); xmax = max(xx0);
   if ((xmin - a) > (b-xmax)), xe = [a;xx0];
   else xe = [xx0;b]; end
 elseif (nx == d+2)
   xe = xx0; 
 else 
   fx = feval(f,xx0);
   del = abs(polyval(c(1:d+1),xx0) - fx);
   [dummy,ind] = sort(del);
   xe = xx0(ind(end-d-1:end));
 end

 % Deviation in sampling points and approximate alternants
 fxe=feval(f,xe);
 del = [polyval(c(1:d+1),a) - ftab(1);...
	polyval(c(1:d+1),xe) - fxe;...
	polyval(c(1:d+1),b) - ftab(end)];
 % Approximation of supremum norm of approximation error
 dev = max(abs(del));
 % \com{Termination} of Remez iteration
 if ( dev < tol*fsupn), break; end
end;
c = c(1:d+1);
