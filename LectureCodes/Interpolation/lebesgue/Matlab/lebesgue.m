function l = lebesgue(t,N)
% Computation of \Hyperlink{LEBESGUE}{Lebesgue constant} of polynomial interpolation with
% knots \Blue{$t_i$} passed in the row vector \texttt{t} based on
% \eqref{eq:IPN1}. \texttt{N} specifies the number of sampling points for the approximate
% computation of the maximum norm of the Lagrange polynomial on the interval
% \Blue{$[-1,1]$}.
n = length(t);
den = []; % denominators of normalized Lagrange polynomials relative to the nodes t
for i=1:n
  den = [den; prod(t(i)-t([1:i-1,i+1:end]))];
end

% Default argument
if (nargin < 2), N = 1E5; end
l = 0; % return value
for x=-1:2/N:1 % sampling points for approximate computation of \Blue{$\NxLinf{\cdot}{[-1,1]}$}
  s = 0;
  for i = 1:n
    % \texttt{v} provides value of the normalized Lagrange polynomials
    v = prod(x-t([1:i-1,i+1:end]))/den(i);
    s = s+abs(v); % sum over the modulus of the polynomials 
  end
  l = max(l,s);   % maximum of sampled values
end
