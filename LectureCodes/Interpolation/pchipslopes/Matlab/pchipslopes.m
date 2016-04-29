function c = pchipslopes(t,y)
% Calculation of local slopes \Blue{$c_i$} for shape preserving cubic Hermite interpolation, see \eqref{mteq:lim}, \eqref{mteq:hm}
% \texttt{t}, \texttt{y} are row vectors passing the data points
n = length(t); h = diff(t); delta = diff(y)./h; % linear slopes
c = zeros(size(h));
k = find(sign(delta(1:n-2)).*sign(delta(2:n-1))>0)+1;
% Compute reconstruction slope according to \eqref{mteq:hm}
w1 = 2*h(k)+h(k-1); w2 = h(k)+2*h(k-1);
c(k) = (w1+w2)./(w1./delta(k-1) + w2./delta(k));
%  Special slopes at endpoints, beyond \eqref{mteq:hm}
c(1) = pchipend(h(1),h(2),delta(1),delta(2));
c(n) = pchipend(h(n-1),h(n-2),delta(n-1),delta(n-2));

function d = pchipend(h1,h2,del1,del2)
%  Noncentered, shape-preserving, three-point formula.
d = ((2*h1+h2)*del1 - h1*del2)/(h1+h2);
if sign(d) ~= sign(del1), d = 0;
elseif (sign(del1)~=sign(del2))&(abs(d)>abs(3*del1)), d = 3*del1; end 
