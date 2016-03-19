function I = adaptquad(f,M,rtol,abstol)
% adaptive numerical quadrature: \texttt{f} is a function handle to integrand
h = diff(M);                     % distances of quadrature nodes \label{aq:1}
mp = 0.5*(M(1:end-1)+M(2:end));  % positions of midpoints \label{aq:2}
fx = f(M); fm = f(mp);           % \label{aq:3}
trp_loc  = h.*(fx(1:end-1)+2*fm+fx(2:end))/4; % trapezoidal rule \eqref{eq:comptrap}\label{aq:4}
simp_loc = h.*(fx(1:end-1)+4*fm+fx(2:end))/6; % Simpson rule \eqref{eq:compsimp}\label{aq:5}
I = sum(simp_loc);               % Simpson approximation of integral value\label{aq:6}
est_loc  = abs(simp_loc -trp_loc);  % local error estimate \eqref{eq:aqest}\label{aq:7}
err_tot  = sum(est_loc);            % estimate for quadrature error\label{aq:8}
% Termination based on \eqref{eq:aqstop}
if ((err_tot > rtol*abs(I)) && (err_tot > abstol)) % \label{aq:term}
  refcells = find(est_loc > 0.9*sum(est_loc)/length(h));
  I = adaptquad(f,sort([M,mp(refcells)]),rtol,abstol); % \label{aq:10} 
end
