function [t,y] = odeintssctrl(Psilow,p,Psihigh,T,y0,h0,reltol,abstol,hmin)
t = 0; y = y0; h = h0;             % \label{ssctrl:1}
while ((t(end) < T) && (h > hmin)) % \label{ssctrl:2}
  yh = Psihigh(h,y0);  % high order discrete evolution \Blue{$\widetilde{\Psibf}^h$}\label{ssctrl:3}
  yH = Psilow(h,y0);   % low order discrete evolution \Blue{${\Psibf}^h$}\label{ssctrl:4}
  est = norm(yH-yh);   % $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\label{ssctrl:5}
  
  tol = max(reltol*norm(y(:,end)),abstol); % \label{ssctrl:6a}
  h = h*max(0.5,min(2,(tol/est)^(1/(p+1)))); % Optimal stepsize according to \eqref{eq:ssc}\label{ssctrl:6b}
  if (est < tol)                                           % \label{ssctrl:6}
    y0 = yh; y = [y,y0]; t = [t,t(end) + min(T-t(end),h)]; % step \Magenta{accepted}\label{ssctrl:7}
  end
end

