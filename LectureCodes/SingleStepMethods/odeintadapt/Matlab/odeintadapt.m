function [t,y] = odeintadapt(Psilow,Psihigh,T,y0,h0,reltol,abstol,hmin)
t = 0; y = y0; h = h0;             % \label{odeintadapt:1}
while ((t(end) < T) && (h > hmin)) % \label{odeintadapt:2}
  yh = Psihigh(h,y0);  % high order discrete evolution \Blue{$\widetilde{\Psibf}^h$} \label{odeintadapt:3}
  yH = Psilow(h,y0);   % low order discrete evolution \Blue{${\Psibf}^h$} \label{odeintadapt:4}
  est = norm(yH-yh);   % $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\label{odeintadapt:5}
  
  if (est < max(reltol*norm(y0),abstol))                   % \label{odeintadapt:6}
    y0 = yh; y = [y,y0]; t = [t,t(end) + min(T-t(end),h)]; % \label{odeintadapt:7}
    h = 1.1*h;  % step \Magenta{accepted}, try with increased stepsize \label{odeintadapt:8}
  else, h = h/2; end   % step \Magenta{rejected}, try with half the stepsize \label{odeintadapt:9}
end
