function [t,y,rej,ee] = odeintadapt_ext(Psilow,Psihigh,T,y0,h0,reltol,abstol,hmin)
t = 0; y = y0; h = h0; rej = [];  ee = 0;    % \label{odeintadapt:1}
while ((t(end) < T) && (h > hmin)) % \label{odeintadapt:2}
  yh = Psihigh(h,y0);  % \label{odeintadapt:3}
  yH = Psilow(h,y0);   % \label{odeintadapt:4}
  est = norm(yH-yh);   % \label{odeintadapt:5}
  
  if (est < max(reltol*norm(y0),abstol))                   % \label{odeintadapt:6}
%    fprintf('accept step, h = %f, est = %f. t = %f\n',h,est,t(end));
    y0 = yh; y = [y,y0]; t = [t,t(end) + min(T-t(end),h)]; % \label{odeintadapt:7}
    h = 1.1*h;  ee = [ee,est];                               % \label{odeintadapt:8}
  else
%    fprintf('reject step, h = %f, est = %f. t = %f\n',h,est,t(end));
    rej = [rej,t]; h = h/2; end                                       % \label{odeintadapt:9}
end
