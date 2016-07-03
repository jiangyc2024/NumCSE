function rescirc(R,L,C)
% Ex.~\ref{ex:networlewp}: Numerical nodal analysis of the resonant circuit
% $R$, $L$, $C$ $\hat{=}$ network component parameters

Z = 1/R; K = 1/L;

% Matrices \Blue{$\VW$}, \Blue{$\VC$}, \Blue{$\VS$} for nodal analysis of circuit
Wmat = [0 0 0; 0 Z 0; 0 0 0];
Cmat = [C 0 0; 0 C 0; 0 0 C];
Smat = [K -K 0; -K 2*K -K; 0 -K K];
% System matrix from nodal analysis
Amat = @(w) (Wmat+i*w*Cmat+Smat/(i*w));

% Scanning source currents
res = [];
for w=0.01:0.01:2
  res = [res; w, abs(Amat(w)\[C;0;0])'];
end

figure('name','resonant circuit');
plot(res(:,1),res(:,2),'r-',res(:,1),res(:,3),'m-',res(:,1),res(:,4),'b-');
xlabel('{\bf angular frequency \omega of source voltage U}','fontsize',14);
ylabel('{\bf maximum nodal potential}','fontsize',14);
title(sprintf('R = %d, C= %d, L= %d',R,L,C));
legend('|u_1|','|u_2|','|u_3|');

print -depsc2 '../PICTURES/rescircpot.eps'

% Solving generalized eigenvalue problem \eqref{eq:nwevp1}
Zmat = zeros(3,3); Imat = eye(3,3);
% Assemble \Blue{$6\times 6$}-matrices \Blue{$\VM$} and \Blue{$\VB$}
Mmat = [Wmat,Smat; Imat, Zmat];
Bmat = [-i*Cmat, Zmat; Zmat , i*Imat];
% Solve \emph{generalized eigenvalue problem}, \emph{cf.} \eqref{eq:nevpgen}
omega = eig(Mmat,Bmat);

figure('name','resonances');
plot(real(omega),imag(omega),'r*'); hold on;
ax = axis;
plot([ax(1) ax(2)],[0 0],'k-');
plot([ 0 0],[ax(3) ax(4)],'k-');
grid on;
xlabel('{\bf Re(\omega)}','fontsize',14);
ylabel('{\bf Im(\omega)}','fontsize',14);
title(sprintf('R = %d, C= %d, L= %d',R,L,C));
legend('\omega');

print -depsc2 '../PICTURES/rescircomega.eps'
