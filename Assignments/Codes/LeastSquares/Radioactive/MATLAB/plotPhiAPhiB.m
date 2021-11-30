% MATLAB script for plotting analytic solution of sub-problem 1 of HW problem "Radioactive"
l1 = 1;
l2 = 0.5;
PhiA = @(t,A0,B0) (A0*exp(-l1*t));
PhiB = @(t,A0,B0) (B0*exp(-l2*t) + A0*l1/(l2-l1)*(exp(-l1*t)-exp(-l2*t)));
% Time range and initial values
t = 0:0.01:5;
A0 = 1.0;
B0 = 0.0;

figure('name','Radioactive decay');
plot(t,PhiA(t,A0,B0),'b-',t,PhiB(t,A0,B0),'m-');
xlabel('{\bf time t}');
ylabel('{\bf amount of substances}');
title('Radioactive decay chain, \lambda_1 = 1, \lambda_2 = 0.5');
legend('\Phi_{A}','\Phi_{B}','location','best');
print -depsc2 'PhiApHiB.eps';
