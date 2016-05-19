function chebpolplot(nmax)
% plots Chebychev polybomials up to degree \texttt{nmax} on \Blue{$[-1,1]$}
N = 1000; x = -1:2/N:1;
V = chebpolmult(nmax,x); % compute values of Chebychev polynomials
for j=0:nmax, leg{j} = sprintf('n=%i',j); end
plot(x,V','-');
legend(leg);
xlabel('{\bf t}','Fontsize',14);
ylabel('{\bf T_n(t)}','Fontsize',14);
axis([-1.1 1.1 -1.1 1.1]);
