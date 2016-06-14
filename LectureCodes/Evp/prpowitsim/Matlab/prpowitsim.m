function prpowitsim(d,Nsteps)
% MATLAB way of specifying Default arguments 
if (nargin < 2), Nsteps = 5; end
if (nargin < 1), d = 0.15; end
% load connectivity matrix and build transition matrix
load harvard500.mat; A = prbuildA(G,d); 
N = size(A,1); x = ones(N,1)/N; 

figure('position',[0 0 1200 1000]);
plot(1:N,x,'r+'); axis([0 N+1 0 0.1]); 
% Plain power iteration for stochastic matrix \Blue{$\VA$}
for l=1:Nsteps
  pause; x = A*x; plot(1:N,x,'r+'); axis([0 N+1 0 0.1]); 
  title(sprintf('{\\bf step %d}',l),'fontsize',14);
  xlabel('{\bf harvard500: no. of page}','fontsize',14);
  ylabel('{\bf page rank}','fontsize',14); drawnow;
end
