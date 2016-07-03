function prstochsim(Nhops)
% Load web graph data stored in \Blue{$N\times N$}-matrix \Blue{$\VG$}
load harvard500.mat;
N = size(G,1); d = 0.15;
count = zeros(1,N); cp = 1;
figure('position',[0 0 1200 1000]); pause;
for n=1:Nhops
  % Find links from current page \texttt{cp}
  idx = find(G(:,cp)); l = size(idx,1); rn = rand(); % \label{prst:1}
  % If no links, jump to any other pages with equal probability
  if (isempty(idx)),  cp = floor(rn*N)+1;
  % With probabilty \Blue{$d$} jump to any other page  
  elseif (rn < d), cp = floor(rn/d*N)+1;
  % Follow outgoing links with equal probabilty
  else cp = idx(floor((rn-d)/(1-d)*l)+1,1);
  end
  count(cp) = count(cp) + 1;
  plot(1:N,count/n,'r.'); axis([0 N+1 0 0.1]); 
  xlabel('{\bf harvard500: no. of page}','fontsize',14);
  ylabel('{\bf page rank}','fontsize',14);
  title(sprintf('{\\bf page rank, harvard500: %d hops}',n),'fontsize',14);
  drawnow;
end
