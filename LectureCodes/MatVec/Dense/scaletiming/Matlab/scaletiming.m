% MATLAB script for timing a smart and foolish way to carry out multiplication
% with a scaling matrix in MATLAB, see Rem.~\ref{rem:scaling}.
nruns = 3; timings = [];
for n=2.^(2:14)
  d = rand(n,1); x = rand(n,1);
  tbad = realmax; tgood = realmax;
  for j=1:nruns
    tic; y = diag(d)*x; tbad = min(tbad,toc); 
    tic; y = d.*x; tgood = min(tgood,toc);  
  end
  timings = [timings; n, tgood, tbad]; 
end
    
figure('name','scaletimings');
loglog(timings(:,1),timings(:,2),'r*',...
       timings(:,1),timings(:,3),'b+');
xlabel('{\bf vector length n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
title('Timings for different ways to do scaling');
legend('using .*','using diag','location','best');

print -depsc2 '../PICTURES/scaletiming.eps';
