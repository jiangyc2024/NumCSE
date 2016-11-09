% Visualize limit \Blue{$m\to\infty$} for a \Blue{$2m+!$}-periodic signal and
% its discrete Fourier transform ``squeezed'' into \Blue{$[0,1]$}.

% range of plot for visualization of discrete signal
Npow = 8; N = 3*(2^Npow+1);
% function defining discrete signal
yfn = @(k) 1./(1+k.*k);
% loop over different periods \Blue{$2^{l}+1$}
for mpow = [4 5 6 7]
  m = 2^mpow;
  ybas = yfn([(-m:-1),0,(1:m)]);
  Ncp = floor(N/(2*m+1));
  y = repmat(ybas,1,Ncp);
  
  figure; hy = stem((1:length(y))-(length(y)+1)/2,y,'r'); 
  ax = axis; axis([ax(1) ax(2) 0 1.1]); hold on;
  set(hy,'markersize',0);
  xlabel('{\bf index i of sampling instance}','fontsize',14);
  ylabel('{\bf y_{i}}','fontsize',14);
  title(sprintf('period of signal (y_{i}) = %d',2*m+1));
  print('-depsc2',sprintf('../PICTURES/persig%d.eps',mpow));
  
  % DFT of wrapped signal (one period)
  c = fft([ybas(m+1:end),ybas(1:m)]);
  hc = 1/(2*m+1);
  figure; hc = stem(0:hc:1-hc,c,'m');
  set(hc,'markersize',0);  
  xlabel('{\bf t_{k}}','fontsize',14);
  ylabel('{\bf c(t_{k})}','fontsize',14);
  title(sprintf('DFT of periodic signal (y_{i}) with period %d',2*m+1));
  print('-depsc2',sprintf('../PICTURES/persigdft%d.eps',mpow));
end
