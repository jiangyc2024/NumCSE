function trusssim(pos,top,mode)
% Visualizing nonzero eigenmodes of an elastic truss
n = size(pos,1); mode = min(n,mode);
% Compute eigenmodes
[l,V] = trussvib(pos,top);

% Plot resonances
figure('name','truss resonances');
plot(1:2*n,l,'r+'); hold on;
plot(1:3,[0 0 0],'mo');
xlabel('{\bf no. of eigenvalue}','fontsize',14);
ylabel('{\bf eigenvalue}','fontsize',14);
title('{\bf truss resonant frequencies}','fontsize',14);
print -depsc2 '../PICTURES/trussfreq.eps';

figure('name','truss mode');
trussplot(pos,top,'mo','k--');
u = [V(1:2:2*n-1,mode),V(2:2:2*n,mode)];
trussplot(pos+0.5*u,top,'r*','b-');
title(sprintf('mode %i: frequency = %d',mode,l(mode)),'fontsize',14);

print('-depsc2',sprintf('../PICTURES/trussmode%i.eps',mode));

