% Example: Convergence of Broydens method
F = @(x) [x(1)^2-x(2)^4; x(1) - x(2)^3];
DF = @(x) [2*x(1),-4*x(2)^3; 1 , -3*x(2)^2];
x0 = [0.7;0.7];

[x,res] = broyden(F,x0,DF(x0),0.0000001,20,[1;1]);
[x,res] = updbroyd(F,x0,DF(x0),0.0000001,20,[1;1]);
% return;

[x,nres] = newton(x0,F,DF,size(res,1),[1;1]);
[x,sres] = simplenewt(x0,F,DF,size(res,1),[1;1]);
figure('name','Convergence of Broydens method');
semilogy(res(:,1),res(:,3),'m.',...
	 res(:,1),res(:,4),'r-*',...
	 nres(:,1),nres(:,3),'k.',...
	 nres(:,1),nres(:,4),'b-*',...
	 sres(:,1),sres(:,4),'c--*');
axis([0 size(res,1) 1.0e-15 10]);
set(gca,'fontsize',12);
grid on;
xlabel('Iterationsschritt');
ylabel('Normen');
legend('Broyden: ||F(x^{(k)})||','Broyden: Fehlernorm',...
       'Newton: ||F(x^{(k)})||','Newton: Fehlernorm',...
       'Newton (vereinfacht)',4);
print -dpsc2 'broydencvg.eps';

figure('name','Convergence monitor');
res(:,5) = [res(2:end,5);res(end,5)];
[ax,h1,h2] = plotyy(res(:,1),res(:,4),res(:,1),res(:,5),'semilogy','semilogy');
hold on;
set(ax(1),'fontsize',12);
set(ax(1),'ycolor',[0 0.75 0]);
set(ax(1),'ylim',[1.0e-12 2]);
set(ax(1),'ytick',10.^[-10 -8 -6 -4 -2 0]);
set(ax(2),'fontsize',12);
set(ax(2),'xcolor',[ 0 0 0 ]);
set(ax(2),'ycolor',[ 1 0 0 ]);
set(ax(2),'ygrid','on');
set(ax(2),'ytick',[0.01 0.1 1 10]);
set(h1,'linestyle','-','color',[0 0.75 0],'linewidth',1,'marker','*');
set(h2,'linestyle','-','color','r','linewidth',1,'marker','p');
set(get(ax(1),'ylabel'),'string','Fehlernorm','color',[0 0 0],'fontsize',14);
set(get(ax(1),'xlabel'),'string','Iterationsschritt','color',[0 0 0],'fontsize',14);
set(get(ax(2),'ylabel'),'string','Konvergenzmonitor','color',[0 0 0],'fontsize',14);
%plot(ax(2),res(:,1),ones(size(res,1),1),'r--');
hold off;

print -dpsc2 'broydencvgmon.eps';

% F = @(x) x*exp(x)-1;
% DF = @(x) exp(x)*(1+x);
% J0 = (F(5)-F(0))/5;
% x0 = 5;
% x = broyden(F,x0,J0,0.000001);
% x = secant(F,x0,J0,0.000001);
