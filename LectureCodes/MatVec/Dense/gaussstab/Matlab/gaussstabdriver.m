% Driver for gaussstab

gaussstab;

figure('name','gaussstab');
loglog(result(:,1),result(:,2),'r+-',...
        result(:,1),result(:,3),'b-*');
xlabel('{\bf \epsilon}','fontsize',14);
legend('relative error','relative residual');

print -depsc2 '../PICTURES/gaussstab.eps';

figure('name','gaussstabcond');
[ax,h1,h2] = plotyy(result(:,1),result(:,4),result(:,1),result(:,2),@loglog);
set(h1,'linestyle','-','color','m','marker','*');
set(h2,'linestyle','-','color','r','marker','+');
set(get(ax(1),'Ylabel'),'String','{\bf cond(A)}','fontsize',14); 
set(get(ax(2),'Ylabel'),'String','relative error','fontsize',14); 
xlabel('{\bf \epsilon}','fontsize',14);
ylabel('{\bf cond(A)}','fontsize',14);
legend(h1,'cond(A)','location','northeast')
legend(h2,'relative error','location','southwest');

print -depsc2 '../PICTURES/gaussstabcond.eps';



