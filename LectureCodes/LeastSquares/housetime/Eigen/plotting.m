delimiterIn = ' ';
headerlinesIn = 0;
r = importdata('housetime.dat', delimiterIn, headerlinesIn);

figure('name','housetiming');
loglog(r(:,1),r(:,2),'r+', r(:,1),r(:,3),'m*', r(:,1),r(:,4),'bo',...
       r(:,1),r(1,2)*r(:,1).^4/(r(1,1)^4),'k-',...
       r(:,1),r(1,2)*r(:,1).^6/(r(1,1)^6),'k--');
xlabel('{\bf n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
legend('HouseholderQR','qr\_decomp\_full()','qr\_decomp\_eco()',...
       'O(n^4)','O(n^6)','location','northwest');

print -depsc2 './housetime_m.eps';



