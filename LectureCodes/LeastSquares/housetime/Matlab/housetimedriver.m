housetime;

figure;
loglog(r(:,1),r(:,3),'r+', r(:,1),r(:,4),'m*', r(:,1),r(:,5),'bo',...
       r(:,1),r(1,3)*r(:,1).^4/(r(1,1)^4),'k-',...
       r(:,1),r(1,3)*r(:,1).^6/(r(1,1)^6),'k--');
xlabel('{\bf n}','fontsize',14);
ylabel('{\bf time [s]}','fontsize',14);
legend('[Q,R] = qr(A)','[Q,R] = qr(A,0)','R = qr(A)',...
       'O(n^4)','O(n^6)','location','northwest');

print -depsc2 '../PICTURES/housetime.eps';
