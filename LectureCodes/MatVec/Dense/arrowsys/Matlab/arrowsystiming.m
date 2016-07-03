t = []; K = 3;
for l=3:12
    n = 2^l; alpha = 2;
    b = ones(n,1); c = (1:n)';
    d = -ones(n,1); 
    y = (-1).^(1:(n+1))'; 
    t1 = realmax; t2 = realmax;
    for k=1:K    
      tic; x1 = arrowsys_slow(d,c,b,alpha,y); t1 = min(t1,toc); 
      tic; x2 = arrowsys_fast(d,c,b,alpha,y); t2 = min(t2,toc); 
    end  
    t = [t; n t1 t2];
end
loglog(t(:,1),t(:,2),'b-*',t(:,1),t(:,3),'r-+');
xlabel('{\bf n}','fontsize',14);
ylabel('{\bf tic-toc runtime [s]}','fontsize',14);
legend('arrowsys_slow','arrowsys_fast','location','best');
print -depsc '../PICTURES/arrowsystiming.eps';