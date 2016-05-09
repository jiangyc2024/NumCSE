r = [];
for phi=pi/200:pi/100:pi/2
    A = [1,cos(phi); 0,sin(phi)];
    r = [r; phi, cond(A),cond(A,'inf')];
end
plot(r(:,1),r(:,2),'r-', r(:,1),r(:,3),'b--');
xlabel('{\bf angle of n_1, n_2}','fontsize',14);
ylabel('{\bf condition numbers}','fontsize',14);
legend('2-norm','max-norm');
print -depsc2 '../PICTURES/linesec.eps';