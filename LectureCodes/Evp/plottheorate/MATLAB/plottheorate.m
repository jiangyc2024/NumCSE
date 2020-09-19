function plottheorate
% Python: ../PYTHON/ConjugateGradient.py plottheorate()

[X,Y] = meshgrid(1:10,1:100); R = zeros(100,10);
for I=1:100
  t = 1/I;
  for j=1:10
    R(I,j) = 2*(1-t)^j/((1+t)^(2*j)+(1-t)^(2*j));
  end
end

figure; view([-45,28]); mesh(X,Y,R); colormap hsv;
xlabel('{\bf CG step l}','Fontsize',14);
ylabel('{\bf \kappa(A)^{1/2}}','Fontsize',14);
zlabel('{\bf error reduction (energy norm)}','Fontsize',14);

print -depsc2 '../PICTURES/theorate1.eps';

figure; [C,h] = contour(X,Y,R); clabel(C,h);
xlabel('{\bf CG step l}','Fontsize',14);
ylabel('{\bf \kappa(A)^{1/2}}','Fontsize',14);

print -depsc2 '../PICTURES/theorate2.eps';