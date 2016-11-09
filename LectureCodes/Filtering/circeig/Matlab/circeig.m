function circeig
% Computing the eigenvectors of a \Magenta{random} circulant matrix
n = 8;
% The MATLAB way to generate a circulant matrix
C = gallery('circul',rand(n,1)); [V1,D1] = eig(C);

for j=1:n
  figure; bar(1:n,[real(V1(:,j)),imag(V1(:,j))],1,'grouped');
  title(sprintf('Circulant matrix 1, eigenvector %d',j));
  xlabel('{\bf vector component index}','fontsize',14);
  ylabel('{\bf vector component value}','fontsize',14);   
  legend('real part','imaginary part','location','southwest'); 
  print('-depsc2',sprintf('../PICTURES/circeig1ev%d.eps',j));
end

C = gallery('circul',rand(n,1)); [V2,D2] = eig(C);

for j=1:n
  figure; bar(1:n,[real(V2(:,j)),imag(V2(:,j))],1,'grouped');
  title(sprintf('Circulant matrix 2, eigenvector %d',j));
  xlabel('{\bf vector component index}','fontsize',14);
  ylabel('{\bf vector component value}','fontsize',14);   
  legend('real part','imaginary part','location','southwest'); 
  print('-depsc2',sprintf('../PICTURES/circeig2ev%d.eps',j));
end

figure; plot(1:n,real(diag(D1)),'r+',1:n,imag(diag(D1)),'b+',...
	     1:n,real(diag(D2)),'m*',1:n,imag(diag(D2)),'k*');	
ax = axis; axis([0 n+1 ax(3) ax(4)]);
xlabel('{\bf index of eigenvalue}','fontsize',14);
ylabel('{\bf eigenvalue}','fontsize',14);   
legend('C_1: real(ev)', 'C_1: imag(ev)', 'C_2: real(ev)', 'C_2: imag(ev)', 'location','northeast'); 

print -depsc2 '../PICTURES/circeigev.eps';

