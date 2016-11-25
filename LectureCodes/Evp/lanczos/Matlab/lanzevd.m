% ex: lanzcosev
n=100; M=gallery('tridiag',-0.5*ones(n-1,1),2*ones(n,1),-1.5*ones(n-1,1));
[Q,R]=qr(M); A=Q'*diag(1:n)*Q; % eigenvalues \Blue{$1,2,3,\ldots,100$}
ev = sort(eig(A));

[V,res] = lanczosev(A,30,5,ones(100,1));
figure; semilogy(1:30,abs(res(1:30,end-0)-ev(end-0)),'r-+',...
		 2:30,abs(res(2:30,end-1)-ev(end-1)),'m-^',...
		 3:30,abs(res(3:30,end-2)-ev(end-2)),'b-o',...
		 4:30,abs(res(4:30,end-3)-ev(end-3)),'k-*');
xlabel('{\bf step of Lanzcos process}','Fontsize',14);
ylabel('{\bf |Ritz value-eigenvalue|}','fontsize',14);
legend('\lambda_n','\lambda_{n-1}','\lambda_{n-2}','\lambda_{n-3}',3);

print -depsc2 '../PICTURES/lanczosev1.eps';
clear all;

A = gallery('minij',100); ev = sort(eig(A));
[V,res] = lanczosev(A,15,5,ones(100,1));
figure; semilogy(1:15,abs(res(1:15,end-0)-ev(end-0)),'r-+',...
		 2:15,abs(res(2:15,end-1)-ev(end-1)),'m-^',...
		 3:15,abs(res(3:15,end-2)-ev(end-2)),'b-o',...
		 4:15,abs(res(4:15,end-3)-ev(end-3)),'k-*');
xlabel('{\bf step of Lanzcos process}','Fontsize',14);
ylabel('{\bf error in Ritz values}','fontsize',14);
legend('\lambda_n','\lambda_{n-1}','\lambda_{n-2}','\lambda_{n-3}',3);

print -depsc2 '../PICTURES/lanczosev.eps';
