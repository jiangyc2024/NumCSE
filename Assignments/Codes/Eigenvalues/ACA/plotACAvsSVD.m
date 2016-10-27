function plotACAvsSVD ()
n = 1000;
[J I] = meshgrid(1:n,1:n);
A1 = 1 ./ (1 + (I - J).^2 / n.^2 );
A2 = exp ( - abs(I - J) / n );
errT1 = []; errT2 = [];
errS1 = []; errS2 = [];
B1 = ACA(A1,1e-14);
B2 = ACA(A2,5e-2);
Ls1 = 1:size(B1,2);
Ls2 = 1:size(B2,2);
[U1,S1,V1] = svd(A1);
[U2,S2,V2] = svd(A2);
for L = Ls1
    TL = computeT(B1,L);
    SL = U1(:,1:L) * S1(1:L,1:L) * V1(:,1:L)';
    errT1 = [errT1 norm(A1 - TL)];
    errS1 = [errS1 norm(A1 - SL)];
end
for L = Ls2
    TL = computeT(B2,L);
    SL = U2(:,1:L) * S2(1:L,1:L) * V2(:,1:L)';
    errT2 = [errT2 norm(A2 - TL)];
    errS2 = [errS2 norm(A2 - SL)];
end
semilogy(Ls1,errT1,'r-x',Ls1,errS1,'b-o', ...
         Ls2,errT2,'g-x',Ls2,errS2,'m-o','Linewidth',2,'MarkerSize',8);
title('{\bf approximation error for ACA and SVD}');
xlabel('{\bf L}');
ylabel('{\bf 2-norm of the approximation error}');
legend('ACA1','SVD1','ACA2','SVD2');
print -depsc '../PICTURES/plotACAvsSVD.eps'
