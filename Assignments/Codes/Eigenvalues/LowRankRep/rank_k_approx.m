function [Az,Bz] = rank_k_approx(Ax,Ay,Bx,By)
% Az*Bz' = best rank-k approx of AxBx'+AyBy'
% {Ax,Ay,AZ} m x k;   {Bx,By,BZ} n x k;   

[U,S,V] = svd_ab([Ax,Ay],[Bx,By]);  
% U: m x 2k;   S: 2k x 2k;   V: n x 2k 
k = size(Ax,2);
Az = U(:,1:k) * S(1:k,1:k);
Bz = V(:,1:k);