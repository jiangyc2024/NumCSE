function [U, S, V] = svd_ab(A,B)     % A,B:   n x k
[QA, RA] = qr(A,0);                  % QA,QB: n x k; RA,RB k x k
[QB, RB] = qr(B,0);
[U, S, V] = svd(RA*RB');             % U,S,V: k x k
U = QA*U;                            % U,V:   n x k
V = QB*V;