function [A,B] = factorize_X_AB(X,k)
% factorize X mxn of rank at most k in X=AB';  A mxk, Bnxk
MyTol = 1e-6;
[U,S,V] = svd(X,0);
if (k+1 <= min(size(X))) 
if (S(k+1,k+1) > MyTol)
warning ('The rank of the matrix is greater than the required rank');
end
end
A = U(:,1:min(k,end))*S(1:min(k,end),1:min(k,end));
B = V(:,1:min(k,end));