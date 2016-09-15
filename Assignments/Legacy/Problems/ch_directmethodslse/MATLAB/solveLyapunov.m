clear
n=5;
I=eye(n);
A= [10 2 3 4 5;
    6 20 8 9 1;
    1 2 30 4 5;
    6 7 8 20 0; 
    1 2 3 4 10];
C=buildC(A)
b=I(:);
vecX=C\b;
X=reshape(vecX',n,n)
error=norm(A*X+X*A'-I)
