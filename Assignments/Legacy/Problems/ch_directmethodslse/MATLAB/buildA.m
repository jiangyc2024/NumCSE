clear
n=7;
a=sym('a',[n-1 1]);
b=sym('b',[n-2 1]);
N=8;
C=100000;
deter=zeros(N-1,1)+realmax;
parfor n=2:N
for i=1:C
    a=rand([n-1 1]);
    b=rand([n-2 1]);

A=zeros(n,n);

for i=1:n
    A(i,i)=2;
end

for i=1:n-1
    A(i,i+1)=a(i);
end

for i=1:n-2
    A(i+2,i)=b(i);
end
deter(n-1)=min(deter(n-1),det(A));
end
end
deter
