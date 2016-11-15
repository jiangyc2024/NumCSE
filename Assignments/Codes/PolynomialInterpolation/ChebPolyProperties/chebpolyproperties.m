clear

n=20;
X=[-1:.01:1];
Y=zeros(size(X));

alpha=[ 0.196209
 6.50521e-18
   -0.263812
           0
    0.177416
-7.37257e-17
   -0.119406
  5.0307e-17
   0.0805004
-2.29851e-17
  -0.0544753
-3.85976e-17
   0.0371663
-1.30104e-17
  -0.0258038
 1.62413e-16
     0.01857
-7.25114e-16
  -0.0143074
 6.12574e-17
    0.012334];

for j=0:n
    Y=Y+alpha(j+1)*chebyshevT(j,X);
end
plot(X,Y,X,1./((5*X).^2+1))
legend('Chebyshev best approximation','Exact function')