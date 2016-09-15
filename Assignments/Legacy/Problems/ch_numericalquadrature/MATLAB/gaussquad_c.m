% Input:  p, p-Punkt-Quadratur, exakt bis Polynom (2p-1)-ten Grades
% Output: x, Nullstellen dh Stuetzstellen
%         w, Gewichte
function [x,w]=gaussQuad(p)
b = zeros(p-1,1);
for i=1:(p-1)
   b(i)=i/sqrt(4*i*i-1);
end
J=diag(b,-1)+diag(b,1);
[ev,ew]=eig(J);
for i=1:p 
  ev(:,i) = ev(:,i)./norm(ev(:,i)); 
end
x=diag(ew);
w=(2*(ev(1,:).*ev(1,:)))';
