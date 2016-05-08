function x = gn(x,F,J,tol)
s = J(x)\F(x);                %\label{gn:2}
x = x-s;
while (norm(s) > tol*norm(x)) %\label{gn:term}
  s = J(x)\F(x);              %\label{gn:5} 
  x = x-s;
end
