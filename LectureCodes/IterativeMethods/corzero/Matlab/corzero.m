y  = [];
x = 0.4;
for i = 1:15
  x = x + (cos(x)+1)/sin(x); 
  y = [y,x]; 
end 
err = y - x; 
rate = err(2:15)./err(1:14);
