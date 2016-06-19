function x = solvearrow2(alpha,b,c,d,y)
z = b./d;
if (max(abs(z)) > 1/eps)
  error('Possible Instability!'); end;
xi = (y(1) - dot(z,y(2:end)))/(alpha-dot(z,c));
x = (y(2:end)-xi*c)./d;
x = [xi; x];
