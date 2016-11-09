function c = fftrec(y)
n = length(y);
if (n == 1), c = y; return;
else
  c1 = fftrec(y(1:2:n)); 
  c2 = fftrec(y(2:2:n));
  c = [c1;c1] + (exp(-2*pi*i/n).^((0:n-1)')) .*[c2;c2];
end
