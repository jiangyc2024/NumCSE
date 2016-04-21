x=1.1; h=2.^[-1:-5:-36];
atanerr = abs(dirnumdiff(atan,x,h)-1/(1+x^2))*(1+x^2);
sqrterr = abs(dirnumdiff(sqrt,x,h)-1/(2*sqrt(x)))*(2*sqrt(x));
experr  = abs(dirnumdiff(exp,x,h)-exp(x))/exp(x);

function [df]=dirnumdiff(f,x,h)
df=(f(x+h)-f(x))./h;
end
