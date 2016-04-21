% 16.06.2009            horner.m
% evaluation of a polynomial in monomial representation using Horner scheme
% p is a line vector 1x(deg+1),
% x and y are matrices of the same size

function y = horner(p,x)
y = p(1);
for i=2:length(p)
    y = x.*y+p(i);
end
