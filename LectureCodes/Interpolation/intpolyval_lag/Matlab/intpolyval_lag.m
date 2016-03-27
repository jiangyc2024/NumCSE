% 16.06.2009            intpolyval_leg.m
% Evaluation of the interpolation polynomials with Lagrange polynomials
%
% input:    t   nodes
%           y   values in t
%           x   evaluation points (scalar vector or matrix)

function p = intpolyval_lag(t,y,x)
p=zeros(size(x));

for k=1:length(t);
    p=p + y(k)*lagrangepoly(x, k-1, t);
end
    