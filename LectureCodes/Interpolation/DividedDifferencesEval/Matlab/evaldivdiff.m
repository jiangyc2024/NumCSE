% 16.06.2009            evaldivdiff.m
%
% Evaluation of polynomial in Newton form = divided differences
% Backward algorithm
% input:    t   nodes
%           y   values in t
%           x   evaluation points (scalar vector or matrix)

function p = evaldivdiff(t,y,x)
n = length(y)-1;

dd=divdiff(t,y);

p=dd(n+1);
for j=n:-1:1
    p=   ( x-t(j) ).*p + dd(j);
end

