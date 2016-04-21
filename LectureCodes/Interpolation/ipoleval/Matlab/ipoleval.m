% 16.06.2009            ipoleval.m
% Evaluation of the interpolation polynomials with Matlab polyfit+polyval
%
% input:    t   nodes
%           y   values in t
%           x   evaluation points


function v=ipoleval(t,y,x)
p = polyfit(t,y,length(y)-1);
v=polyval(p,x);
