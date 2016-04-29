% 16.06.2009            lagrangepoly.m
%
% direct evaluation in x of the index^th Lagrange polynomial relative to nodes
%
% lagrangepoly(nodes(j), k-1, nodes) = delta_{j,k}
%
% inputs:   x       in R
%           index   integer in [0, length(nodes)]
%           nodes   row vector of integers

function L=lagrangepoly(x, index, nodes)
L=1;
for j=[0:index-1, index+1:length(nodes)-1];
    L = L .* (x-nodes(j+1)) ./ ( nodes(index+1)-nodes(j+1) );
end
