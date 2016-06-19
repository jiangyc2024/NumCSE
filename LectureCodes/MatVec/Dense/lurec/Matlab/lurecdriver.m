function [L,U] = lurecdriver(A)
A = lurec(A);
% post-processing: 
% extract \Blue{$\VL$} and \Blue{$\VU$}
U = triu(A);
L = tril(A,-1) + eye(size(A));
