
% Problem 5d)-e)
% INPUT
% X - original, non-transformed, points of the image
% P - data about a figure from the tree
% OUTPUT
% A - approximation of the matrix A
% err - error in approximation
function [A, err] = complinmap(X, P)

% Compute the matrix of the overdetermined system
B = shapeidentmat(X);
% Adjust P so that according to the system
P = reshape(P, 30, 1);

% Find the solution of the said system and adjust its form
A = reshape( B\P, 2, 2)';

% Compute the error
err = norm(reshape(P, 2, 15) - A*X);