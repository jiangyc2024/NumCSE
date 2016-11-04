function [C,idx] = pointcluster(X,n)
% n-quantization of point set in k-dimensional space based on minimizing the mean square
% error of Euclidean distances. The columns of the matrix X contain the point
% coordinates, n specifies the desired number of clusters.

N = size(X,2); % no. of points
k = size(X,1); % dimension of space

% Start with two clusters obtained by principal axis separation
nc = 1;     % Current number of clusters
Ibig = 1:N; % Initial single cluser encompassing all points
nbig = 1;   % Index of largest cluster
C = sum(X')'/N; % center of gravity
idx = ones(1,N);

while (nc < n)
  % Split largest cluster into two using the principal axis separation
  % algorithm
  [i1,i2] = princaxissep(X(:,Ibig));
  i1 = Ibig(i1); i2 = Ibig(i2);
  n1 = length(i1); n2 = length(i2);
  c1 = sum(X(:,i1)')'/n1; c2 = sum(X(:,i2)')'/n2; 
  C(:,nbig) = c1; C = [C,c2];
  nc = nc+1;
  % Improve clusters by Lloyd-Max iteration
  [C,idx,cds] = lloydmax(X,C);
  % Identify cluster with biggest contribution to mean square error
  [cdm,nbig] = max(cds);
  Ibig = find(idx == nbig);
end