P = imread('root2.bmp'); [m,n] = size(P); [A,D] = imgsegmat(P);
[V,E] = eig(full(A+D),full(D)); % grossly inefficient \Red{!}
xs = reshape(V(:,2),m,n); Xidx = find(xs>(sum(sum(xs))/(n*m)));