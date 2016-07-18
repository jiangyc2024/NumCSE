M = imread('root.pbm');
[m,n] = size(M);
fprintf('%dx%d grey scale pixel image\n',m,n);
figure; image(M); title('ETH view');
col = [0:1/215:1]'*[1,1,1]; colormap(col);
