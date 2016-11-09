% Generate artifical ``image''
P = kron(magic(3),ones(30,40))*31;
col = [0:1/254:1]'*[1,1,1];
figure; image(P); colormap(col); title('Original');
print -depsc2 '../PICTURES/dborigimage.eps';
% Generate point spread function
L = 5; S = psf(L); 
% Blur image
C = blur(P,S);
figure; image(floor(C));  colormap(col); title('Blurred image');
print -depsc2 '../PICTURES/dbblurredimage.eps';
% Deblur image
D = deblur(C,S);
figure; image(floor(real(D)));  colormap(col); 
fprintf('Difference of images (Frobenius norm): %f\n',norm(P-D));

