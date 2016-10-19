function S = psf(L)
[X,Y] = meshgrid(-L:1:L,-L:1:L);
S = 1./(1+X.*X+Y.*Y);
S = S/sum(sum(S));
