
/*
void sinft2d(Y)
{
	[m,n] = size(Y);
	C = fft([zeros(1,n); Y;...
			 zeros(1,n);...
			 -Y(end:-1:1,:)]);
	C = i*C(2:m+1,:)'/2;
	C = fft([zeros(1,m); C;...
			 zeros(1,m);...
			 -C(end:-1:1,:)]);
	C= i*C(2:n+1,:)'/2;

	return C;
}*/
