function s = sinetrans(y)
	n = length(y)+1;
	sinevals = imag(exp(i*pi/n).^(1:n-1))';
	yt = [0; (sinevals.*(y+y(end:-1:1)) + 0.5*(y-y(end:-1:1)))];
	c = fftreal(yt);
	s(1) = dot(sinevals,y);
	for k=2:n-1
		if (mod(k,2) == 0), s(k) = -imag(c(k/2+1));
		else, s(k) = s(k-2) + real(c((k-1)/2+1)); 
	end
end
