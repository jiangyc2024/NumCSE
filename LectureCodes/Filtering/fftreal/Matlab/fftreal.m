function c = fftreal(y)
n = length(y); m = n/2;
if (mod(n,2) ~= 0), error('n must be even'); end
y = y(1:2:n)+i*y(2:2:n); h = fft(y); h = [h;h(1)];
c = 0.5*(h+conj(h(m+1:-1:1))) - ...
    (0.5*i*exp(-2*pi*i/n).^((0:m)')).*...
    (h-conj(h(m+1:-1:1)));
c = [c;conj(c(m:-1:2))];
