function [low,high] = freqfilter(y,k)
m = length(y)/2; c = fft(y); 
clow = c; clow(m+1-k:m+1+k) = 0; 
chigh = c-clow;
low = real(ifft(clow)); 
high = real(ifft(chigh));
