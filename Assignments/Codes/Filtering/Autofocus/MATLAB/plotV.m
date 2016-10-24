function plotV ()
N = 100;
fs = linspace(0,5,N);
V = zeros(1,length(fs));
for i = 1:N
    Bh = fft2(setFocus(fs(i)));
    V(i) = computeV(Bh);
end
plot(fs,V,'r-o');
title('{\bf The 2nd moment of the FFT transformed images}');
xlabel('{\bf f}');
ylabel('{\bf V(f)}');
print -depsc '../PICTURES/plotV.eps'
