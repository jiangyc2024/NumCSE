function c = gaussfit (m, d)
n = length(d);
RHS = fft(d);
RHS = RHS(1:m);
RHS = real(RHS);
c(2:m) = 2*RHS(2:m)/n;
c(1) = RHS(1)/n;

