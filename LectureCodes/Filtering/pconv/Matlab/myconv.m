function y = myconv(h,x)
n = length(h); 
% Zero padding, cf. \eqref{eq:zeropad}
h = [h;zeros(n-1,1)]; 
x = [x;zeros(n-1,1)]; 
% Periodic discrete convolution of length \Blue{$2n-1$}
y = pconvfft(h,x);
