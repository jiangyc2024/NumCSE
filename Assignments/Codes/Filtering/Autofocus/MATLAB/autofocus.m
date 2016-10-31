function f0 = autofocus ()
iters = 6;
MAXF = 5;
df = MAXF/1e2;
f0 = MAXF/2;
step = MAXF/2;
for i=1:iters
    % numerical differentiation
    Vl = computeV(fft2(setFocus(f0)));
    Vr = computeV(fft2(setFocus(f0+df)));
    dV = Vr - Vl;
    % bisection method
    step = step/2;
    f0 = f0 + sign(dV) * step;
end

