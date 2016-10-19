t = 0:63; x = sin(2*pi*t/64)+sin(7*2*pi*t/64); 
y = x + randn(size(t)); %distortion
