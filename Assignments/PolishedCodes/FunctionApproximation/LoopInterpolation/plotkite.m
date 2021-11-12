t = linspace(0,1,100);
cx = cos(2*pi*t) + 2/3 * cos(4*pi*t);
cy = 3/2 * sin(2*pi*t);

figure;
plot(cx,cy);
hold on;
title("Kite-shaped curve");

grid on;

