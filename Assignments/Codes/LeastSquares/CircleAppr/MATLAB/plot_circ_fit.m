function plot_circ_fit ()
close all; tic;
% data:
x = [0.7; 3.3; 5.6; 7.5; 6.4; 4.4; 0.3; -1.1];
y = [4; 4.7; 4; 1.3; -1.1; -3; -2.5; 1.3];
% test data with known center and radius + perturbation:
%  mtest = [2,1];  rtest = 2;  angs = 2*pi*rand(10,1);
%  x = mtest(1) + rtest*cos(angs)+0.3*randn(size(angs));
%  y = mtest(2) + rtest*sin(angs)+0.3*randn(size(angs));

% initial guess
r0 = 1;
m0 = [mean(x), mean(y)];

figure; subplot(1,2,1);
plot(x,y,'m+', 'linewidth',2);
axis equal; hold on;

[m2, r2, err_GN] = circ_geo_fit_gauss(x, y, m0, r0);
draw_circle(m2, r2, 'g--');

[m3, r3, err_DN] = circ_geo_fit_newton(x, y, m0, r0);
draw_circle(m3, r3, 'r-.');

% print the centers and the radii
centers_and_radii = [m2,m3;r2,r3]

% plot the centers of the circles:
plot(m2(1), m2(2), 'g+', m3(1), m3(2), 'r+');
legend('data', 'Gauss-Newton', 'Newton', 'location', 'so');

% plot convergence:
subplot(1,2,2);
semilogy(1:length(err_GN), err_GN, 'bo-',...
    1:length(err_DN), err_DN, 'ro-', 'linewidth',2)
xlabel('iteration');  ylabel('error');
legend('Gauss-Newton','Newton');
print -depsc2 '../PICTURES/plot_circ_fit.eps';
end

function draw_circle(C,r,s)
% plot the circle given by the center C = [c1 c2]
% and the radius r in the color s.
theta = (0:0.02:2)*pi;
plot(C(1)+r*cos(theta), C(2)+r*sin(theta), s, 'linewidth',2);
end
