function sdirkConv


% initialization
gamma = (3+sqrt(3))/6;
y0 = [1;0];
T = 10;
N = 20*2.^[0:9];

%analytic solution
sol=@(t,y0) 1/3*exp(-t/2) .* ( 3*y0(1) * cos( sqrt(3)*t/2 ) + ...
						sqrt(3)*y0(1) * sin( sqrt(3)*t/2 ) + ...
						2*sqrt(3)*y0(2) * sin( sqrt(3)*t/2 ) );

% allocate memory for results
err = nan(size(N));

% compute error for different step lengths
for j=1:length(N)
    
    % step size
    h = T/N(j);
    
    % set initial value
    y = y0;
    
    % compute solution
    for i=1:N(j)
        y = sdirkStep(y,h,gamma);
    end
       
    %compute error
	err(j) = norm( y(1) - sol(T,y0) );
end

%loglog
figure;
loglog(N,err,'bx:')
grid on;
title(['Error vs. Step Size for SDIRK, \gamma =',...
	 num2str(gamma)]);
xlabel('N');
ylabel('error');

% calculate order of convergence
fit = polyfit(log(N),log(err),1);
fprintf('numerically approximated order of convergence: %f\n',fit(1));
