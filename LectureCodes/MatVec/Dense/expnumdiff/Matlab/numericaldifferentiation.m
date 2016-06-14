function experr = numericaldifferentiation()
% Numerical differentiation of exponential function with extended precision arithmetic
% Uses the Advanpix extended precision library for MATLAB, wwww.advanpix.com

experr = []; l = 1;
for ndigits=[5 10 15 20 25 30 35]
  mp.Digits(ndigits), % Set number of digits
  x = mp('1.1'); % Evaluation point in extended precision 
  h = mp(2).^[-1:-5:-61]; % width of difference quotient in extended precision
  experr = [experr; abs(((exp(x+h)-exp(x))./h) - exp(x))]; % compute (absolute) error
  leg{l} = sprintf('%2.0d digits',ndigits);
  l = l+1;
end
  
% Graphical output of error
figure('name','numdiff');
loglog(h,experr,'+-');
title('One-sided difference quotient approximation of derivative of e^x');
xlabel('h','fontsize',14);
ylabel('error','fontsize',14);
legend(leg,'Location','best');
print -depsc2 'expnumdiffmultiprecision.eps';
