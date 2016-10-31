% Tracking of periodicity in data
% Data obtained from \textsf{Google Trends}, keyword ``Vorlesungsverzeichnis''
% Exported as .csv-file, non data lines removed, preprocessed by command \texttt{cut -f 2 -d,}

% read ASCI data from file
y = dlmread('trend.dat'); n = length(y);

figure('name','data');
plot(y,'r-'); grid on;
title('{\bf Google searches for keyword Vorlesungsverzeichnis}','fontsize',14);
xlabel('{\bf week (1.1.2004-31.12.2010)sa}','fontsize',14);
ylabel('{\bf relative no. of searches}','fontsize',14);
print -depsc2 '../PICTURES/searchdata.eps';

% Periodicity analysis by means of DFT
c = fft(y); 
p = abs(c(2:floor((n+1)/2))).^2; % Power spectrum
figure('name','Fourier spectrum');
plot(2:floor((n+1)/2),p,'m-'); grid on; hold on;
[mx,idx] = sort(p,'descend');
plot(1+idx(1:4),p(idx(1:4)),'rp');
xlabel('{\bf index j of Fourier component}','fontsize',14);
ylabel('{\bf |c_j|^2}','fontsize',14);
title('{\bf Energy spectrum of web search data}','fontsize',14);
print -depsc2 '../PICTURES/fourierdata.eps';
