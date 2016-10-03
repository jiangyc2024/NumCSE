clear all; close all; tic;
nruns = 5;                   % we average over a few runs
nn = 2.^(4:12);              % dimensions used
ttS  = zeros(1,length(nn));  % times for efficient solution
ttSp = zeros(1,length(nn));  % times for sparse matrix solution
ttF  = zeros(1,length(nn));  % times for full matrix solution
ttM  = zeros(1,length(nn));  % times for Matlab solution

for i = 1:length(nn)
    n = nn(i);
    a = rand(n,1);
    b = rand(n,1);
    
    % 1: smart implementation
    tic;
    for run = 1:nruns
        xS = (b-[0;b(1:n-1)]) ./ a;
    end;    ttS(i) = toc;
    
    % 2: decomposition A=E*D with sparse matrix
    tic;
    for run = 1:nruns
        EE = speye(n) - spdiags(ones(n,1),-1,n,n);
        xSp = a.^(-1) .* (EE*b);
    end;    ttSp(i) = toc;
    
    % 3: decomposition A=E*D with full matrix
    tic;
    for run = 1:nruns
        E = eye(n) - diag (ones(1,n-1),-1);
        xF = a.^(-1) .* (E*b);
    end;    ttF(i) = toc;
 
    % 4: Matlab backslash
    tic;
    A = tril(meshgrid(a,a));
    for run = 1:nruns
        xM = A\b;
    end;    ttM(i) = toc;
    
    % display dimension n and errors
    n_and_errors = [n, norm(xS-xM), norm(xSp-xM), norm(xF-xM)]    
end

figure('name','Structured LSE timings');
loglog(nn, ttM,'m+',  nn, ttF,'ro',  nn, ttSp,'g*',  nn, ttS,'bd',... 
    nn, nn * min(ttS)/10, 'k',    nn, nn.^2 * min(ttS)/100, ...
    'k--', 'linewidth', 2);
xlabel('{\bf dimension n}','fontsize',14);
ylabel('{\bf average runtime (s)}','fontsize',14);
title(sprintf('tic-toc timing averaged over %d runs', nruns),...
    'fontsize',14);
legend('Matlab backslash', 'full matrix', 'sparse matrix',...
    'smart implementation', 'O(n)', 'O(n^2)','location','northwest');
grid on;
print -depsc2 'structeredLSE_timings.eps';
total_time = toc
