% MATLAB test script: assessing the gain from using linsolve

K = 3; % number of runs (to offset OS activity)
T = []; % Matrix for recording results
for n=2.^(4:12)
 % Create test matrix   
 A = triu(diag(1:n) + ones(n,1)*ones(1,n));
 % Slight perturbation below the diagonal
 A = A + flipud(diag(eps*rand(n,1)));  % \Label{lisot:1}
 b = rand(n,1); 
 t_b = realmax; t_l = realmax;
 opts.UT = true;
 for k=1:K
    tic; xb = A\b; t = toc; t_b = min(t_b,t);
    tic; xl = linsolve(A,b,opts); t = toc; t_l = min(t_l,t);
 end
 T = [T; n t_b t_l norm(xb-xl)];
end

figure; loglog(T(:,1),T(:,2),'r+',T(:,1),T(:,3),'m*');
xlabel('{\bf matrix size n}','fontsize',14); 
ylabel('{\bf runtime for direct solver [s]}','fontsize',14);
legend('backslash','linsolve(UT)','location','best');
print -depsc2 '../PICTURES/linsolvetest.eps';