clear all; close all; tic;
n = 100;                       % dimension
NumExp = 10;
alphas = 2.^(1:NumExp);        % different parameters
NumSV  = 2;                    % # s.v. to be computed

svs    = zeros(NumExp,NumSV); svsmin = zeros(NumExp,NumSV);
normAx = zeros(NumExp,1);
for j = 1:NumExp
    alpha = alphas(j);
    A = speye(n);   A(1,2:n) = alpha; 
    % compute and store the NumSV largest s.v. of A:
    svs(j,:)  = svds(A, NumSV);
    % compute the NumSV smallest s.v. of A:
    svsmin(j,:) = svds(A, NumSV,0);
    
    % test our vector:
    x = [1/alpha; sqrt((1-1/alpha^2)/(n-1) )*ones(n-1,1)];
    normAx(j) = norm(A*x)/norm(x);
end
figure; 
loglog(alphas,normAx,'k', alphas,svs(:,1),'ro',alphas,svsmin(:,2),'+b',...
    alphas,sqrt((alphas.^2-1)*(n-1)),'g--',...
    alphas,svs(:,2),'ro',alphas,svsmin(:,1),'+b','linewidth',2);
legend('||Ax||/||x||', '2 largest s.vs.',...
    '2 smallest s.vs.','sqrt((a^2-1)(n-1))','location','nw');
xlabel('alpha','fontsize',14);
print -depsc2 'ex_SVDalpha.eps'
toc