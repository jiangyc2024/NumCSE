Rank1Circuit_A;         %generate Q,R,b
u = [zeros(14,1); 1; -1];
v = zeros(16,1);
max_res = 0;

Rx = 49.9:-0.1:1;

I = [I_15_16, zeros(1, length(Rx))];
for n = 1:length(Rx)
    Delta  = 1/Rx(n) - 1/50;
    v(15)   =  Delta;
    v(16)   = -Delta;
    [Q1,R1] = qrupdate(Q,R,u,v);
    U_all   = R1\(Q1'*b);
    I(n+1)  = (U_all(15)-U_all(16))/Rx(n);
% test your code with the relative residual ||A_mod*U-b|| / ||b||:
    max_res = max(norm((A+u*v')*U_all-b)/norm(b),max_res);
    
    % cheap alternative version:
    % U_red = R1(15:16,15:16)\(Q1(:,15:16)'* b);
    % I_red(n+1) = (U_red(1) - U_red(2))/Rx(n);
end
sprintf('Maximal relative residual = %g\n', max_res)
% error of the cheap version:
% norm(I(2:end)-I_red(2:end))

close all;
plot([50 Rx],real(I),'r', [50 Rx],imag(I),'b','Linewidth',2);
legend('Real(I_{15,16})', 'Imag(I_{15,16})','location','southeast');
xlabel('Resistance Rx','fontsize',14);
ylabel('Current between nodes 15 and 16','fontsize',14);
print -depsc2 'ex_Rank1Circuit_plotI.eps'