c = fft(y); p = (c.*conj(c))/64;

figure('Name','power spectrum');
bar(0:31,p(1:32),'r');
set(gca,'fontsize',14);
axis([-1 32 0 max(p)+1]);
xlabel('{\bf index k of Fourier coefficient}','Fontsize',14);
ylabel('{\bf |c_k|^2}','Fontsize',14);
