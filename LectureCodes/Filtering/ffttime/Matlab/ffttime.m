res = [];
for n=1:1:3000,  y = rand(n,1); c = zeros(n,1);
  t1 = realmax; for k=1:5, tic; 
    omega = exp(-2*pi*i/n); c(1) = sum(y); s = omega;
    for j=2:n, c(j) = y(n);
      for k=n-1:-1:1, c(j) = c(j)*s+y(k);  end
      s = s*omega;
    end
    t1 = min(t1,toc);
  end 
  [I,J] = meshgrid(0:n-1,0:n-1); F = exp(-2*pi*i*I.*J/n);
  t2 = realmax; for k=1:5, tic; c = F*y; t2 = min(t2,toc); end
  t3 = realmax; for k=1:5, tic; d = fft(y); t3 = min(t3,toc); end
  res = [res; n t1 t2 t3];
end

figure('name','FFT timing');
semilogy(res(:,1),res(:,2),'b-',res(:,1),res(:,3),'k-', res(:,1),res(:,4),'r-');
ylabel('{\bf run time [s]}','Fontsize',14);
xlabel('{\bf vector length n}','Fontsize',14);
legend('loop based computation','direct matrix multiplication','MATLAB fft() function',1);
print -deps2c '../PICTURES/ffttime.eps'
