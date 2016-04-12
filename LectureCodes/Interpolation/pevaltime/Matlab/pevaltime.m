% 05.11.2009            pevaltime.m
% comparison of computational time for the evaluation of the interpolating
% polynomial: 
% _Aitken-Neville scheme
% _Matlab polyfit + polyval
% _barycentric formula
% _explicit Legendre polynomials
%
% interpolation with polynomials of degrees 3:200,
% evaluations in 100 points

time=zeros(1,4);
% function to interpolate:
f=@(x) sqrt(x);

%set 'warning off' because polyfit complains
warning off;

% do the same many times and choose the best result:
for k=1:100     %100
  res = [];
  
  % n = increasing polynomial degree
  for n=3:1:200
        fprintf('Degree = %d\n',n);
    t = (1:n);
    y = f(t);
    %x=linspace(1,n,100);                 %vector version
    x=n*rand;
    
    %evaluation in x
    tic;   v1 = ANipoleval(t,y,x);        time(1) = toc;
    tic;   v2 = ipoleval(t,y,x);          time(2) = toc;
    tic;   v3 = intpolyval(t,y,x);        time(3) = toc;
    tic;   v4 = intpolyval_lag(t,y,x);    time(4) = toc;
    res = [res; n,time];
  end
  if (k == 1), finres = res;
  else, finres = min(finres,res); end
end
figure
semilogy(finres(:,1),finres(:,2),'r-',finres(:,1),finres(:,3),'m-',finres(:,1),finres(:,4),'k-',finres(:,1),finres(:,5),'b-');
xlabel('Polynomial degree','Fontsize',14);
ylabel('Computational time [s]','Fontsize',14);
legend('Aitken-Neville scheme','MATLAB polyfit', 'Barycentric formula', 'Lagrange polynomials',2);
warning on

print -depsc2 '../PICTURES/pevaltime.eps';