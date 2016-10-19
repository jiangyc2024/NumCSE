function z=pconv(u,x)
n = length(x); z = zeros(n,1);
for i=1:n, z(i)=dot(conj(u), x([i:-1:1,n:-1:i+1]));
end
