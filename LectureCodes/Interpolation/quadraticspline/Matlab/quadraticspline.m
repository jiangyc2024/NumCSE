p = (t(1)-1)*ones(1, length(t)-1);
for j=1:n-1
  if (c(j) ~= c(j+1))
    p(j)=(y(j+1)-y(j)+t(j)*c(j)-t(j+1)*c(j+1))/(c(j)-c(j+1));
  end
  if ((p(j)<t(j))|(p(j)>t(j+1)))
    p(j) = 0.5*(t(j)+t(j+1));
  end;
end
