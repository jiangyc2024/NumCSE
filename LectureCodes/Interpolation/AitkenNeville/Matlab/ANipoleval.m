% 16.06.2009            ANipoleval.m
% Evaluation of the interpolation polynomials with Aitken-Neville scheme
%
% input:    t   nodes
%           y   values in t
%           x   evaluation point (scalar)


 function v = ANipoleval(t,y,x)
% version for scalar y:
for i=1:length(y)
    for k=i-1:-1:1
        y(k) = y(k+1)+(y(k+1)-y(k))*(x-t(i))/(t(i)-t(k));
    end
end
v = y(1);


% % version for vector y:

% yy=ones(length(x), length(y)) * diag(y);        %each line is a copy of y
% for i=1:length(y)
%     for k=i-1:-1:1
%         yy(:,k) = yy(:,k+1)+(yy(:,k+1)-yy(:,k)).*(x(:)-t(i))/(t(i)-t(k));
%     end
% end
% v = yy(:,1);
