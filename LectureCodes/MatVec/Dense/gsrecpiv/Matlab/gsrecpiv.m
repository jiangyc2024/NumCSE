function A = gsrecpiv(A) 
n = size(A,1);
if (n > 1)
  [p,j] = max(abs(A(:,1))./max(abs(A)')'); %\label{gsrp:1}
  if (p < eps*norm(A(:,1:n),1)), disp('A nearly singular'); end %\label{gsrp:2} 
    A([1,j],:) = A([j,1],:);                     % \label{gsrp:3} 
    fac = A(2:end,1)/A(1,1);                     % \label{gsrp:f}
    C = gsrecpiv(A(2:end,2:end)-fac*A(1,2:end)); % \label{gsrp:4}  
    A = [A(1,:) ; -fac, C ];                     % \label{gsrp:5}  
end
