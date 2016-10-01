function [Q,A] = qrgivens(A)
% in situ QR decomposition of \emph{square} matrix \texttt{A}, same semantics as
% \matlab{} built \texttt{qr()} function
n=size(A,1);
Q=eye(n); % Assemble rotations in matrix, alternative see Rem.~\ref{rem:orthstore}
for i=1:(n-1)
  for j=n:-1:(i+1)
    G=planerot([A(j-1,i);A(j,i)]); % see Code~\ref{planerot}
    A([j-1,j],:)=G*A([j-1,j],:);
    Q(:,[j-1,j])=Q(:,[j-1,j])*G';
  end
end


