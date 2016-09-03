% Read image and build matrices, see \cref{mc:imgsegmat} and \eqref{eq:segmat}
P = imread('root2.bmp'); [m,n] = size(P); [A,D] = imgsegmat(P);
% Build scaling matrics 
N = size(A,1); dv = sqrt(spdiags(A,0));
Dm = spdiags(1./dv,[0],N,N); % \Blue{$\VD^{-\nicefrac{1}{2}}$}
Dp = spdiags(dv,[0],N,N);    % \Blue{$\VD^{-\nicefrac{1}{2}}$}
% Build (densely populated !) matrix \Blue{$\wh{\VA}$}
c = Dp*ones(N,1); Ah = Dm*A*Dm + 2*c*c';
% Compute and sort eigenvalues; grossly inefficient \Red{!}
[W,E] = eig(full(Ah)); [ev,idx] = sort(diag(E)); W(:,idx) = W;
% Obtain eigenvector \Blue{$\Vx^{\ast}$} belonging to 2nd smallest generalized
% eigenvalue of \Blue{$\VA$} and \Blue{$\VD$}
v = W(:,1); % v/norm(v);
x = Dm*v;
% Extract segmented image
xs = reshape(x,m,n); Xidx = find(xs>(sum(sum(xs))/(n*m)));