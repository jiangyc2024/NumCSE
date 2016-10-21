function pts = Xmastree(fn)
% MATLAB function creating the points of a tree
pts = [1 0; 4 -3 ; 1 -3 ;1 -4 ; -1 -4; -1 -3; -4 -3; ...
       -1 0; -3 0; -1 2; -2 2; 0 4; 2 2; 1 2; 3 0]';
% Close polygon
ppts = [pts,pts(:,1)];
% Draw tree
figure('name','tree');
plot(ppts(1,:),ppts(2,:),'b-',ppts(1,:),ppts(2,:),'rp'); 
grid on; axis([-5 5 -5 5]); axis equal;

% Print, if filename was specified
if (nargin > 0),  print('-depsc2',[fn , '.eps']); end

