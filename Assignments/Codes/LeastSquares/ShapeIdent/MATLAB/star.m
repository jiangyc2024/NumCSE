function pts = star(fn)
% MATLAB function creating the points of a star
% points of rays
outpts = [cos((0:4)*2*pi/5);sin((0:4)*2*pi/5)];
% inner corners
inpts = 0.5*[cos((0:4)*2*pi/5+pi/5);sin((0:4)*2*pi/5+pi/5)];
% merge both point sets
pts = zeros(2,10); pts(:,1:2:9) = outpts; pts(:,2:2:10) = inpts;
% centers of rays are stored in mp
inpts = [inpts(:,end),inpts];
mp = (outpts+inpts(:,1:end-1)+inpts(:,2:end))/3;

% Draw star
figure('name','star');
plot([pts(1,:),pts(1,1)],[pts(2,:),pts(2,1)],'b-',...
     [pts(1,:),pts(1,1)],[pts(2,:),pts(2,1)],'rp',...
     inpts(1,:),inpts(2,:),'m--',...
     mp(1,:),mp(2,:),'rp');
grid on; axis equal; axis([-1 1 -1 1]);

% Print, if filename was specified
if (nargin > 0),  print('-depsc2',[fn , '.eps']); end
pts = [pts mp];