function trussplot(pos,top,ptstyle,rodstyle)
% Draws a truss structure with mass point positions passed in the $n\times 2$-matrix
% \texttt{poss} and the connectivity encoded in the sparse matrix \texttt{top}
if (nargin < 4), rodstyle = 'b-'; end
if (nargin < 3), ptstyle = 'ro'; end
plot(pos(:,1),pos(:,2),ptstyle,'markersize',5); 
ax = axis; l = [ax(2)-ax(1); ax(4)-ax(3)];
axis([ax(1)-0.1*l(1),ax(2)+0.1*l(1),ax(3)-0.1*l(2),ax(4)+0.1*l(2)]);
hold on;
[Iidx,Jidx] = find(top); idx = [Iidx,Jidx];
for ij=idx'
  plot([pos(ij(1),1) pos(ij(2),1)],[pos(ij(1),2) pos(ij(2),2)],rodstyle);
end
