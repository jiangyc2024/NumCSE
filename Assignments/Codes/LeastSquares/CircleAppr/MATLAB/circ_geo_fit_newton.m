function [m,r,Ncorr] = circ_geo_fit_newton(x,y,m,r)
% geometric L.S. fitting for circles with Newton method
% inputs:  x,y    column vectors with point coordinates
%          m,r    initial guess for center (m,2x1) and radius (r)
% outputs: m,r    fitted center (m,2x1) and radius (r,1x1)
%          Ncorr  L^infinity norm of Newt corrections (update)

z = [m(:); r];   % convert circles params in a 3-vector
running = 1;
Ncorr = [];
while running
    % gradient and Hessian of Phi, through auxiliary function:
    [GradPhi, HessPhi] = GradHess(x,y,z);
    % compute Newton correction and store its norm:
    s = HessPhi\GradPhi;
    z = z-s;       
    Ncorr = [Ncorr; norm(s,inf)];
    running = (Ncorr(end) > 1e-12);
end
m = z(1:2);
r = z(3);

% compute & print norm of the correction, rate, order of convergence:
rateN = Ncorr(2:end)./Ncorr(1:end-1);
lN = log(Ncorr);
orderN = (lN(3:end)-lN(2:end-1))./(lN(2:end-1)-lN(1:end-2));
Newt_corr_rate_order = [Ncorr(:),[0;rateN(:)],[0;0;orderN(:)]]



function [GradPhi, HessPhi] = GradHess(x,y,z)
% inputs:    N-vectors x,y and a 3-vector z
% outputs:   3-vector GradPhi,   3x3 matrix HessPhi
N  = length(x);
xm = x - z(1);
ym = y - z(2);
R  = sqrt(xm.^2 + ym.^2);
% gradPhi is a column vector of length 3
GradPhi = [dot(z(3)./R-1, xm);  dot(z(3)./R-1, ym);   N*z(3)-sum(R)];

%MixedSum2 = sum(xm.*ym./R.^2);
SumX = sum(xm./R);
SumY = sum(ym./R);
MixedSum = z(3)* sum(xm.*ym./R.^3);
InvSum   = z(3)*sum(1./R);

HessPhi = [N-InvSum+z(3)*sum(xm.^2 ./ R.^3),  MixedSum, SumX;
    MixedSum,  N - InvSum + z(3)*sum(ym.^2 ./ R.^3), SumY;
    SumX,  SumY,  N]; 