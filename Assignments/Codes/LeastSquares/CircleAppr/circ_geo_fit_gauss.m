function [m,r,GNcorr] = circ_geo_fit_gauss(x,y,m,r)
% geometric L.S. fitting for circles with Gauss-Newton
% inputs:  x,y     column vectors with point coordinates
%          m,r     initial guess for center (m,2x1) and radius (r)
% outputs: m,r     fitted center (m,2x1) and radius (r,1x1)
%          GNcorr  L^infinity norm of GN correction (update)

z = [m(:); r];
running = 1;
GNcorr = [];
while running
    % compute R, function F and Jacobian
    R = sqrt((x-z(1)).^2 + (y-z(2)).^2);
    Fz = R - z(3);
    DFz = [-(x-z(1))./R, -(y-z(2))./R, -ones(size(x))];
    % solve linear least squares:
    s = DFz \ Fz;
    %[Q,R] = qr(DFz, 0);   s = R\(Q'*Fz);
    z = z-s;          % compute GN correction and store the norm
    GNcorr = [GNcorr; norm(s,inf)];
    running = (GNcorr(end) > 1e-12);
end
m = z(1:2);
r = z(3);

% compute & print norm of the correction, rate, order of convergence:
rate = GNcorr(2:end)./GNcorr(1:end-1);
lGN = log(GNcorr);
order = (lGN(3:end)-lGN(2:end-1))./(lGN(2:end-1)-lGN(1:end-2));
GN_corr_rate_order = [GNcorr(:),[0;rate(:)],[0;0;order(:)]]