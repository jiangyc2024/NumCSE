function y = Taylor_step(odefun, t, y0, h, varargin)

y1 = odefun(t, y0, 'f' ,  [], varargin{:});
y2 = odefun(t, y0, 'df',  y1, varargin{:});
y3 = odefun(t, y0, 'df',  y2, varargin{:}) + ...
     odefun(t, y0, 'd2f', y1, varargin{:});
 
y  = y0 + h*y1 + h^2*y2/2 + h^3*y3/6;