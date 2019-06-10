function [h] = PlotCircle(x,y,r)
%% PlotCircle:
%   INPUTS:
%       x   :   x-position center
%       y   :   y-position center
%       r   :   radius
%   OUTPUTS:
%       h   :   plot handle
%---------------------------------------------------------------------------------------------------------------------------------
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
end