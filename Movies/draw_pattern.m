function [h] = draw_pattern(xpos,ypos,pattern,center,radius,width)
%% draw_pattern: draw pattern ring at given channel positions
%   INPUT:
%       xpos       	: pattern x-channel position (1 to pattern.x_num)
%       ypos       	: pattern y-channel position (1 to pattern.y_num)
%       patttern    : pattern structure (requires at least .Pats & .gs_val)
%       center      : center point [x,y]
%       radius      : radius for display
%       width       : width for display
%   OUTPUT:
%       h           : handles for figure
%---------------------------------------------------------------------------------------------------------------------------------
rin  = radius;
rout = radius + width;
x1 = center(1);
y1 = center(2);
sA = 3.75*pi/180; % angle pixel subtends

pat = pattern.Pats(1,:,xpos,ypos); % top row of pattern
patS = circshift(pat,[0 0]); % shift pattern to fly reference frame

I = find(patS~=0);
theta = (I.*3.75)*(2*pi/360); % lit pixels
theta_ALL = deg2rad(3.75*(1:96));

% Make pattern ring
hold on
h = cell(length(theta_ALL),1);
for jj = 1:length(theta_ALL)
    xin     = x1 + rin*cos(theta_ALL(jj));
    xout    = x1 + rout*cos(theta_ALL(jj));
    xinN    = x1 + rin*cos(theta_ALL(jj) + sA);
    xoutN   = x1 + rout*cos(theta_ALL(jj) + sA);
    yin     = y1 + rin*sin(theta_ALL(jj));
    yout    = y1 + rout*sin(theta_ALL(jj));
    yinN    = y1 + rin*sin(theta_ALL(jj) + sA);
    youtN   = y1 + rout*sin(theta_ALL(jj) + sA);

    if sum(ismember(theta, theta_ALL(jj))) == 1 % if lit
        h{jj} = patch([xout, xoutN, xinN, xin], [yout, youtN,yinN, yin],'g','linestyle','none',...
            'FaceAlpha',pat(jj)*(1/(2^(pattern.gs_val)-1)));
    else % if dark
        h{jj} = patch([xout, xoutN, xinN, xin], [yout, youtN,yinN, yin],'k','linestyle','none');
    end
end
hold off
end