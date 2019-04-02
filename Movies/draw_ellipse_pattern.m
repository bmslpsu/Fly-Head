function [h,top,bot] = draw_ellipse_pattern(center,L,ratio,ecc,spatFreq,ang)
%% draw_circle_pattern: draw a circle on figure window
%   INPUT:
%       center      : center point [x,y]
%       L           : length of ellipse
%       ratio       : portion of ellipse length below center
%       ecc       	: eccentricity
%       ang       	: angle of ellispe in degrees
%   OUTPUT:
%       h           : handles for ellipse figure data
%       top         : top point of ellipse
%       bottom    	: bottom point of ellipse
%---------------------------------------------------------------------------------------------------------------------------------
% Example Input %
% center = [0,0];
% L = 5;
% ratio = 0.5;
% ecc = 0;
% spatFreq = 7.5;
% ang = 45;
%---------------------------------------------------------------------------------------------------------------------------------
% Top point
top(1) = center(1) + L*(1-ratio)*sind(ang);
top(2) = center(2) + L*(1-ratio)*cosd(ang);

% Bottom point
bot(1) = center(1) - L*ratio*sind(ang);
bot(2) = center(2) - L*ratio*cosd(ang);

% All [x,y] points defining ellipse
xx = [center(1) top(1) bot(1)];
yy = [center(2) top(2) bot(2)];

hold on
h = {};
% h{end+1} = plot(xx,yy,'-k','LineWidth',2); % vector
% h{end+1} = plot(center(1),center(2),'-ko','MarkerSize',5,'MarkerFaceColor','k'); % center point
% h{end+1} = plot(top(1),top(2),'-ko','MarkerSize',5,'MarkerFaceColor','k'); % top point
% h{end+1} = plot(bot(1),bot(2),'-ko','MarkerSize',5,'MarkerFaceColor','k'); % bottom point

% Calculate border points
a = 1/2*sqrt((top(1)-bot(1))^2+(top(2)-bot(2))^2);
b = a*sqrt(1-ecc^2);
w = atan2(bot(2)-top(2),bot(1)-top(1));

t = linspace(0,2*pi);
X = a*cos(t);
Y = b*sin(t);
x = (top(1)+bot(1))/2 + X*cos(w) - Y*sin(w);
y = (top(2)+bot(2))/2 + X*sin(w) + Y*cos(w);
h{end+1} = plot(x,y,'-k','LineWidth',2); % border

tt = linspace(0,2*pi,360/(spatFreq));
X = a*cos(tt);
Y = b*sin(tt);
x = (top(1)+bot(1))/2 + X*cos(w) - Y*sin(w);
y = (top(2)+bot(2))/2 + X*sin(w) + Y*cos(w);
h{end+1} = plot(x,y,'sg','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','g'); % border

hold off
end