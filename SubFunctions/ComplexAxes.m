function [ax,h] = ComplexAxes(mag,magOff,invertColor)
%% ComplexAxes:
%   INPUTS:
%       mag             :   magntidues rings to display (default = 1:10)
%       invertColor     :   invert background & axes colors if true
%   OUTPUTS:
%       ax              :   axis handle
%       h               :   graphics handle structure
%---------------------------------------------------------------------------------------------------------------------------------
if nargin<=2
   invertColor = false; % default
   if nargin<=1
       magOff = 160; % default
        if ~nargin
            mag = 1:10;
        end
   end
end

if invertColor
    backColor = 'k';
    axColor = 'w';
else
    backColor = 'w';
    axColor = 'k';
end

maxMag = max(mag);

ax = gca; 
hold on
axis square
box on
ax.Color = backColor;
% ax.Colormap = cmap;
ax.FontSize = 12;
% ax.Title.String = 'Complex Gain';
ax.YLabel.String = 'Img';
ax.YLabel.FontSize = 14;
ax.XLabel.String = 'Real';
ax.XLabel.FontSize = ax.YLabel.FontSize;
ax.XColor = axColor;
ax.YColor = axColor;

ang = 90*[0 1 2 3];
for kk = 1:length(ang)
    xx = [0, 2*maxMag*cosd(ang(kk))];
    yy = [0, 2*maxMag*sind(ang(kk))];
    h.origin(kk) = plot(xx,yy,'--','Color',axColor,'LineWidth',1);
end

magR = 1.0*mag;
magX = magR*cosd(magOff);
magY = magR*sind(magOff);
for kk = 1:length(mag)
    h.circle(kk) = PlotCircle(0,0,mag(kk));
    h.circle(kk).Color = [0.5 0.5 0.5 1];
    h.circle(kk).LineStyle = '-';
    h.circle(kk).LineWidth = 1;
    h.text(kk) = text(magX(kk),magY(kk),num2str(mag(kk)),'Color',axColor);
    h.text(kk).FontWeight = 'bold';
%     if mag(kk)==1
%         h.circle(kk).Color = [0.8 0 0 0.5];
%     end
end
ax.XLim = (1.1*maxMag)*[-1 1];
ax.YLim = (1.1*maxMag)*[-1 1];
ax.XTick = ax.YTick;
