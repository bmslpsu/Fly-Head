function [ax,h] = ComplexAxes(mag)
%% ComplexAxes:
%   INPUTS:
%       mag     :   magntidues rings to display (default = 1:10)
%   OUTPUTS:
%       ax      :   axis handle
%       h       :   graphics handle structure
%---------------------------------------------------------------------------------------------------------------------------------
if ~nargin
    mag = 1:10;
end

maxMag = max(mag);

ax = gca; 
hold on
axis square
box on
ax.Color = 'w';
% ax.Colormap = cmap;
ax.FontSize = 12;
% ax.Title.String = 'Complex Gain';
ax.YLabel.String = 'Img';
ax.YLabel.FontSize = 14;
ax.XLabel.String = 'Real';
ax.XLabel.FontSize = ax.YLabel.FontSize;

axVector = 2*(-maxMag:0.1:maxMag);
axColor = [1 0 0];
h.xorigin = plot(axVector,0*axVector,'--','Color',axColor,'LineWidth',1);
h.yorigin = plot(0*axVector,axVector,'--','Color',axColor,'LineWidth',1);
magOff = 160;
magR = 1.0*mag;
magX = magR*cosd(magOff);
magY = magR*sind(magOff);
for kk = 1:length(mag)
    h.circle = PlotCircle(0,0,mag(kk));
    h.circle.Color = [0.5 0.5 0.5 1];
    h.circle.LineStyle = '-';
    h.circle.LineWidth = 1;
    h.text = text(magX(kk),magY(kk),num2str(mag(kk)));
    h.text.FontWeight = 'bold';
    if mag(kk)==1
        h.circle.Color = [0.8 0 0 0.5];
    end
end
ax.XLim = (1.1*maxMag)*[-1 1];
ax.YLim = (1.1*maxMag)*[-1 1];
ax.XTick = ax.YTick;
