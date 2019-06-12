function [FIG] = MakeFig_SOS_HeadFree_pat2head_ComplexGain()
%% MakeFig_SOS_HeadFree_pat2head_ComplexGain:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');
%%
figNum = 1;
catIdx = 5;
xIdx = 1;

filename = 'SOS_HeadFree_pat2head_ComplexGain';

Freq = HeadFree.TRIAL{1}{1,catIdx}.IOFreq;
nFreq = length(Freq);
fLabel = cellfun(@(x) [num2str(x) ' Hz'], num2cell(Freq),'UniformOutput',false);

pp = 1;
for kk = 1:HeadFree.N{1,1}
    for ii = 1:size(HeadFree.TRIAL{kk},1)
        CmplxGain(pp,:) = HeadFree.TRIAL{kk}{ii,catIdx}.IOCmplxGain(:,xIdx);
        pp = pp +1;
    end
end

Real = real(CmplxGain);
Imag = imag(CmplxGain);
CmplxGain = Real + 1i*Imag;

cnd = (Real<0) & (Imag>0);
Real(cnd) = -Real(cnd);
Imag(cnd) = -Imag(cnd);

REAL        = mean(Real,1);
IMAG        = mean(Imag,1);
CMPLXGAIN   = REAL + 1i*IMAG;

GAIN    = 2*abs(CMPLXGAIN);
PHASE   = rad2deg(angle((CMPLXGAIN)));

cList= jet(nFreq);

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 800];
FIG.Name = filename;
movegui(FIG,'center')
hold on

gains = 0.2:0.2:1;
maxGain = max(gains);

ax = gca; hold on
axis square
box on
ax.Color = 'w';
ax.FontSize = 12;
ax.Title.String = 'Complex Gain';
ax.YLabel.String = 'Img';
ax.YLabel.FontSize = 14;
ax.XLabel.String = 'Real';
ax.XLabel.FontSize = ax.YLabel.FontSize;

axVector = 2*(-maxGain:0.1:maxGain);
axColor = [1 0 0];
h.xorigin = plot(axVector,0*axVector,'--','Color',axColor,'LineWidth',1);
h.yorigin = plot(0*axVector,axVector,'--','Color',axColor,'LineWidth',1);
gainOff = 160;
gainR = 1.0*gains;
gainX = gainR*cosd(gainOff);
gainY = gainR*sind(gainOff);
for kk = 1:length(gains)
    h.circle = PlotCircle(0,0,gains(kk));
    h.circle.Color = [0.5 0.5 0.5 1];
    h.circle.LineStyle = '-';
    h.circle.LineWidth = 1;
    h.text = text(gainX(kk),gainY(kk),num2str(gains(kk)));
    h.text.FontWeight = 'bold';
    if gains(kk)==1
        h.circle.Color = [0.8 0 0 0.5];
    end
end
ax.XLim = (1.1*maxGain)*[-1 1];
ax.YLim = (1.1*maxGain)*[-1 1];
ax.XTick = ax.YTick;

for jj = 1:nFreq
    h.trial = scatter(Real(:,jj), Imag(:,jj), 40, 'o','MarkerEdgeColor','k',...
      'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.65, 'LineWidth', 0.5);
end

for jj = 1:nFreq
  	h.rr = plot([0 REAL(jj)],[0 IMAG(jj)],'Color',[0 0 0 1],'LineWidth',1);
    
    rSTD = PolarSTD(Real(:,jj),Imag(:,jj),[REAL(jj) IMAG(jj)]);
    
	[h.std] = draw_ellipse([REAL(jj) IMAG(jj)], 3*rSTD, 0.5, 0, 90, cList(jj,:)); hold on
    h.std{1}.FaceAlpha = 0.4;
    for kk = 3:length(h.std)
       delete(h.std{kk})
    end
    
    h.grand = scatter(REAL(jj),IMAG(jj),1,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
        'MarkerFaceAlpha',1,'LineWidth',1.5);
    
    h.leg = scatter(REAL(jj),IMAG(jj),40,'o','MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerFaceAlpha',1,'LineWidth',1.5);

    hh(jj) = h.grand;
end

leg = legend(hh,fLabel);
leg.Title.String = 'Frequency';
leg.Location = 'northwest';
legend boxoff







end