function [FIG] = MakeFig_SOS_HeadFree_pat2head_ComplexGain_NEW()
%% MakeFig_SOS_HeadFree_pat2head_ComplexGain_NEW:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%

root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');

%%
clearvars -except HeadFree

catIdx = 5;
xIdx = 1;

filename = 'SOS_HeadFree_pat2head_ComplexGain';

Freq = HeadFree.TRIAL{1}{1,catIdx}.IOFreq;
nFreq = length(Freq);
fLabel = string(Freq) + " Hz";

pp = 1;
for kk = 1:HeadFree.N{1,1}
    for ii = 1:size(HeadFree.TRIAL{kk},1)
        CmplxGain(pp,:) = HeadFree.TRIAL{kk}{ii,catIdx}.IOCmplxGain(:,xIdx);
        pp = pp +1;
    end
end

Real        = real(CmplxGain);
Imag        = imag(CmplxGain);

for jj = 1:nFreq
    for kk = 1:size(Real,1)
        if jj~=1 && ( Real(kk,jj)<0 && Imag(kk,jj)>0 )
            Imag(kk,jj) = -Imag(kk,jj);
            Real(kk,jj) = -Real(kk,jj);
        end
    end
end

CmplxGain   = Real + 1i*Imag;

Gain        = abs(CmplxGain);
Phase       = rad2deg(angle((CmplxGain)));

REAL        = mean(Real,1);
IMAG        = mean(Imag,1);

GAIN_Mean 	= mean(Gain,1);
PHASE_Mean 	= rad2deg(circ_mean(deg2rad(Phase),[],1));
GAIN_STD    = std(Gain,[],1);
[~,PHASE_STD] = circ_std(deg2rad(Phase),[],[],1);
PHASE_STD = rad2deg(PHASE_STD);

CMPLXGAIN   = REAL + 1i*IMAG;
GAIN_Mean  	= abs(CMPLXGAIN);
PHASE_Mean 	= rad2deg(angle((CMPLXGAIN)));

cList = hsv(nFreq);

gains = 0.2:0.2:1;

%% Complex Gain
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
FIG.Name = filename;
movegui(FIG,'center')
hold on

[ax,h.label] = ComplexAxes(gains);
ax.FontSize = 8;
ax.XLabel.FontSize = 8;
ax.YLabel.FontSize = 8;

for jj = 1:nFreq
    h.trial = scatter(Real(:,jj), Imag(:,jj), 15, 'o','MarkerEdgeColor','k',...
      'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.65, 'LineWidth', 0.5);
end

for jj = 1:nFreq
  	h.rr = plot([0 REAL(jj)],[0 IMAG(jj)],'Color',[0 0 0 1],'LineWidth',1);
    
    rSTD = PolarSTD(Real(:,jj),Imag(:,jj),[REAL(jj) IMAG(jj)]);
    
	[h.std] = draw_ellipse([REAL(jj) IMAG(jj)], 3*rSTD, 0.5, 0, 90 - PHASE_Mean(jj), cList(jj,:)); hold on
    h.std{1}.FaceAlpha = 0.4;
    for kk = 3:length(h.std)
       delete(h.std{kk})
    end
    h.std{2}.LineWidth = 1;
    
    h.leg = scatter(REAL(jj),IMAG(jj),1,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
        'MarkerFaceAlpha',1,'LineWidth',1.5);
    
    h.grand = scatter(REAL(jj),IMAG(jj),10,'o','MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerFaceAlpha',1,'LineWidth',1.5);

    hh(jj) = h.leg;
end
leg = legend(hh,string(Freq));
leg.Title.String = 'Frequency (Hz)';
leg.Location = 'northwest';
legend boxoff

%% BODE
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
movegui(FIG,'center')
hold on

ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 12;
    ax1.XLim = [0 10]; 
    ax1.YLim = [0 1];
    ax1.XTickLabel = '';
    ax1.XLabel.FontSize = 14;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
  	[h.std,h.gain] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx), HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4}, 1, HeadFree.N{1,1}, 'k', [0.4 0.4 0.6], 0.5, 1);
    h.std.EdgeColor = 'none';
 	h.gain.Marker = '.';
    h.gain.MarkerSize = 20;
    
	[h.std,h.gain] = PlotPatch(GAIN_Mean, GAIN_STD, Freq, 1, HeadFree.N{1,1}, 'g', [0.4 0.4 0.6], 0.5, 1);
    h.std.EdgeColor = 'none';
	h.gain.Marker = '.';
    h.gain.MarkerSize = 20;
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-300 120];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 14;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -300:60:120;
        
 	[~,h.phase1] = PlotPatch(rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx)), rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx)), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4}, 1, 1, 'k', [0.4 0.4 0.6], 0.5, 2);
    h.phase1.Marker = '.';
    h.phase1.MarkerSize = 20;
    
	[~,h.phase2] = PlotPatch(PHASE_Mean, PHASE_STD, Freq, 1, 1, 'g', [0.4 0.4 0.6], 0.5, 2);
    h.phase2.Marker = '.';
    h.phase2.MarkerSize = 20;
    
    leg = legend([h.phase1 h.phase2],'Gain/Phase','Complex');
    leg.Box = 'off';
    
    plot(ax1.XLim,[0 0],'--k','LineWidth',1);

end