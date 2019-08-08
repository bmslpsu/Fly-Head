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
clearvars -except HeadFree

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

Real        = real(CmplxGain);
Imag        = imag(CmplxGain);

for jj = 1:nFreq
    for kk = 1:size(Real,1)
        if (jj==5) && ( Real(kk,jj)<0 && Imag(kk,jj)>0 )
            Imag(kk,jj) = -Imag(kk,jj);
            Real(kk,jj) = -Real(kk,jj);
        end
    end
end

CmplxGain   = Real + 1i*Imag;

Gain        = abs(CmplxGain);
Phase       = rad2deg(angle((CmplxGain)));

GAIN        = mean(Gain,1);
PHASE       = rad2deg(circ_mean(deg2rad(Phase),[],1));
GAIN_STD    = std(Gain,[],1);
PHASE_STD   = std(Phase,[],1);

REAL        = mean(Real,1);
IMAG        = mean(Imag,1);

% CMPLXGAIN   = REAL + 1i*IMAG;
% GAIN        = abs(CMPLXGAIN);
% PHASE       = rad2deg(angle((CMPLXGAIN)));

cList = prism(nFreq);

gains = 0.2:0.2:1;

%% Complex Gain
%---------------------------------------------------------------------------------------------------------------------------------
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

% ax.XTick = '';
% ax.YTick = '';

for jj = 1:nFreq
    h.trial = scatter(Real(:,jj), Imag(:,jj), 15, 'o','MarkerEdgeColor','k',...
      'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.65, 'LineWidth', 0.5);
end

for jj = 1:nFreq
  	h.rr = plot([0 REAL(jj)],[0 IMAG(jj)],'Color',[0 0 0 1],'LineWidth',1);
    
    rSTD = PolarSTD(Real(:,jj),Imag(:,jj),[REAL(jj) IMAG(jj)]);
    
	[h.std] = draw_ellipse([REAL(jj) IMAG(jj)], 3*rSTD, 0.5, 0, 90 - PHASE(jj), cList(jj,:)); hold on
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
leg = legend(hh,num2strcell(Freq));
leg.Title.String = 'Frequency (Hz)';
leg.Location = 'northwest';
legend boxoff

%% BODE
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
FIG.Name = filename;
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
    
    errorbar(Freq,GAIN,GAIN_STD,'-ob','LineWidth',3);
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-300 60];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 14;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -300:60:60;
    
    errorbar(Freq,PHASE,PHASE_STD,'-ob','LineWidth',3);
    
    plot(ax1.XLim,[0 0],'--g','LineWidth',2);

end