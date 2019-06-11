function [] = MakeFig_Sine_HeadFree_ComplexGain()
%% MakeFig_Sine_HeadFree_ComplexGain:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILES,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nAmp = length(FILES);
Amp = nan(nAmp,1);
for ww = 1:nAmp
filedata = textscan(FILES{ww}, '%s', 'delimiter', '_');
    Amp(ww) = str2double(filedata{1}{3});
end

HeadFree = cell(nAmp,1);
for ww = 1:nAmp
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N');
end

filename = 'Sine_HeadFree_ComplexGain';

catIdx = 5;
xIdx = 1;
figNum = 1;

Freq = HeadFree{ww}.U{1,3}{1};
nFreq = HeadFree{ww}.N{1,3};
vel = nan(nFreq,1);
velIdx = [4,4,3,2,1,1];
Vel = repmat(2*pi*3.75*HeadFree{1}.U{1,3}{1},1,nAmp);
for kk = 1:nAmp
   Vel(:,kk) = round(kk*Vel(:,kk));
end
for jj = 1:nFreq
    vel(jj) = Vel(jj,velIdx(jj));
end
vLabel = cellfun(@(x) strcat(num2str(x),[' ' char(176) '/s']), num2cell(vel),'UniformOutput',false);
cList = jet(nFreq);
cLabel = cellfun(@(x) [num2str(x) ' Hz  : '], num2cell(Freq),'UniformOutput',false);

legLabel = cellfun(@(x,y) [x ' ' y], cLabel, vLabel,'UniformOutput',false);

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end
hold on

ax = gca;
box on
axis equal
ax.Color = 'w';
ax.YLabel.String = 'Img';
ax.YLabel.FontSize = 14;
ax.XLabel.String = 'Real';
ax.XLabel.FontSize = ax.YLabel.FontSize;

if catIdx==5
    gains = 0.2:0.2:1;
else
    gains = 0.05:0.05:0.25;
end
maxGain = max(gains);

axisVector = 2*(-maxGain:0.1:maxGain);
h.xorigin = plot(axisVector,0*axisVector,'--','Color',[0.5 0.5 0.5],'LineWidth',1);
h.yorigin = plot(0*axisVector,axisVector,'--','Color',[0.5 0.5 0.5],'LineWidth',1);
for kk = 1:length(gains)
    h.circle = PlotCircle(0,0,gains(kk));
    h.circle.Color = [0.5 0.5 0.5 1];
    h.circle.LineStyle = '-';
    h.circle.LineWidth = 1;
    h.text = text(0,1.06*gains(kk),num2str(gains(kk)));
    h.text.FontWeight = 'bold';
    if gains(kk)==1
        h.circle.Color = [0.8 0 0 0.5];
    end
end
ax.XLim = (1.1*maxGain)*[-1 1];
ax.YLim = (1.1*maxGain)*[-1 1];


REAL = cell(nAmp,1);
IMAG = REAL;
for ww = 1:nAmp
    REAL{ww} = nan(100,nFreq);
    for jj = 1:nFreq
        pp = 1;
        for kk = 1:HeadFree{ww}.N{1,1}
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1)
                R = real(HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.IOCmplxGain(:,xIdx));
                I = imag(HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.IOCmplxGain(:,xIdx));
                
                if jj==1 && I<0
                    I = abs(I);
                    R = abs(R);
                elseif jj==5 && (I>0 || R<0)
                 	I = -abs(I);
                    R = abs(R);
                end
                
                REAL{ww}(pp,jj) = R;
                IMAG{ww}(pp,jj) = I;
                
                if ww==velIdx(jj)
                    h.trial = scatter(R,I,40,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
                                'MarkerFaceAlpha',0.65,'LineWidth',0.5);
                end
                
%                   plot([0 R],[0 I],'Color',[0.5 0.5 0.5 0.2])
%                   h.trial = plot(R,I,'.','Color',cList(jj,:),'MarkerSize',20);
%                   h.trial.Color(4) = 0.2;
                pp = pp + 1;
            end
        end
    end
end

REAL = cellfun(@(x) nanmean(x,1), REAL, 'UniformOutput', false);
IMAG = cellfun(@(x) nanmean(x,1), IMAG, 'UniformOutput', false);
REAL = cat(1,REAL{:});
IMAG = cat(1,IMAG{:});
for jj = 1:nFreq
    RR(jj) = REAL(velIdx(jj),jj);
    II(jj) = IMAG(velIdx(jj),jj);
end
REAL = RR;
IMAG = II;
for ww = 1:1
    for jj = 1:nFreq
%         R = real(HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{18}(:,xIdx));
%         I = imag(HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{18}(:,xIdx));
%         R = REAL{ww}(jj);
%         I = IMAG{ww}(jj);
        R = REAL(jj);
     	I = IMAG(jj);
        
       	plot([0 R],[0 I],'Color',[0 0 0 1],'LineWidth',1)
      	h.grand = scatter(R,I,120,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
                        'MarkerFaceAlpha',1,'LineWidth',1.5);
                    
        hh(jj) = h.grand;
%         h.text = text(R,I,vLabel{jj});
%         h.text.FontWeight = 'bold';
    end
end

leg = legend(hh,legLabel);
leg.Title.String = 'Frequency / Velocity';
leg.Location = 'northeast';
legend boxoff

CmplxGain = REAL + 1i*IMAG;
GAIN = 2*abs(CmplxGain);
PHASE = rad2deg(angle(CmplxGain));

figure (2) ; clf
plot(Freq,PHASE)

end