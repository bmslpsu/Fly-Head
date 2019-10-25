function [FIG] = MakeFig_Ramp_HeadWing_CrossCorr()
%% MakeFig_Ramp_HeadWing_CrossCorr:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'SACD','U','N','TRIAL','GRAND');

%%
clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND

clear FIG ax
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 3];
movegui(FIG,'center')
FIG.Name = 'Head vs Wing';

headIdx = 2;
wingIdx = 3;
h2w = 6;
pp = 1;
TD = cell(N{1,3},1);
ax = gobjects(N{1,3},1);
for jj = 1:N{1,3}
    ax(jj) = subplot(2,3,jj); hold on ; title(['Vel = ' num2str(3.75*U{1,3}{1}(jj)) ' (' char(176) '/s)'])
    xlabel('Time (ms)')
    ylabel('Cross Correlation')
    for kk = 1:N{1,1}
        for ii = 1:size(TRIAL{kk,jj},1)
            td = 1000*TRIAL{kk,jj}{ii,h2w}.TimeLags;
            cc = TRIAL{kk,jj}{ii,h2w}.CrossCorr;
            
            t_diff = 1000*TRIAL{kk,jj}{ii,h2w}.TimeDiff;
            max_cc = TRIAL{kk,jj}{ii,h2w}.MaxCC;
            
            TD{jj}(end+1,1) = t_diff;
            
%             plot(td, cc, 'k', 'LineWidth', 0.25)
%             plot(t_diff, max_cc, 'r.', 'MarkerSize', 10)
            
            pp = pp + 1;
        end
    end
    td = 1000*GRAND{jj,h2w}.Mean{2}{11};
    cc = GRAND{jj,h2w}.Mean{2}{10};
    cc_std = GRAND{jj,h2w}.STD{2}{10};
    [cc_peak,idx_peak] = max(cc);
    
    PlotPatch(cc, cc_std, td, 1, N{1,1}, 'b', [0.4 0.4 0.6], 0.2, 2);
    plot([td(idx_peak) , td(idx_peak)], [cc_peak, 0], 'r--')
	plot(td(idx_peak), cc_peak, 'r.', 'MarkerSize', 15)
end
set(ax,'XLim',1000*0.2*[-1 1])
set(ax,'YTickLabels',{})
% set(ax,'XTick',-150:50:150)
linkaxes(ax,'xy')

%%
clear FIG ax
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 6 3];
movegui(FIG,'center')
FIG.Name = 'Head vs Wing';
FIG.Visible = 'off';

ax.L = gca ; cla ; hold on
ax.L.FontSize = 12;
ax.L.Color = 'none';
ax.L.YColor = [0 0 0.7];
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel(['Head (' char(176) '/s)'],'FontSize',12,'FontWeight','bold')
ax.L.YLim = 1000*[-1 1];
ax.L.XLim = [0 10];

ax.R = axes;
ax.R.YAxisLocation = 'right';
ax.R.Color = 'none';
ax.R.YColor = [0.7 0 0];
ax.R.YLabel.String = 'Wing (V/s)';
ax.R.YLabel.FontSize = 12;
ax.R.YLabel.FontWeight = 'bold';
ax.R.XLim = ax.L.XLim;
ax.R.YLim = 200*[-1 1];
ax.R.XTickLabels = '';
ax.R.XColor = 'none';
ax.R.Position = ax.L.Position;
linkaxes([ax.L ax.R],'x');

PP = nan(size(SACD.Head,1),3); pp = 1;
debug = false;
window = 0.1; % search wiondow [s]
FIG.Visible = 'off';
for jj = 1:N{1,3}
    for kk = 1:N{1,1}
        for ii = 1:size(TRIAL{kk,jj},1)
         	HEAD = TRIAL{kk,jj}{ii,headIdx};
            WING = TRIAL{kk,jj}{ii,wingIdx};
            
          	[Head.SACD,Head.thresh,Head.count] = SacdDetect(HEAD.X(:,1),HEAD.Time,2.5,debug);
            Head.svel = HEAD.X(:,2);
            Head.svel( (Head.svel<Head.thresh)  & (Head.svel>0))  = 0;
            Head.svel( (Head.svel>-Head.thresh) & (Head.svel<0))  = 0;
            
         	[Wing.SACD,Wing.thresh,Wing.count] = SacdDetect(WING.X(:,1),WING.Time,1.75,debug);
            Wing.svel = WING.X(:,2);
            Wing.svel( (Wing.svel<Wing.thresh)  & (Wing.svel>0) )  = 0;
            Wing.svel( (Wing.svel>-Wing.thresh) & (Wing.svel<0) )  = 0;
            
            FIG.Visible = 'off';
            axes(ax.L) ; cla ; hold on
            plot(HEAD.Time, Head.svel,'Color', ax.L.YColor)
            plot(Head.SACD.PeakTime, Head.SACD.PeakVel, '*', 'Color', 'b', 'MarkerSize', 15)

            axes(ax.R) ; cla ; hold on
            plot(WING.Time, Wing.svel,'Color', ax.R.YColor)
            plot(Wing.SACD.PeakTime, Wing.SACD.PeakVel, '*', 'Color', 'r', 'MarkerSize', 15)
            
            
            if ~isnan(Head.count)
                for ww = 1:Head.count
                    tpeak = Head.SACD.StartTime(ww);
                    tstart = Head.SACD.StartTime(ww);
                    cnd = (Wing.SACD.PeakTime>=(tpeak-window)) & (Wing.SACD.PeakTime<=(tpeak+window));

                    if sum(cnd)==0
                        PP(pp,1) = false;
                    elseif sum(cnd)==1
                        PP(pp,1) = true;
                        plot(Wing.SACD.PeakTime(cnd), Wing.SACD.PeakVel(cnd), 'o', 'Color', 'g', 'MarkerSize', 20)
                        tdiff = Wing.SACD.StartTime(cnd) - tstart;
                        match = sign(Wing.SACD.PeakVel(cnd)) == sign(Head.SACD.PeakVel(ww));
                        PP(pp,2) = tdiff;
                        PP(pp,3) = match;
                        if abs(tdiff)>0.1
                            % pause()
                        end
                    elseif sum(cnd)>1
                        PP(pp,1) = false;
                        disp('More than one saccade in window')
                        pause()
                    else
                        error('WUTT')
                    end

                pp = pp + 1;
                end
            end
            % pause()
        end
    end
end

disp('Linked saccade ratio:')
disp(sum(PP(:,1))./size(PP,1))
disp('Median time difference (s):')
disp(nanmedian(PP(:,2)))
disp('Mathed direction saccade ratio:')
disp(nansum(PP(:,3))./size(PP,1))

clear FIG ax
FIG = figure (2) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 5 3];
movegui(FIG,'center')
FIG.Name = 'Head vs Wing Time Difference';
FIG.Visible = 'on';

h = histogram(1000*PP(:,2), 1000*(-0.15:0.005:0.15),'FaceColor',[0.1 0 0.7],'FaceAlpha',0.6);
ax = gca;
ax.FontSize = 8;
xlabel('Time Difference (ms)','FontSize',8)
ylabel('Count','FontSize',8)

end