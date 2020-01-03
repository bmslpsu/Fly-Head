function [FIG] = MakeFig_Ramp_HeadWing_Saccade_Manual()
%% MakeFig_Ramp_HeadWing_Saccade:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\30\Vid\tracked\SACD';

% Select files
[files, PATH] = uigetfile({'*mat', 'files'}, 'Select files',root, 'MultiSelect', 'on');
FILES = cellstr(files)';

[D,I,N,U,T,~,~,basename] = GetFileData(FILES,false,'fly','trial','vel','wave');

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

PP = nan(N.file,1);
pp = 1;
debug = true;
window = 0.1; % search wiondow [s]
FIG.Visible = 'off';
for jj = 1:N.file
    sacd = load(fullfile(PATH,FILES{jj}),'SACD');
    
    HEAD = sacd.SACD.Head;
    WING = sacd.SACD.Wing;
    head_count = size(HEAD,1);

%     [Head.SACD,Head.thresh,Head.count] = SacdDetect(HEAD.X(:,1),HEAD.Time,350,debug);
%     Head.svel = HEAD.X(:,2);
%     Head.svel( (Head.svel<Head.thresh)  & (Head.svel>0))  = 0;
%     Head.svel( (Head.svel>-Head.thresh) & (Head.svel<0))  = 0;
% 
%     [Wing.SACD,Wing.thresh,Wing.count] = SacdDetect(WING.X(:,1),WING.Time,1.75,debug);
%     Wing.svel = WING.X(:,2);
%     Wing.svel( (Wing.svel<Wing.thresh)  & (Wing.svel>0) )  = 0;
%     Wing.svel( (Wing.svel>-Wing.thresh) & (Wing.svel<0) )  = 0;
% 
%     FIG.Visible = 'off';
%     axes(ax.L) ; cla ; hold on
%     plot(HEAD.Time, Head.svel,'Color', ax.L.YColor)
%     plot(Head.SACD.PeakTime, Head.SACD.PeakVel, '*', 'Color', 'b', 'MarkerSize', 15)
% 
%     axes(ax.R) ; cla ; hold on
%     plot(WING.Time, Wing.svel,'Color', ax.R.YColor)
%     plot(Wing.SACD.PeakTime, Wing.SACD.PeakVel, '*', 'Color', 'r', 'MarkerSize', 15)


    if ~isnan(HEAD{1,1})
        for ww = 1:head_count
            tpeak = HEAD.StartTime(ww);
            tstart = HEAD.StartTime(ww);
            cnd = (WING.PeakTime>=(tpeak-window)) & (WING.PeakTime<=(tpeak+window));

            if sum(cnd)==0
                PP(pp,1) = false;
            elseif sum(cnd)==1
                PP(pp,1) = true;
                plot(WING.PeakTime(cnd),WING.PeakVel(cnd), 'o', 'Color', 'g', 'MarkerSize', 20)
                tdiff = WING.StartTime(cnd) - tstart;
                match = sign(WING.PeakVel(cnd)) == sign(HEAD.PeakVel(ww));
                PP(pp,2) = tdiff;
                PP(pp,3) = match;
                if abs(tdiff)>0.1
                    % pause()
                end
            elseif sum(cnd)>1
                PP(pp,1) = false;
                disp('More than one saccade in window')
%                 pause()
            else
                error('WUTT')
            end

        pp = pp + 1;
        end
    end
%     pause()
    close all
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