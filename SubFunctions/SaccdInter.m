function [Saccade,Interval,Stimulus,Error,IntError] = SaccdInter(xx,tt,SACD,match,STIM,debug)
%% SaccdInter: extracts saccade & inter-saccade intervals
%   INPUTS:
%       xx          :   position
%       tt          :   time vector
%       SACD      	:   saccade table
%       match      	:   remove saccades of a specified direction
%       Stim      	:   stimulus data (same size as "xx")
%       debug       :   showplot boolean
%   OUTPUTS:
%       Saccade     :   saccade intervals
%       Interval   	:   inter-saccade intervals
%       Stimulus   	:   stimulus intervals
%       Error       :   error intervals
%       IntError  	:   integrated error intervals
%---------------------------------------------------------------------------------------------------------------------------------
% Check if there are any saccades in the data
if isnan(SACD.Duration(1))
    SACD(1,:) = [];
    warning('No saccades detected')
end

if nargin>3
  	matchFlag = false;
    % Remove saccades in specified direction
    if ~isnan(match)
        [anti] = find(SACD.Match~=match);
        if ~isempty(anti)
            SACD(anti,:) = [];
        end
        
        if isempty(SACD)
            warning('No saccades fit the match condition')
            matchFlag = true;
        end
        
    end
end

xx = xx(:);
tt = tt(:);
Ts  = mean(diff(tt));           % sampling time
vel(:,1) = [0 ; diff(xx)/Ts];  	% derivative of data

if ~isempty(STIM) % we know the stimulus
    nn = 4;
    vel(:,2) = [0 ; diff(STIM)/Ts];	% derivative of stimulus
else
    nn = 2;
end

n = size(SACD,1);               % # of saccades
Inter = cell(n,nn);          	% store saccde & inter-sacade intervals
Time = cell(n,2);            	% store saccde & inter-sacade time intervals
PeakTime = cell(n,1);           % peak velocity time for each saccade
StartTime = cell(n,1);          % start times for each inter-saccade interval
EndTime = cell(n,1);            % start times for each inter-saccade interval
for kk = 1:n % all saccades
    % Saccade Intervals
    SI = SACD.StartIdx(kk);             % saccade interval start index
    EI = SACD.EndIdx(kk);               % saccade interval end index
    Inter{kk,1}(:,1) = xx(SI:EI);       % saccade interval
 	Inter{kk,1}(:,2) = vel(SI:EI);      % inter-saccade interval
    Time{kk,1} = tt(SI:EI);             % saccade interval times
    PeakTime{kk} = SACD.PeakTime(kk);   % peak velocity time
    StartTime{kk} = SACD.StartTime(kk);	% start times for each inter-saccade interval (based on 
    EndTime{kk} = SACD.EndTime(kk);     % start times for each inter-saccade interval
    if ~isempty(STIM) % we know the stimulus
        Inter{kk,3}(:,1) = STIM(SI:EI);
        Inter{kk,3}(:,2) = vel(SI:EI,2);
    end
    
    % Inter-Saccde Intervals
    if kk==1 % 1st saccade does not have valid inter-saccade interval
    	Time{kk,2}(:,1) = tt(SI:EI);
      	Inter{kk,2}(:,1:2) = nan(length(Time{kk,2}(:,1)),2);
        if ~isempty(STIM) % we know the stimulus
            Inter{kk,4}(:,1:2) = nan(length(Time{kk,2}(:,1)),2);
        end
    else
        SI = SACD.EndIdx(kk-1);             % inter-saccade interval start index
        EI = SACD.StartIdx(kk);             % inter-saccade interval end index
        Inter{kk,2}(:,1) = xx(SI:EI);       % inter-saccade interval
    	Inter{kk,2}(:,2) = vel(SI:EI,1);    % inter-saccade velocity interval
    	Time{kk,2}(:,1) = tt(SI:EI);        % inter-saccade interval time
        if ~isempty(STIM) % we know the stimulus
            Inter{kk,4}(:,1) = STIM(SI:EI);
            Inter{kk,4}(:,2) = vel(SI:EI,2);
        end        
    end
end

TimeNorm = cellfun(@(x) x - x(1), Time, 'UniformOutput', false); % normalized times
TimePeak = cellfun(@(x,y) x - y, Time(:,1), PeakTime, 'UniformOutput', false); % aligned to peak times of saccades
TimeEnd = cellfun(@(x,y) x - y, Time(:,2), StartTime, ... % aligned to end time of inter-saccde intervals (saccde start)
    'UniformOutput', false); 

[Saccade.Time,~,~,~,dR,~] = nancat_center(TimePeak,0,1);

maxInter = cellfun(@(x) length(x), Inter(:,2), 'UniformOutput', true);
[Interval.Time,~,~,~,~,~] = nancat_center(TimeNorm(:,2),[],1,1);

Saccade.Pos    = nan(size(Saccade.Time));
Saccade.Vel    = nan(size(Saccade.Time));
Interval.Pos   = nan(size(Interval.Time));
Interval.Vel   = nan(size(Interval.Time));

Stimulus.Saccade.Pos 	= nan(size(Saccade.Time));
Stimulus.Saccade.Vel 	= nan(size(Saccade.Time));
Stimulus.Interval.Pos 	= nan(size(Interval.Time));
Stimulus.Interval.Vel 	= nan(size(Interval.Time));

for kk = 1:n
    Saccade.Pos(:,kk)      = cat_pad(Inter{kk,1}(:,1), dR{kk},nan);
	Saccade.Vel(:,kk)      = cat_pad(Inter{kk,1}(:,2), dR{kk},nan);
    
    interOff = [0,max(maxInter)] - [0,maxInter(kk)];
    
    Interval.Pos(:,kk)     = cat_pad(Inter{kk,2}(:,1), interOff, nan);
	Interval.Vel(:,kk)     = cat_pad(Inter{kk,2}(:,2), interOff, nan);
    
    Interval.Pos(:,kk) = Interval.Pos(:,kk) - Interval.Pos(1,kk);
    
    if ~isempty(STIM) % we know the stimulus
        Stimulus.Saccade.Pos(:,kk)      = cat_pad(Inter{kk,3}(:,1), dR{kk},nan);
        Stimulus.Saccade.Vel(:,kk)      = cat_pad(Inter{kk,3}(:,2), dR{kk},nan); 
        Stimulus.Interval.Pos(:,kk)     = cat_pad(Inter{kk,4}(:,1), interOff, nan);
        Stimulus.Interval.Vel(:,kk)   	= cat_pad(Inter{kk,4}(:,2), interOff, nan);
        
        Stim_Sacd_Pos_Norm = Stimulus.Saccade.Pos(~isnan(Stimulus.Saccade.Pos(:,kk)),kk);
        Stim_Sacd_Pos_Norm = Stim_Sacd_Pos_Norm(1);
        Stim_Inter_Pos_Norm = Stimulus.Interval.Pos(~isnan(Stimulus.Interval.Pos(:,kk)),kk);
        if all(isnan(Stim_Inter_Pos_Norm))
            Stim_Inter_Pos_Norm = nan;
        else
            Stim_Inter_Pos_Norm = Stim_Inter_Pos_Norm(1);
        end
        
        Stimulus.Saccade.Pos(:,kk)     = Stimulus.Saccade.Pos(:,kk) - Stim_Sacd_Pos_Norm;
        Stimulus.Interval.Pos(:,kk)    = Stimulus.Interval.Pos(:,kk) - Stim_Inter_Pos_Norm;
    end
    
end

Error.Saccade.Pos = Stimulus.Saccade.Pos - Saccade.Pos;
Error.Saccade.Vel = Stimulus.Saccade.Vel - Saccade.Vel;
Error.Interval.Pos = Stimulus.Interval.Pos - Interval.Pos;
Error.Interval.Vel = Stimulus.Interval.Vel - Interval.Vel;

IntError.Saccade.Pos = nancumsum(Error.Saccade.Pos,1,2);
IntError.Saccade.Vel = nancumsum(Error.Saccade.Vel,1,2);
IntError.Interval.Pos = nancumsum(Error.Interval.Pos,1,2);
IntError.Interval.Vel = nancumsum(Error.Interval.Vel,1,2);

Sacd_Pos_Mean = nanmean(Saccade.Pos,2);
Sacd_Vel_Mean = nanmean(Saccade.Vel,2);
Sacd_Time_Mean = nanmean(Saccade.Time,2);
Inter_Pos_Mean = nanmean(Interval.Pos,2);
Inter_Vel_Mean = nanmean(Interval.Vel,2);
Inter_Time_Mean = nanmean(Interval.Time,2);

if matchFlag
    return
end

if debug
    FIG = figure ; clf % data with saccade & inter-saccde intervals highlighted
    FIG.Color = 'k';
        clear ax
        ax(1) = subplot(2,1,1) ; hold on
        ax(1).XLabel.String = 'Time';
        ax(1).YLabel.String = 'Data';
        
        plot(tt,xx,'--w','LineWidth',0.5)
        for kk = 1:n
            h.sacd  = plot(Time{kk,1},Inter{kk,1}(:,1),'-','LineWidth',2);
            h.inter = plot(Time{kk,2},Inter{kk,2}(:,1),'-w','LineWidth',0.5);
        end
%         legend([h.sacd h.inter],'Saccade','Interval')

        ax(2) = subplot(2,1,2) ; hold on
        ax(2).XLabel.String = 'Time';
        ax(2).YLabel.String = 'Diff Data';
        plot(tt,vel,'--w','LineWidth',0.5)
        for kk = 1:n
            h.sacd  = plot(Time{kk,1},Inter{kk,1}(:,2),'-','LineWidth',2);
            h.inter = plot(Time{kk,2},Inter{kk,2}(:,2),'-w','LineWidth',0.5);
        end
        
        cc = 'w';
        set(ax,'Color','k','YColor',cc,'XColor',cc)
    
	FIG = figure ; clf % normalized saccade & inter-saccde intervals
    FIG.Color = 'k';
        clear ax
        ax(1) = subplot(2,2,1) ; hold on
        ax(1).YLabel.String = 'Saccade';
     	h.sacd  = plot(Sacd_Time_Mean,Saccade.Pos,'-','LineWidth',1);
        plot(Sacd_Time_Mean,Sacd_Pos_Mean,'-w','LineWidth',3);
        plot([0 0],ax(1).YLim,'--w','LineWidth',1.5);
                
        ax(2) = subplot(2,2,2) ; hold on
        ax(2).YLabel.String = 'Interval-Sacadde Interval';
       	h.inter = plot(Interval.Time,Interval.Pos,'-','LineWidth',1);
        plot(Inter_Time_Mean,Inter_Pos_Mean,'-w','LineWidth',3);
        
        if ~isempty(STIM)
            plot(Interval.Time,Stimulus.Interval.Pos,'--r','LineWidth',2)
        end
        
        ax(3) = subplot(2,2,3) ; hold on
        ax(3).XLabel.String = 'Normalized Time';
        ax(3).YLabel.String = 'Saccade Diff';
       	h.sacd  = plot(Saccade.Time,Saccade.Vel,'-','LineWidth',1);
        plot(Sacd_Time_Mean,Sacd_Vel_Mean,'w','LineWidth',3);
        plot([0 0],ax(3).YLim,'--w','LineWidth',1.5);

        ax(4) = subplot(2,2,4) ; hold on
        ax(4).XLabel.String = 'Normalized Time';
        ax(4).YLabel.String = 'Interval-Sacadde Interval Diff';
       	h.inter = plot(Interval.Time,Interval.Vel,'-','LineWidth',1);
       	plot(Inter_Time_Mean,Inter_Vel_Mean,'-w','LineWidth',3);
        
        if ~isempty(STIM)
            plot(Interval.Time,Stimulus.Interval.Vel,'--r','LineWidth',2)
        end

        cc = 'w';
        set(ax,'Color','k','YColor',cc,'XColor',cc)
end

end