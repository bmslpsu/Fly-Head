function [SACD,thresh,count,rate] = SacdDetect(xx,tt,N,debug)
%% SacdDetect: Reads in all raw trials, transforms data, and saves in organized structure
%   INPUTS:
%       xx    	: position
%       tt  	: time vector
%       N   	: STD detection threshold (default: N=2.75)
%   OUTPUTS:
%       SACD    : saccade data table
%       thresh 	: saccade detetcion threshold
%       count 	: # of saccades
%       rate 	: saccades/s
%---------------------------------------------------------------------------------------------------------------------------------
n = length(xx); % # of data points
IDX = (1:n)'; % index vector
if nargin<3
    N = 2.75; % default STD detection threshold
    if nargin<2
        tt = IDX; % use index if no time data
    end
end

time.data       = tt;                               % time vector
Ts              = mean(diff(time.data));        	% sampling time
time.new        = (min(time.data):Ts:max(time.data))'; % new time vector for interpolation
pos.data        = xx;                               % position
vel.data        = [0 ; diff(pos.data)/Ts];          % velocity
avel.data       = abs(vel.data);                    % absolute velocity
vel.mean      	= mean(avel.data);                  % mean absolute velocity
vel.std       	= std(vel.data);                    % std absolute velocity
thresh          = vel.mean + N*vel.std;             % saccade detetcion threshold
svel            = avel.data;                        % preallocate velocity data above threshold

svel(avel.data<thresh) = 0;                         % velocity data above threshold

% Variable names
varnames = {'Duration','Amplitude','Direction','StartIdx','PeakIdx','EndIdx','StartTime','PeakTime','EndTime','StartPos',...
                'PeakPos','EndPos','StartVel','PeakVel','EndVel','StartAbsVel','PeakAbsVel','EndAbsVel'};

% Find velocity spikes
[avel.pks, loc.pks] = findpeaks(svel,'MINPEAKDISTANCE',40); % find local maxima
window      = n*0.015;  % window length (in samples)
I           = find(loc.pks > (window) & loc.pks < ((length(vel.data) - (window+1)))); % ignore saccades at start & end
loc.pks   	= loc.pks(I);
avel.pks  	= avel.pks(I);
vel.pks    	= avel.pks.*sign(vel.data(loc.pks));
pos.pks     = pos.data(loc.pks);
time.pks    = time.data(loc.pks);
count       = length(loc.pks);
rate        = count/(time.data(end)-time.data(1));
boundThresh = 1/8;
if ~isempty(I)
    for ww = 1:count % every saccade
        % Find start & end velocity
        loc.start(ww,1) = find(avel.data(1:loc.pks(ww)) <= avel.pks(ww)/4,1,'last'); % saccade start index               
        Eind = find(avel.data <= avel.pks(ww)*boundThresh); % all values below 1/4 peak
        Es = find(Eind > loc.pks(ww,1),1,'first'); % first value after the start index is the end index
        if ~isempty(Es) % make sure data did not start above threshold
            loc.end(ww,1) = Eind(Es); % saccade end index
        else
            loc.end(ww,1) = 1;
        end

    end
    % Stats
    avel.start  	= avel.data(loc.start);
    avel.end        = avel.data(loc.end);
    vel.start       = vel.data(loc.start);
    vel.end         = vel.data(loc.end);
    pos.start       = pos.data(loc.start);
    pos.end         = pos.data(loc.end);
    time.start   	= time.data(loc.start);
    time.end      	= time.data(loc.end);
    pos.amp         = pos.end - pos.start;
    time.dur        = time.end - time.start;
    dir             = sign(vel.pks);
    
    % All data
    DATA = [time.dur , pos.amp , dir, loc.start , loc.pks , loc.end , time.start , time.pks , time.end , pos.start , pos.pks, ...
                pos.end ,vel.start , vel.pks , vel.end , avel.start , avel.pks , avel.end];
	
else % no saccades
    DATA    = nan(1,length(varnames));
    count   = nan;
    rate    = nan;
end

% Saccade table
SACD = splitvars(table(DATA));
SACD.Properties.VariableNames = varnames;

% Remove Saccades
time.rmv = time.data;
pos.rmv = pos.data;
vel.rmv = vel.data;
if ~isempty(I)
    loc.cell_idx    = cell(count,1);
    loc.cell_all    = cell(count,1);
    time.cell_rmv   = cell(count,1);
 	pos.cell_rmv    = cell(count,1);
 	vel.cell_rmv    = cell(count,1);
    for ww = 1:count % every saccade
        loc.cell_idx{ww}    = (loc.start(ww):loc.end(ww))';
        loc.cell_all{ww}    = IDX( ismember( IDX, loc.cell_idx{ww} ) );
        time.cell_rmv{ww}   = time.data(loc.cell_all{ww});
       	pos.cell_rmv{ww}    = pos.data(loc.cell_all{ww});
       	vel.cell_rmv{ww}    = vel.data(loc.cell_all{ww});
    end
    loc.idx = cat(1,loc.cell_idx{:});
    loc.all = ismember(IDX, loc.idx );
    
    time.rmv(loc.idx) = nan;
    pos.rmv(loc.idx)  = nan;
    vel.rmv(loc.idx)  = nan;
    
    time.nan        = time.rmv;
    
%     time.rmv(loc.idx) = time.data(loc.all);
%     pos.rmv(loc.idx)  = pos.data(loc.all);
%     vel.rmv(loc.idx)  = vel.data(loc.all);
    
    pos.rmv_n = pos.rmv;
    
    for ww = 1:count
        pos.rmv_n(loc.end(ww)-1:end) = pos.rmv_n(loc.end(ww)-1:end) - (pos.end(ww) - pos.start(ww));
    end
  	pos.nan         = pos.rmv_n;
    pos.rmv_n       = pos.rmv_n(~isnan(pos.rmv_n));
    time.rmv        = time.rmv(~isnan(time.rmv));
    pos.rmv_intrp   = interp1(time.rmv,pos.rmv_n,time.new,'linear');

else
    pos.rmv_n = pos.data;
end

% Debug plots
if debug
    FIG = figure ; clf
    FIG.Color = 'w';
    
    ax1 = subplot(3,1,1) ; hold on
        ylabel('Velocity')
        plot(time.data,vel.data,'k')
        plot(time.data, zeros(n,1)*thresh,'--','Color',[0.5 0.5 0.5])
        plot(time.data, ones(n,1)*thresh,'--m')
        plot(time.data,-ones(n,1)*thresh,'--m')
        for kk = 1:count
            plot(time.cell_rmv{kk},vel.cell_rmv{kk},'-c','MarkerSize',5)
        end
        ax1.YLim = max(abs(ax1.YLim))*[-1 1];

        plot(time.pks,vel.pks,'*b')
        plot(time.start,vel.start,'*g')
        plot(time.end,vel.end,'*r')
    
    ax2 = subplot(3,1,2) ; hold on
        ylabel('Position')
        plot(time.data,pos.data,'k')
        plot(time.data, zeros(n,1)*thresh,'--','Color',[0.5 0.5 0.5])
        for kk = 1:count
            plot(time.cell_rmv{kk},pos.cell_rmv{kk},'-c','MarkerSize',5)
        end
        plot(time.pks,pos.pks,'*b')
        plot(time.start,pos.start,'*g')
        plot(time.end,pos.end,'*r')
        ax2.YLim = max(abs(ax2.YLim))*[-1 1];
        
    ax3 = subplot(3,1,3) ; hold on
        plot(time.data,pos.rmv_intrp,'c','LineWidth',1)
        plot(time.nan ,pos.nan,'k','LineWidth',2)
        ylabel('Position (Removed Saccades)')
        xlabel('Time')
end


%         % Inter-saccade interval
%         if ww ~= 1 % disregard first saccade interval
% %             disp(ww)
%             % Raw
%             loc.inter           = loc.end(ww-1):loc.start(ww);
%             time.inter          = time.data(loc.inter);
%             pos.inter           = pos.data(loc.inter);
%             vel.inter           = vel.data(loc.inter);
%             avel.inter          = avel.data(loc.inter);
%             % Normalized
%             loc.Ninter          = loc.inter - loc.inter(1);
%             time.Ninter         = time.inter - time.inter(1);
%             pos.Ninter          = pos.inter  - pos.inter(1);
%             vel.Ninter          = vel.inter  - vel.inter(1);
%             avel.ninter         = avel.inter  - avel.inter(1);
%         end
end