function [SACD,thresh,count,rate] = SacdDetect(fly,tt,N,debug)
%% SacdDetect: Reads in all raw trials, transforms data, and saves in organized structure
%   INPUTS:
%       tt  	: time vector
%       fly    	: fly angle
%       pat   	: pattern angle
%       N   	: STD detection threshold
%   OUTPUTS:
%       SACD    : saccade data table
%       thresh 	: saccade detetcion threshold
%       count 	: # of saccades
%       rate 	: saccades/s
%---------------------------------------------------------------------------------------------------------------------------------
%%
% N = 2.5;
% debug = true;
%%
%---------------------------------------------------------------------------------------------------------------------------------
time.data       = tt;                               % time vector
n               = length(time.data);                % # of data points
Ts              = mean(diff(time.data));        	% sampling time

pos.data        = fly;                              % position
vel.data        = [0 ; diff(pos.data)/Ts];          % velocity
avel.data       = abs(vel.data);                    % absolute velocity
vel.mean      	= mean(avel.data);                  % mean absolute velocity
vel.std       	= std(vel.data);                    % std absolute velocity
thresh          = vel.mean + N*vel.std;             % saccade detetcion threshold
svel            = avel.data;                        % preallocate velocity data above threshold

svel(avel.data<thresh) = 0;                         % velocity data above threshold

% Find velcocity spikes
[avel.pks, loc.pks] = findpeaks(svel,'MINPEAKDISTANCE',40); % find local maxima
window      = n*0.015;  % window length (in samples)
I           = find(loc.pks > (window) & loc.pks < ((length(vel.data) - (window+1)))); % ignore saccades at beginning and end
loc.pks   	= loc.pks(I);
avel.pks  	= avel.pks(I);
vel.pks    	= avel.pks.*sign(vel.data(loc.pks));
pos.pks     = pos.data(loc.pks);
time.pks    = time.data(loc.pks);
count       = length(loc.pks);
rate        = count/(time.data(end)-time.data(1));
if ~isempty(I)
    for ww = 1:length(loc.pks) % every saccade
        % Find start & end velocity
        loc.start(ww,1) = find(avel.data(1:loc.pks(ww)) <= avel.pks(ww)/4,1,'last'); % saccade start index               
        Eind = find(avel.data <= avel.pks(ww)/4); % all values below 1/4 peak
        Es = find(Eind > loc.pks(ww,1),1,'first'); % first value after the start index is the end index
        if ~isempty(Es) % make sure data did not start above threshold
            loc.end(ww,1) = Eind(Es); % saccade end index
        end

        % Inter-saccade interval
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
else
    % Stats
    avel.start  	= [];
    avel.end        = [];
    vel.start       = [];
    vel.end         = [];
    pos.start       = [];
    pos.end         = [];
    time.start   	= [];
    time.end      	= [];
    pos.amp         = [];
    time.dur        = [];
    dir             = [];   
    
    avel.pks        = [];
 	vel.pks         = [];
    pos.pks         = [];
    time.pks        = [];
end


% All data
DATA = [time.dur , pos.amp , dir, loc.start , loc.pks , loc.end , time.start , time.pks , time.end , pos.start , pos.pks , pos.end ,...
        vel.start , vel.pks , vel.end , avel.start , avel.pks , avel.end];

% Variable names
varnames = {'Duration','Amplitude','Direction','StartIdx','PeakIdx','EndIdx','StartTime','PeakTime','EndTime','StartPos','PeakPos','EndPos'...
           	'StartVel','PeakVel','EndVel','StartAbsVel','PeakAbsVel','EndAbsVel'};

% Saccade table
SACD = splitvars(table(DATA));
SACD.Properties.VariableNames = varnames;


% Debug plots
if debug
    figure ; clf
    subplot(2,1,1) ; hold on
    ylabel('Velocity')
    plot(time.data,vel.data,'k')
    plot(time.data, zeros(n,1)*thresh,'--','Color',[0.5 0.5 0.5])
    plot(time.data, ones(n,1)*thresh,'--m')
    plot(time.data,-ones(n,1)*thresh,'--m')
    plot(time.pks,vel.pks,'*b')
    plot(time.start,vel.start,'*g')
    plot(time.end,vel.end,'*r')

    subplot(2,1,2) ; hold on
    plot(time.data,pos.data,'k')
    plot(time.data, zeros(n,1)*thresh,'--','Color',[0.5 0.5 0.5])
    plot(time.pks,pos.pks,'*b')
    plot(time.start,pos.start,'*g')
    plot(time.end,pos.end,'*r')
    ylabel('Position')
    xlabel('Time')
end
end