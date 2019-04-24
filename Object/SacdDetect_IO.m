function [] = SacdDetect(time,fly,pat,N)
%% SacdDetect: Reads in all raw trials, transforms data, and saves in organized structure
%   INPUTS:
%       time  	: time vector
%       fly    	: fly angle
%       pat   	: pattern angle
%       N   	: STD detection threshold
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
%%
% N = 2.5;

%%
%---------------------------------------------------------------------------------------------------------------------------------
time.data       = time;                          	% time vector
n               = length(time.data);                % # of data points
Ts              = mean(diff(time));             	% sampling time
Fs              = 1/Ts;                             % sampling frequency
time.norm    	= (0:Ts:time.data(end))';       	% normalized time

pat.pos.data 	= 60*time.data;                     % stimulus
pat.vel.data  	= [0 ; diff(pat.pos.data)/Ts];    	% pattern velocity
pat.avel.data  	= abs(pat.vel.data);                % pattern absolute velocity

[b,a]           = butter(2,25/(Fs/2));              % low-pass filter coefficents
pos.data        = filtfilt(b,a,fly);                % position
vel.data        = [0 ; diff(pos.data)/Ts];          % velocity
avel.data       = abs(vel.data);                    % absolute velocity
vel.mean      	= mean(avel.data);                  % mean absolute velocity
vel.std       	= std(vel.data);                    % std absolute velocity
thresh          = vel.mean + N*vel.std;             % saccade detetcion threshold
svel            = avel.data;                        % preallocate velocity data above threshold

svel(avel.data<thresh) = 0;                         % velocity data above threshold

% Find velcocity spikes
[avel.pks, loc.pks] = findpeaks(svel,'MINPEAKDISTANCE',25); % find local maxima
window      = n*0.015;  % window length (in samples)
I           = find(loc.pks > (window) & loc.pks < ((length(vel.data) - (window+1)))); % ignore saccades at beginning and end
loc.pks   	= loc.pks(I);
avel.pks  	= avel.pks(I);
vel.pks    	= avel.pks.*sign(vel.data(loc.pks));
pos.pks     = pos.data(loc.pks);
time.pks    = time.data(loc.pks);

% SACD = [];
for ww = 1:length(loc.pks) % every saccade
    % Find start & end velocity
    loc.start(ww) = find(avel.data(1:loc.pks(ww)) <= avel.pks(ww)/4,1,'last'); % saccade start index               
    Eind = find(avel.data <= avel.pks(ww)/4); % all values below 1/4 peak
    Es = find(Eind > loc.pks(ww),1,'first'); % first value after the start index is the end index
    if ~isempty(Es) % make sure data did not start above threshold
        loc.end(ww) = Eind(Es); % saccade end index
    end

    % Inter-saccade interval
    if ww ~= 1 % disregard first saccade interval
        % Raw
        loc.inter          = loc.end(ww-1):loc.start(ww);
        time.inter          = time.data(loc.inter);
        pos.inter           = pos.data(loc.inter);
        vel.inter           = vel.data(loc.inter);
        avel.inter          = avel.data(loc.inter);
        pat.pos.inter     	= pat.pos.data(loc.inter);
        pat.vel.inter     	= pat.vel.data(loc.inter);
     	pat.avel.inter     	= pat.avel.data(loc.inter);
        % Normalized
     	loc.Ninter         = loc.inter - loc.inter(1);
        time.Ninter         = time.inter - time.inter(1);
        pos.Ninter          = pos.inter  - pos.inter(1);
        vel.Ninter          = vel.inter  - vel.inter(1);
        avel.ninter         = avel.inter  - avel.inter(1);
        pat.pos.Ninter     	= pat.pos.inter - pat.pos.inter(1);
        pat.vel.Ninter     	= pat.vel.inter - pat.vel.inter(1);
     	pat.avel.Ninter   	= pat.avel.inter - pat.avel.inter(1);
        % Error
        err.pos             = pat.pos.Ninter - pos.Ninter;
        err.vel             = pat.vel.inter - vel.inter;
        % Integrated error
    	err_int.pos         = cumtrapz(time.Ninter,err.pos);
        err_int.vel         = cumtrapz(time.Ninter,err.vel);
        % Derivative of error
    	err_dx.pos      	= diff(err.pos)/Ts;
        err_dx.vel          = diff(err.vel)/Ts;

        % Pre-window error
        pre_err.pos         = mean(err.pos(end-round(n*0.005)));
        pre_err.vel       	= mean(err.vel(end-round(n*0.005)));
        pre_err_int.pos  	= mean(err_int.pos(end-round(n*0.005)));
        pre_err_int.vel 	= mean(err_int.vel(end-round(n*0.005)));
     	pre_err_dx.pos  	= mean(err_dx.pos(end-round(n*0.005)));
        pre_err_dx.vel      = mean(err_dx.vel(end-round(n*0.005)));
    end
end

avel.start  	= avel.data(loc.start);
avel.end        = avel.data(loc.end);
vel.start       = vel.data(loc.start);
vel.end         = vel.data(loc.end);
pos.start       = pos.data(loc.start);
pos.end         = pos.data(loc.end);
time.start   	= time.data(loc.start);
time.end      	= time.data(loc.end);
pos.amp         = pos.end - pos.start;

figure (1) ; clf
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
% plot(time.data,pat.data,'c');
ylabel('Position')
xlabel('Time')


end