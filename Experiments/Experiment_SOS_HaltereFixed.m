function [] = Experiment_SOS_HaltereFixed(Fn)
%% Experiment_SOS_HaltereFixed: runs an experiment using LED arena and fly panel controller
% Written for Panel Controller v3 and NiDAQ seesion mode
%   INPUTS:
%       Fn          :  	fly number
%---------------------------------------------------------------------------------------------------------------------------------
daqreset
imaqreset
%% Set directories & experimental paramters %%
%---------------------------------------------------------------------------------------------------------------------------------
%rootdir = uigetdir({}, 'Select folder to save data'); % define directory to save file
root = 'C:\BC\Rigid_data\Experiment_SOS_v2_HaltereFixBoth';

% EXPERIMENTAL PARAMETERS
n_tracktime = 20 + 1;       % length(func)/fps; seconds for each EXPERIMENT
n_resttime = 1;             % seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_trial = 20;               % # of repetitions
patID = 2;                  % Spatial frequency grating pattern
yPos  = 5;                  % 30 deg spatial frequency
funcX = 1;                  % SOS (20s)
xUpdate = 200;              % function update rate
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Gain = 750;                	% camera gain
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 1:6;                	% Analog input channels
AO = 0;                     % Analog output channels

%% Set up data acquisition & camera %%
%---------------------------------------------------------------------------------------------------------------------------------
[s,~] = NI_USB_6212(Fs,AI,AO);

% Camera Trigger Signal
t = 0:1/s.Rate:n_tracktime;
TRIG = ((1/2)*(square(2*pi*FPS*t,50) - 1)');
TRIG(TRIG==-1) = 4;

[vid,~] = Basler_acA640_120gm(FPS,Gain,nFrame);

%% EXPERIMENT LOOP %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Start Experiment:')
for ii = 1:n_trial   
    disp('Trial')
    disp(num2str(ii));  % print counter to command line
    preview(vid);       % open video preview window
	start(vid)          % start video buffer
    %-----------------------------------------------------------------------------------------------------------------------------
    % CLOSED LOOP BAR TRACKING %
    Arena_CL(1,'X',-15)
    pause(n_resttime)
    Panel_com('stop')
    %-----------------------------------------------------------------------------------------------------------------------------
    pause(1) % pause between closed-loop & experiment 
    % EXPERIMENT SETUP %
    disp('Play Stimulus: ')
    Panel_com('set_pattern_id', patID); pause(n_pause)           	% set pattern
    Panel_com('set_position',[15 , yPos]); pause(n_pause)           % set starting position (xpos,ypos)
    Panel_com('set_posfunc_id',[funcX, 1]); pause(n_pause)       	% arg1 = channel (x=1,y=2); arg2 = funcid
	Panel_com('set_funcX_freq', xUpdate); pause(n_pause)            % update rate for x-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)              	% update rate for y-channel
    Panel_com('set_mode', [4,0]); pause(n_pause)                    % 0=open,1=closed,2=fgen,3=vmode,4=pmode
	%-----------------------------------------------------------------------------------------------------------------------------
    % START EXPERIMENT & DATA COLLECTION %
    queueOutputData(s,TRIG) % set trigger AO signal
    T = timer('StartDelay',0.1,'TimerFcn',@(src,evt) Panel_com('start'));
    start(T)
    tic
        [data, t_p ] = s.startForeground; % data collection
        Panel_com('stop') % stop stimulus
        stop(vid) % stop video buffer
        [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    toc
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    trig = logical(round(data(:,1)));
    syncIdx  = find(trig==true,1,'first');
    syncTime = t_p(syncIdx);
    t_sync = t_p - syncTime;
    
    %-----------------------------------------------------------------------------------------------------------------------------
    % CLOSED LOOP BAR TRACKING %
    Arena_CL(1,'X',-15)
    %-----------------------------------------------------------------------------------------------------------------------------
    % SAVE DATA %
    disp('Saving...') ; disp('----------------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_SOS_HaltereFixed_Both.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v','Fs','syncTime','t_sync');
end

delete(vid)
disp('Done');
daqreset
imaqreset
PControl
end