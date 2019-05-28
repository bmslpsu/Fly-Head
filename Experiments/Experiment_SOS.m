function [] = Experiment_SOS(Fn)
%% Experiment_SOS: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
%% Set directories & experimental paramters %%
%---------------------------------------------------------------------------------------------------------------------------------
root = 'D:\EXPERIMENTS\Experiment_SOS_v2\';
% viddir = [root 'Vid\'];

% EXPERIMENTAL PARAMETERS
n_tracktime = 21;	% length(func)/fps; seconds for each EXPERIMENT
n_resttime = 1;     % seconds for each REST
n_pause = 0.2;      % seconds for each pause between panel commands
n_trial = 20;     	% # of repetitions
n_AI = 6;          	% # of analog input channels

%% Set up data acquisition on NiDAQ (session mode) %%
%---------------------------------------------------------------------------------------------------------------------------------
devices = daq.getDevices; % get DAQ ID
s = daq.createSession('ni'); % create session

% Add analog input channels
ch.AI = addAnalogInputChannel(s,devices.ID, 1:n_AI, 'Voltage');
% Assign AI channel names
chNames = {'Trigger','Pat x-pos','Pat y-pos','L-WBA','R-WBA','WB-Freq'};
for kk = 1:length(ch.AI)
    ch.AI(kk).Name = chNames{kk};
        % AI1 = AO_0: trigger
        % AI2 = DAC1: x-position of stimulus
        % AI3 = DAC2: y-position of stimulus
        % AI4 = L: left wing beat amplitude
        % AI5 = R: right wing beat amplitude
        % AI6 = FREQ: frequency
end

% Add analog output channels
ch.AO = addAnalogOutputChannel(s,devices.ID,'ao0','Voltage');
ch.AO.Name = 'Trigger';

% Setup Sampling
s.Rate = 5000; % samples per second
s.IsContinuous = false;	% continuous data collection until stopped

FrameRate = 100; % camera frame rate
nFrame = FrameRate * n_tracktime; % # of frames to log

t = 0:1/s.Rate:n_tracktime;
TriggerSignal = (square(2*pi*FrameRate*t,90) + 1)';
% plot(t,TriggerSignal)

%% Setup camera input to sample at Fs = 200Hz %%
%---------------------------------------------------------------------------------------------------------------------------------
adaptorName = 'gige';
deviceID = 1;
vidFormat = 'Mono8';
vid = videoinput(adaptorName, deviceID, vidFormat);

% Configure vidobj properties.
set(vid, 'ErrorFcn', @imaqcallback);
set(vid, 'LoggingMode', 'memory');
set(vid,'FrameGrabInterval',1);
ROI.x = 500;
ROI.y = 250;
ROI.xoff = (round(659 - ROI.x)/2);
ROI.yoff = (round(494 - ROI.y)/2);
vid.ROIPosition = [ROI.xoff ROI.yoff ROI.x ROI.y];
set(vid, 'FramesPerTrigger', 1); % frames to log for each trigger
vid.TriggerRepeat = nFrame-1; % # triggers

% Configure vidobj source properties.
srcObj1 = get(vid, 'Source');
srcObj1.Gamma = 0.386367797851563;
srcObj1.GainRaw = 750;
% srcObj1.ExposureTimeAbs = 8000;
% srcObj1(1).ExposureMode = 'Timed'; % exposure time controlled by pulse width

% Trigger config
triggerconfig(vid, 'hardware','DeviceSpecific','DeviceSpecific')
set(srcObj1(1),'LineSelector','Line1');
set(srcObj1(1),'TriggerActivation','RisingEdge');
set(srcObj1(1),'TriggerMode','on');
set(srcObj1(1),'TriggerSelector','FrameStart');

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
    disp('rest');
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_pattern_id', 1);pause(n_pause)           % set output to p_rest_pat (Pattern_Fourier_bar_barwidth=8)
    Panel_com('set_position',[1, 1]); pause(n_pause)        % set starting position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n_pause)             % closed loop tracking [xpos,ypos] (NOTE: 0=open, 1=closed)
    Panel_com('set_posfunc_id',[1,0]); pause(n_pause)       % no position function for x-channel
    Panel_com('set_funcX_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('set_posfunc_id',[2,0]); pause(n_pause)       % no position function for y-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n_pause)	% [xgain,xoffset,ygain,yoffset]
    Panel_com('start');                                     % start closed-loop rest
    pause(n_resttime)                                       % rest 5 seconds
    Panel_com('stop');                                      % stop closed-loop rest
    %-----------------------------------------------------------------------------------------------------------------------------
    pause(1) % pause between closed-loop & experiment 
    % EXPERIMENT SETUP %
    disp('Play Stimulus: ')
    Panel_com('set_pattern_id', 2); pause(n_pause)                	% set pattern to "Pattern_spatFreq_22_30_60"
    Panel_com('set_position',[15 , 4]); pause(n_pause)             	% set starting position (xpos,ypos)
    Panel_com('set_posfunc_id',[1, 1]); pause(n_pause)              % arg1 = channel (x=1,y=2); arg2 = funcid
	Panel_com('set_funcX_freq', 200); pause(n_pause)                % 100Hz update rate for x-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)              	% 50Hz update rate for y-channel
    Panel_com('set_mode', [4,0]); pause(n_pause)                    % 0=open,1=closed,2=fgen,3=vmode,4=pmode
	%-----------------------------------------------------------------------------------------------------------------------------
    % START EXPERIMENT & DATA COLLECTION %
    queueOutputData(s,TriggerSignal) % set trigger AO signal
    tic
        Panel_com('start') 	% start stimulus
        [data, t_p ] = s.startForeground; % data collection
        Panel_com('stop')  	% stop stimulus
        stop(vid) % stop video buffer
        [vidData t_v] = getdata(vid, vid.FramesAcquired); % get video data
        disp(['Fs = ' num2str(1/mean(diff(t_v)))]) % check Fs of video
    toc
    %-----------------------------------------------------------------------------------------------------------------------------
    % CLOSED LOOP BAR TRACKING %
    disp('rest');
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_pattern_id', 1);pause(n_pause)           % set output to p_rest_pat (Pattern_Fourier_bar_barwidth=8)
    Panel_com('set_position',[1, 1]); pause(n_pause)        % set starting position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n_pause)             % closed loop tracking [xpos,ypos] (NOTE: 0=open, 1=closed)
    Panel_com('set_posfunc_id',[1,0]); pause(n_pause)       % no position function for x-channel
    Panel_com('set_funcX_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('set_posfunc_id',[2,0]); pause(n_pause)       % no position function for y-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n_pause)	% [xgain,xoffset,ygain,yoffset]
    Panel_com('start');                                     % start closed-loop rest
    %-----------------------------------------------------------------------------------------------------------------------------
    % SAVE DATA %
    disp('Saving...') ; disp('----------------------------------------------------------------------')
    save([root   'fly_' num2str(Fn) '_trial_' num2str(ii)  '_SOS.mat'],'-v7.3','data','t_p','vidData','t_v');
%     save([viddir 'fly_' num2str(Fn) '_trial_' num2str(ii)  '_SOS.mat'],'-v7.3','vidData','t_v');
    %-----------------------------------------------------------------------------------------------------------------------------
end

delete(vid)
disp('Done');

clear all
daqreset
imaqreset
%---------------------------------------------------------------------------------------------------------------------------------
