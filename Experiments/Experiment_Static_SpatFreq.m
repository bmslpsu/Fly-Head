function [] = Experiment_Static_SpatFreq(Fn)
% Experiment_Static_SpatFreq: runs an experiment using LED arena and fly panel controller
% This code is written for Panel Controller v3 and NiDAQ seesion mode
% NOTES:    1. Static background with variable spatial frequency
%           2. User must maunally set the root folder in the fucntion
%   INPUTS:
%       Fn      :  	fly number
%---------------------------------------------------------------------------------------------------------------------------------
daqreset
imaqreset
%% Set Directories & Controller Parameters %%
%---------------------------------------------------------------------------------------------------------------------------------
% rootdir = ['D:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\LowContrast\' num2str(spatFreq) '\'];
rootdir = 'D:\EXPERIMENTS\Experiment_Static_SpatFreq\';
viddir = [rootdir 'Vid\'];
%% EXPERIMENTAL PARAMETERS %%
%---------------------------------------------------------------------------------------------------------------------------------
n.tracktime = 10;  	% seconds for each EXPERIMENT
n.resttime = 2;    	% seconds for each REST
n.pause = 0.2;      % pause between commands
n.AI = 6;        	% # of analog input channels
n.rep = 6;          % number of cycles through velocities for each fly

%% SETUP DATA AQUISITION: NiDAQ (session mode)
%---------------------------------------------------------------------------------------------------------------------------------
% Find DAQ
devices = daq.getDevices;
s = daq.createSession('ni');
addAnalogInputChannel(s,devices.ID, 1:n.AI, 'Voltage'); % Add Analog Input Channels
    % AI1 = AO_0: trigger
    % AI2 = DAC1: x-position of stimulus
    % AI3 = DAC2: y-position of stimulus
    % AI4 = L: left wing beat amplitude
    % AI5 = R: right wing beat amplitude
    % AI6 = FREQ: frequency
    
% Add analog output channels
ch.AO = addAnalogOutputChannel(s,devices.ID,'ao0','Voltage');
ch.AO.Name = 'Trigger';

% Set Sampling Rate
s.Rate = 5000;                   % samples per second
s.IsContinuous = false;          % continuous data collection until stopped
% s.DurationInSeconds = 10;

% Setup Sampling
s.Rate = 5000; % samples per second
s.IsContinuous = false;	% continuous data collection until stopped

FrameRate = 200; % camera frame rate
nFrame = FrameRate * n.tracktime; % # of frames to log

t = 0:1/s.Rate:n.tracktime;
TriggerSignal = (square(2*pi*FrameRate*t,90) + 1)';

disp('DAQ Setup Done...')
%% SETUP CAMERA INPUT %%
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
srcObj1.GainRaw = 964;
srcObj1.ExposureTimeAbs = 10000; % 100 Hz frame rate
% srcObj1(1).ExposureMode = 'Timed'; % exposure time controlled by pulse width

% Trigger config
triggerconfig(vid, 'hardware','DeviceSpecific','DeviceSpecific')
set(srcObj1(1),'LineSelector','Line1');
set(srcObj1(1),'TriggerActivation','RisingEdge');
set(srcObj1(1),'TriggerMode','on');
set(srcObj1(1),'TriggerSelector','FrameStart');

disp('VID Setup Done...')
%% Set variable to control pattern spatial frequency %%
%---------------------------------------------------------------------------------------------------------------------------------
freq = 7.5*[0 3 4 8 inf]';  	% [deg] spatial frequencies
n.freq = length(freq);      	% # of velocities
ypos = [1 4 5 7 12];          	% pattern y-pos corresponding to spatial frequencies

% Create sequence of randomly shuffled frequencies
Freq_all = nan(n.freq*n.rep,1);
pp = 0;
for kk = 1:n.rep
    Freq_rand = freq(randperm(n.freq),:);    % reshuffle randomly
    Freq_all(pp+1:pp+n.freq,1) = Freq_rand;  % add rep
    pp = kk*n.freq;
end

% y-pos index vector
ypos_all = Freq_all;
for kk = 1:length(ypos)
    ypos_all(ypos_all==freq(kk)) = ypos(kk);
end

%% START EXPERIMENT AND DATA COLLECTION %%
%---------------------------------------------------------------------------------------------------------------------------------
tic
for kk = 1:n.rep*n.freq
    disp('-------------------------------------------------------')
    disp(['Trial:  ' num2str(kk)]) % prints counter to command line
   	preview(vid); % open video preview window
	start(vid) % start video buffer
    %-----------------------------------------------------------------------------------------------------------------------------
    % CLOSED LOOP BAR TRACKING
    disp('rest');
    Panel_com('stop'); pause(n.pause)
    Panel_com('set_pattern_id', 1);pause(n.pause)               % set pattern to "Pattern_Fourier_bar_barwidth=8"
    Panel_com('set_position',[1, 1]); pause(n.pause)            % set start position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n.pause)                 % 0=open,1=closed,2=fgen,3=vmode,4=pmode)
    Panel_com('set_funcX_freq', 50); pause(n.pause)             % default X update rate
	Panel_com('set_funcY_freq', 50); pause(n.pause)           	% default Y update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n.pause)     % [xgain,xoffset,ygain,yoffset]
    Panel_com('start');                                         % start closed-loop tracking
    pause(n.resttime)
    Panel_com('stop');
    %-----------------------------------------------------------------------------------------------------------------------------
    % SETUP EXPERIMENT
    pause(1)
    disp(['Spatial Frequency ' num2str(Freq_all(kk)) ' Hz'])
    Panel_com('stop');pause(n.pause);
    Panel_com('set_pattern_id', 2);pause(n.pause)                       % set pattern
    Panel_com('set_position',[randi([1,96]), ypos_all(kk)]);pause(n.pause) 	% set starting position (xpos,ypos) [ypos = spatFreq]
    Panel_com('set_funcX_freq', 50);pause(n.pause);                         % default X update rate
    Panel_com('set_funcY_freq', 50);pause(n.pause);                         % default Y update rate
    Panel_com('set_mode', [0,0]);pause(n.pause)                             % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    Panel_com('send_gain_bias',[0,0,0,0]);pause(n.pause)                    % no gain
    %-----------------------------------------------------------------------------------------------------------------------------
    % RUN EXPERIMENT AND COLLECT DATA
    queueOutputData(s,TriggerSignal) % set trigger AO signal
    Panel_com('start')  % run trial
    [data,t_p] = s.startForeground;
	stop(vid) % stop video buffer
    Panel_com('stop')
    %-----------------------------------------------------------------------------------------------------------------------------
    % GET DATA AND SAVE TO .mat FILE
	[vidData,t_v] = getdata(vid, vid.FramesAcquired); % get video data
      
    % CLOSED LOOP BAR TRACKING
    disp('rest');
    Panel_com('stop'); pause(n.pause)
    Panel_com('set_pattern_id', 1);pause(n.pause)               % set pattern to "Pattern_Fourier_bar_barwidth=8"
    Panel_com('set_position',[1, 1]); pause(n.pause)            % set start position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n.pause)                 % 0=open,1=closed,2=fgen,3=vmode,4=pmode)
    Panel_com('set_funcX_freq', 50); pause(n.pause)             % default X update rate
	Panel_com('set_funcY_freq', 50); pause(n.pause)           	% default Y update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n.pause)     % [xgain,xoffset,ygain,yoffset]
    Panel_com('start');                                         % start closed-loop tracking
    
    % Save data
    disp('Saving...')
    filename  = ['Fly_' num2str(Fn) '_Trial_' num2str(kk) '_SpatFreq_' num2str(Freq_all(kk)) '_Vel_' num2str(0) '.mat'];
    save(fullfile(rootdir,filename),'-v7.3','data','t_p');
    save(fullfile(viddir,filename),'-v7.3','vidData','t_v');
    %-----------------------------------------------------------------------------------------------------------------------------
end
toc
%---------------------------------------------------------------------------------------------------------------------------------
disp('Done')
daqreset
imaqreset
end