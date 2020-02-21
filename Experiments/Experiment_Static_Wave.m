function [] = Experiment_Static_Wave(Fn)
%% Experiment_Static_Wave: runs an experiment using LED arena and fly panel controller
%
%  Written for Panel Controller v3 and NiDAQ seesion mode
%
%   INPUTS:
%       Fn      :  	fly number
%
daqreset
imaqreset

%% Set directories & experimental parameters
root = 'C:\BC\Rigid_data\Experiment_Static_Wave';

% EXPERIMENTAL PARAMETERS
n_tracktime = 10 + 1;       % length(func)/fps; seconds for each EXPERIMENT
n_resttime = 1;             % seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 6;                  % number of cycles through spatial frequencies for each fly
patID = 2;                  % Spatial frequency grating pattern
FPS = 200;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Gain = 750;                	% camera gain
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 1:6;                	% Analog input channels
AO = 0;                     % Analog output channels

%% Set up data acquisition & camera
[s,~] = NI_USB_6212(Fs,AI,AO);

% Camera Trigger Signal
t = 0:1/s.Rate:n_tracktime;
TRIG = ((1/2)*(square(2*pi*FPS*t,50) - 1)');
TRIG(TRIG==-1) = 4;

[vid,~] = Basler_acA640_120gm(FPS,Gain,nFrame);

%% Set variable to control pattern spatial frequency
wave = 7.5*[0 3 4 8 inf nan]';  % [deg] spatial frequencies
n_wave = length(wave);         	% # of velocities
ypos = [1 4 5 7 12 1];       	% pattern y-pos corresponding to spatial frequencies

% Create sequence of randomly shuffled frequencies
Freq_all = nan(n_wave*n_rep,1);
pp = 0;
for kk = 1:n_rep
    Freq_rand = wave(randperm(n_wave),:);    % reshuffle randomly
    Freq_all(pp+1:pp+n_wave,1) = Freq_rand;  % add rep
    pp = kk*n_wave;
end

% y-pos index vector
ypos_all = Freq_all;
for kk = 1:length(ypos)
    ypos_all(ypos_all==wave(kk)) = ypos(kk);
end

n_trial = n_rep * n_wave;

%% EXPERIMENT LOOP
disp('Start Experiment:')
disp('Wavelengths:')
disp(Freq_all')
disp('--------------------------------------------------')
for ii = 1:n_trial
    disp('Trial')
    disp(num2str(ii));  % print counter to command line
    preview(vid);       % open video preview window
	start(vid)          % start video buffer
    
    % CLOSED LOOP BAR TRACKING %
    Arena_CL(1,'X',-15)
    pause(n_resttime)
    Panel_com('stop')
    
    pause(1) % pause between closed-loop & experiment 
    
    % EXPERIMENT SETUP %
    disp('Play Stimulus: ')
    disp(['Wavelength = ' num2str(Freq_all(ii))])
    
    if isnan(Freq_all(kk))
        Panel_com('set_pattern_id', 1);pause(n_pause)	% set pattern
    else
        Panel_com('set_pattern_id', 2);pause(n_pause)	% set pattern
    end
    
    Panel_com('set_pattern_id', patID); pause(n_pause)    	% set pattern
    Panel_com('set_position',[randi([1,96]) , ...        	% set starting position (xpos,ypos)
                  	ypos_all(kk)]); pause(n_pause)           
	Panel_com('set_funcX_freq', 50); pause(n_pause)      	% update rate for x-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)     	% update rate for y-channel
    Panel_com('set_mode', [0,0]); pause(n_pause)         	% 0=open,1=closed,2=fgen,3=vmode,4=pmode
	
    % START EXPERIMENT & DATA COLLECTION %
    queueOutputData(s,TRIG) % set trigger AO signal
    T = timer('StartDelay',0.5,'TimerFcn',@(src,evt) Panel_com('start'));
    start(T)
    tic
        [data, t_p ] = s.startForeground; % data collection
        Panel_com('stop') % stop stimulus
        stop(vid) % stop video buffer
        [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    toc
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    % CLOSED LOOP BAR TRACKING %
    Arena_CL(1,'X',-15)
    
    % SAVE DATA %
    disp('Saving...')
    disp('----------------------------------------------------------------------')
    
    fname = ['Fly_' num2str(Fn) '_Trial_' num2str(ii) '_Wave_' num2str(Freq_all(ii)) ...
        '_Vel_' num2str(0) '.mat'];
    
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
end

delete(vid)
disp('Done');
daqreset
imaqreset
PControl
end