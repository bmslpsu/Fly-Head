function [PAT,WING,HEAD,BODE,FD,T,n,unq] = MakeData_Chirp_HeadFree_TEST(rootdir,filename)
%% MakeData_Chirp_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       PAT     : pattern structure
%       WING   	: wings structure
%       HEAD  	: head structure
%       FD      : file data
%       T       : fly data table
%       n       : field #'s 
%       unq   	: unique fields
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
rootdir = 'F:\EXPERIMENTS\Experiment_ChirpLog_HeadFree\';
% filename = 'Chirp_HeadFree_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;
root.ang = [root.daq '\Vid\Angles\'];

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = FILES';

PATH.daq = root.daq;

[~,I,N,U] = GetFileData(FILES);

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
% Preallocate data cells
% WING.ALL.Pos = cell(n.Amp,1);
% HEAD.ALL.Pos = cell(n.Amp,1);

% Store data in cells
for kk = 1:N{1,4}
    disp(kk)
    % Load HEAD & DAQ data %
    data = [];
	load([PATH.daq   FILES{kk}],'data','t_p'); % load pattern x-position
    load([PATH.ang   FILES{kk}],'hAngles','t_v'); % load head angles % time arrays
    %-----------------------------------------------------------------------------------------------------------------------------
    % Check WBF
	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',FD.Fly(kk),FD.Trial(kk))
        continue
    end
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get head data %
    head.Time = t_v;
    head.Pos = hAngles - mean(hAngles);
    Head = Fly(head.Pos,head.Time,20); % head object
  	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ %
    wing.Time       = t_p; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = 20; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
  	wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
    
   	Wing = Fly(wing.Pos,t_p,20,Head.time); % wing object
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ %
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,Head.time); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.time,0.4*Head.Fs); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between hea d& pattern %
    head.Err.Pos    = Pat.x   - Head.x; % calculate position error between head & pattern [deg]
    head.Err.Vel    = Pat.d1x - Head.d1x; % calculate position Error between head & pattern (retinal slip) [deg/s]
    
    PosErr = Fly(head.Err.Pos,Head.time,0.4*Head.Fs); % pattern object
	VelErr = Fly(head.Err.Vel,Head.time,0.4*Head.Fs); % pattern object

    % Calculate coherence %
    [head.cohr.mag,head.cohr.f] = mscohere(pat.Pos , head.Pos ,[],[] , head.Freq , head.Fs);
	[wing.cohr.mag,wing.cohr.f] = mscohere(pat.Pos , wing.Pos ,[],[] , wing.Freq , wing.Fs);
    %-----------------------------------------------------------------------------------------------------------------------------
	% Calculate BODE gain & phase difference for head & wings %
    head.GAIN   = medfilt1(head.Mag./pat.Mag,5);
    head.PHASE  = medfilt1(-(pat.Phase - head.Phase),5);
 	wing.GAIN   = medfilt1(wing.Mag./head.Err.Mag,5);
    wing.PHASE  = medfilt1(-(head.Err.Phase - wing.Phase),5);
    bode.head2wing.GAIN = medfilt1(wing.Mag./head.Mag,5);
    bode.head2wing.PHASE = medfilt1(-(head.Phase - wing.Phase),5);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
  	% PATTERN
% 	PAT.Time        {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.Time;
	

end
clear jj kk a b t_p t_v hAngles data head wing pat bode tt
disp('LOADING DONE')


end