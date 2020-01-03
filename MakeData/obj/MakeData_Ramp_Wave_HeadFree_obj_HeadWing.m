function [] = MakeData_Ramp_Wave_HeadFree_obj_HeadWing(Fc,wave)
%% MakeData_Ramp_Wave_HeadFree_obj_HeadWing: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       -
%

Fc = 30;
wave = 30;
filename = ['Ramp_HeadFree_SACCD_HeadWing_filt=' num2str(Fc) '_Wave=' num2str(wave)];
rootdir = ['H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

%% Setup Directories %%
root.daq = rootdir;
root.vid = fullfile(root.daq,'\Vid');
root.ang = fullfile(root.daq,'\Vid\tracked');

% Select files
[files, PATH.ang] = uigetfile({'*csv', 'files'}, 'Select files',root.ang, 'MultiSelect', 'on');
FILES = cellstr(files)';

[D,I,N,U,T,~,~,basename] = GetFileData(FILES,false,'fly','trial','vel','wave');

PATH.daq = root.daq;
PATH.vid = root.vid;
PATH.sacd = fullfile(root.ang,'SACD');

%% Get Data %%
disp('Loading...')
clear SACD
badtrial = {};

Vel = U.vel{1};
tt = (0:(1/200):10)'; 
Stim = (Vel*tt')';
bad = 1;
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    disp(basename{kk})
   	
    % Get head data
    head.Time = t_v;
    head.Pos = rad2deg(benifly.Head);
    head.Fc = Fc;
    Head = Fly(head.Pos,head.Time,head.Fc,[],tt); % head object
    if max(abs(Head.X(3:end,2)))>900
        badtrial{bad,1} = basename{kk};
        bad = bad + 1;
        % error('Head Trial Error')
    end
  	
    % Get wing data from DAQ
	wing.f          = medfilt1(100*data(:,6),3); % wing beat frequency [Hz]
    wing.Time       = t_p; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = Fc; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
  	wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
   	Wing            = Fly(wing.Pos,t_p,15,[],Head.Time); % wing object
    
    wing.f          = interp1(wing.Time,wing.f,Head.Time);
   	wing.Left     	= interp1(wing.Time,wing.Left,Head.Time);
   	wing.Right     	= interp1(wing.Time,wing.Right,Head.Time);

    Wing.WBF        = wing.f;
    Wing.WBA        = [wing.Left,wing.Right,wing.Left + wing.Right];
	
  	% Check WBF & WBA
    if min(wing.f)<150 || mean(wing.f)<200 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        % continue
    elseif any(wing.Left>10.6) || any(wing.Right>10.6)
        fprintf('WBA out of range: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        % continue
    end
	% Get Saccade Stats   
    % [head.SACD,head.thresh,head.count,head.rate,head.SACDRmv] = SacdDetect(Head.X(:,1),Head.Time,2.5,true);
    % [wing.SACD,wing.thresh,wing.count,wing.rate,wing.SACDRmv] = SacdDetect(Wing.X(:,1),Wing.Time,1.75,true);
    
    [SACD.Head,~,~,~,~,~,~] = Sacd_Manual(Head.X(:,1),Head.Time,false);
    pause
 	[SACD.Wing,~,~,~,~,~,~] = Sacd_Manual(Wing.X(:,1),Wing.Time,false);
    pause
    close all
    
    save(fullfile(PATH.sacd, [basename{kk} '.mat']), 'SACD')
    
end
disp('DONE')
end