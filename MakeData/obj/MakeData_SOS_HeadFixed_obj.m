function [] = MakeData_SOS_HeadFixed_obj(rootdir)
%% MakeData_SOS_HeadFixed_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       rootdir    	:   root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = 'H:\EXPERIMENTS\Experiment_SOS_v2_HeadFixed';
filename = 'SOS_HeadFixed_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;

% Select files
[FILES, PATH.daq] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.daq, 'MultiSelect','on');
FILES = cellstr(FILES)';

[D,I,N,U,T] = GetFileData(FILES);

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
% IOFreq = 0.1*round(linspace(0.1,8,10)/0.1)';
IOFreq = [1, 3.1, 5.3, 7.4, 9.6];
disp('Loading...')
ALL 	= cell([N{1,end},5]); % cell array to store all data objects
TRIAL  	= cell(N{1,1},1);
n.catg  = size(N,2) - 1;
pp = 0;
span = 1:2100;
tt = linspace(0,20,2000)';
for kk = 1:N{1,end}
    disp(kk)
    % Load DAQ data
    data = [];
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
  	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ
	wing.f          = medfilt1(100*data(:,6),3); % wing beat frequency [Hz]
    wing.Time       = t_p; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = 20; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
  	wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
   	Wing            = Fly(wing.Pos,t_p,40,IOFreq,tt); % wing object
    
    wing.f          = interp1(wing.Time,wing.f,tt);
   	wing.Left     	= interp1(wing.Time,wing.Left,tt);
   	wing.Right     	= interp1(wing.Time,wing.Right,tt);

    Wing.WBF        = wing.f;
    Wing.WBA        = [wing.Left,wing.Right,wing.Left + wing.Right];
    %-----------------------------------------------------------------------------------------------------------------------------
	% Check WBF & WBA
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    elseif any(wing.Left>10.6) || any(wing.Right>10.6)
        fprintf('WBA out of range: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    else
        pp = pp + 1; % set next index to store data
    end
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,tt); % fit panel data
 	Pat      	= Fly(pat.Pos,tt,0.4*100,IOFreq); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between head & pattern
    Err = Pat; % error object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate iput-output relationships
    err2wing = IO_Class(Err,Wing);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    for jj = 1:n.catg
        ALL{pp,jj} = I{kk,jj};
    end

	vars = {Pat,Wing,err2wing};
    for jj = 1:length(vars)
        ALL{pp,n.catg+jj} = vars{jj};
    end
    
	qq = size(TRIAL{I{kk,1}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1}}{qq+1,ww} = vars{ww};
    end
end

ALL( all(cellfun(@isempty, ALL),2), : ) = []; % get rid of emtpty rows becuase of low WBF

clear jj ii kk pp qq ww n a b spant_p t_v hAngles data head wing pat bode tt ...
    Head Pat Wing Err pat2head pat2wing err2wing head2wing vars root t_p span
disp('LOADING DONE')

%% Fly Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
FLY = cell(N{1,1},size(TRIAL{1},2));
for kk = 1:N{1,1}
    for ii = 1:size(TRIAL{kk},2)
        FLY{kk,ii} = FlyStats(TRIAL{kk}(:,ii));
    end
end
clear kk ii
%% Grand Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
GRAND = cell(1,size(FLY,2));
for ii = 1:size(FLY,2)
    GRAND{ii} = GrandStats(FLY(:,ii));
end
clear ii
%% SAVE %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end