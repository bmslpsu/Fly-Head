function [] = MakeData_Chirp_HeadFree_TEST(rootdir,filename)
%% MakeData_Chirp_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       PAT     : pattern structure
%       WING   	: wings structure
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = 'H:\EXPERIMENTS\Experiment_ChirpLog_HeadFree\';
% filename = 'Chirp_HeadFree_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;
root.ang = [root.daq '\Vid\Angles\'];

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = cellstr(FILES)';

PATH.daq = root.daq;

[~,I,N,U] = GetFileData(FILES);

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
ALL 	= [table2cell(I), cell([N{1,end},4])]; % cell array to store all data objects
TRIAL  	= cell(N{1,1},1);
for kk = 1:N{1,1}
	TRIAL{kk} = cell(N{1,3},1);
end
n.catg = size(N,2) - 1;
pp = 0;
for kk = 1:N{1,4}
    disp(kk)
    % Load HEAD & DAQ data
    data = [];
	load([PATH.daq   FILES{kk}],'data','t_p'); % load pattern x-position
    load([PATH.ang   FILES{kk}],'hAngles','t_v'); % load head angles % time arrays
    %-----------------------------------------------------------------------------------------------------------------------------
    % Check WBF
	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',FD.Fly(kk),FD.Trial(kk))
        continue
    else
        pp = pp + 1; % set next index to store data
    end
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get head data
    head.Time = t_v;
    head.Pos = hAngles - mean(hAngles);
    Head = Fly(head.Pos,head.Time,20); % head object
  	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ
    wing.Time       = t_p; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = 20; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
  	wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
   	Wing = Fly(wing.Pos,t_p,20,Head.Time); % wing object
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,Head.Time); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,0.4*Head.Fs); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between head & pattern
    head.Err = Pat.X(:,1) - Head.X(:,1); % calculate position error between head & pattern [deg]
    Err = Fly(head.Err,Head.Time,0.4*Head.Fs); % error object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate iput-output relationships
    pat2head    = IO_Class(Pat,Head);
    err2wing    = IO_Class(Err,Head);
	head2wing   = IO_Class(Head,Wing);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    for jj = 1:n.catg
        ALL{pp,jj} = I{kk,jj};
    end
    ALL{pp,n.catg+1}	= Pat;
	ALL{pp,n.catg+2} 	= Head;
    ALL{pp,n.catg+3}	= Wing;
	ALL{pp,n.catg+4} 	= Err;
    ALL{pp,n.catg+5}	= pat2head;
	ALL{pp,n.catg+6} 	= err2wing;
    ALL{pp,n.catg+7}	= head2wing;
    
    vars = {Pat,Head,Wing,Err,pat2head,err2wing,head2wing};
	qq = size(TRIAL{I{kk,1}}{I{kk,3}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1}}{I{kk,3}}{qq+1,ww} = vars{ww};
    end
end
clear jj ii kk pp qq ww n a b t_p t_v hAngles data head wing pat bode tt Head Pat Wing Err pat2head err2wing head2wing vars root
disp('LOADING DONE')

%% Fly Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
FLY = cell(N{1,1},1);
for kk = 1:N{1,1}
    FLY{kk} = cell(N{1,3},1);
    for jj = 1:N{1,3}
        FLY{kk}{jj} = FlyStats(TRIAL{kk}{jj});
    end
end

FlyGroup = cell(N{1,3},1);
for jj = 1:N{1,3} 
    FlyGroup{jj} = cell(N{1,1},1);
    for kk = 1:N{1,1}
        FlyGroup{jj}{kk} = FLY{kk}{jj};
    end 
end
clear jj ii
%% Grand Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
GRAND = cell(N{1,3},1);
for jj = 1:N{1,3}
    GRAND{jj} = GrandStats(FlyGroup{jj});
end
clear jj
end