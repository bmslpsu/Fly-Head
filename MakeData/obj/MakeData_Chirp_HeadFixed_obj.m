function [] = MakeData_Chirp_HeadFixed_obj(rootdir)
%% MakeData_Chirp_HeadFixed_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root        : root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = 'F:\EXPERIMENTS\Experiment_ChirpLog_HeadFixed';
filename = 'Chirp_HeadFixed_DATA';
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
IOFreq = 1;
disp('Loading...')
ALL 	= cell([N{1,end},3]); % cell array to store all data objects
TRIAL  	= cell(N{1,1},N{1,3});
n.catg  = size(N,2) - 1;
pp = 0;
Head.Time = (0:(1/200):20)';
Head.Fs = 1/mean(diff(Head.Time));
for kk = 1:N{1,end}
    disp(kk)
    % Load DAQ data
    data = [];
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
    %-----------------------------------------------------------------------------------------------------------------------------
    % Check WBF
	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
%         continue
    else
        pp = pp + 1; % set next index to store data
    end
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
   	Wing = Fly(wing.Pos,t_p,40,IOFreq,Head.Time); % wing object
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,Head.Time); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,0.4*Head.Fs,IOFreq); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between head & pattern
    Err = Pat;
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate iput-output relationships
    err2wing = IO_Class(Err,Wing);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    for jj = 1:n.catg
        ALL{pp,jj} = I{kk,jj};
    end
    ALL{pp,n.catg+1}	= Pat;
    ALL{pp,n.catg+2}	= Wing;
    ALL{pp,n.catg+3} 	= err2wing;
    
    vars = {Pat,Wing,err2wing};
	qq = size(TRIAL{I{kk,1},I{kk,3}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1},I{kk,3}}{qq+1,ww} = vars{ww};
    end
end

ALL( all(cellfun(@isempty, ALL),2), : ) = []; % get rid of emtpty rows becuase of low WBF

clear jj ii kk pp qq ww n a b spant_p t_v hAngles data head wing pat bode tt ...
    Head Pat Wing Err pat2head err2wing head2wing vars root t_p span
disp('LOADING DONE')

%% Fly Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
FLY = cell(N{1,3},1);
for kk = 1:N{1,1}
    for jj = 1:N{1,3}
        for ii = 1:size(TRIAL{kk,jj},2)
            FLY{jj}{kk,ii} = FlyStats(TRIAL{kk,jj}(:,ii));
        end
    end
end
clear kk jj ii
%% Grand Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
GRAND = cell(size(FLY,2),size(FLY{1},2));
for jj = 1:N{1,3}
    for ii = 1:size(FLY{jj},2)
        GRAND{jj,ii} = GrandStats(FLY{jj}(:,ii));
    end
end
clear jj ii
%% SAVE %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save(['F:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'ALL','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end