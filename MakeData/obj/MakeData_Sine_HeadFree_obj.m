function [] = MakeData_Sine_HeadFree_obj(rootdir,Amp)
%% MakeData_Sine_HeadFree_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root        : root directory
%       Amp         : amplitude folder
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% rootdir = 'F:\EXPERIMENTS\Experiment_Sinusoid';
filename = ['Sine_HeadFree_' num2str(Amp) '_DATA'];
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = fullfile(rootdir, num2str(Amp));
root.ang = fullfile(root.daq, '\Vid\Angles\');

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = cellstr(FILES)';

PATH.daq = root.daq;

[D,I,N,U,T] = GetFileData(FILES,true);

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
IOFreq = U{1,3}{1};
disp('Loading...')
ALL 	= cell([N{1,end},6]); % cell array to store all data objects
TRIAL  	= cell(N{1,1},N{1,3});
n.catg  = size(N,2) - 1;
pp = 0;
span = 1:2000;
tt = linspace(0,10,2000)';
for kk = 1:N{1,end}
    disp(kk)
    % Load HEAD & DAQ data
    data = [];
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.ang, FILES{kk}),'hAngles','t_v'); % load head angles % time arrays
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
    % Get head data
    head.Time = t_v(span);
    head.Pos = hAngles - mean(hAngles);
    Head = Fly(head.Pos(span),head.Time,40,IOFreq(I{kk,3})); % head object
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
   	Wing = Fly(wing.Pos,t_p,40,IOFreq(I{kk,3}),tt); % wing object
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,tt,false,false); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,0.4*Head.Fs,IOFreq(I{kk,3})); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between head & pattern
    head.Err = Pat.X(:,1) - Head.X(:,1); % calculate position error between head & pattern [deg]
    Err = Fly(head.Err,Head.Time,0.4*Head.Fs,IOFreq(I{kk,3})); % error object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate iput-output relationships
    pat2head    = IO_Class(Pat,Head);
    err2wing    = IO_Class(Err,Wing);
	head2wing   = IO_Class(Head,Wing);
    pat2wing    = IO_Class(Pat,Wing);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    for jj = 1:n.catg
        ALL{pp,jj} = I{kk,jj};
    end
    
	vars = {Pat,Head,Wing,Err,pat2head,err2wing,head2wing,pat2wing};
    for jj = 1:length(vars)
        ALL{pp,n.catg+jj} = vars{jj};
    end
    
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
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'ALL','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end