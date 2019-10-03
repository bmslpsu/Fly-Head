function [] = MakeData_SOS_Walking_obj(rootdir)
%% MakeData_SOS_HeadFree_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       rootdir    	:   root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = 'E:\Walking_Experiments\SOS\mat';
filename = 'SOS_Walking_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;
% root.ang = fullfile(root.daq,'\Vid\Angles\');
root.ang = fullfile(root.daq,'\Retrack\');

% Select files
[FILES, PATH.ang] = uigetfile({'*.csv', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = cellstr(FILES)';

PATH.daq = root.daq;

[D,I,N,U,T] = GetFileData(FILES);

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------

IOFreq = [1, 3.1, 5.3, 7.4, 9.6];
disp('Loading...')
ALL 	= cell([N{1,end},10]); % cell array to store all data objects
TRIAL  	= cell(N{1,1},1);
n.catg  = size(N,2) - 1;
pp = 0;
span = 1:2000;
tt = linspace(0,20,2000)';
for kk = 1:N{1,end}
    disp(kk)
    % Load HEAD & DAQ data
    data = [];
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.ang, FILES{kk}),'hAngles','t_v'); % load head angles % time arrays
	%-----------------------------------------------------------------------------------------------------------------------------
    % Get head data
    head.Time = t_v(span);
    head.Pos = hAngles(span) - mean(hAngles(span));
    Head = Fly(head.Pos,head.Time,40,IOFreq); % head object
    pp = pp + 1; % set next index to store data
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,tt); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,0.4*Head.Fs,IOFreq); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between head & pattern
    head.Err = Pat.X(:,1) - Head.X(:,1); % calculate position error between head & pattern [deg]
    Err = Fly(head.Err,Head.Time,0.4*Head.Fs,IOFreq); % error object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate iput-output relationships
    pat2head    = IO_Class(Pat,Head);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
	vars = {Pat,Head,Err,pat2head};
	qq = size(TRIAL{I{kk,1}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1}}{qq+1,ww} = vars{ww};
    end
end

clear jj ii kk pp qq ww n a b spant_p t_v hAngles data head pat bode tt ...
    Head Pat Err pat2head vars root t_p span
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

%%
save(fullfile(targetDir,fname), 'func');

%% SAVE %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end