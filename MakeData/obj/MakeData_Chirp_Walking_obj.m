function [] = MakeData_Chirp_Walking_obj(rootdir)
%% MakeData_Chirp_Walking_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root        : root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = 'H:\EXPERIMENTS\Experiment_Chirp_Walking\mat';
filename = 'Chirp_Walking_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------

% Select files
[FILES, PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', rootdir, 'MultiSelect','on');
FILES = cellstr(FILES)';

[D,I,N,U,T] = GetFileData(FILES,true);

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
IOFreq = 1;
disp('Loading...')
ALL 	= cell([N{1,end},6]); % cell array to store all data objects
TRIAL  	= cell(N{1,1},N{1,3});
n.catg  = size(N,2) - 1;
tt = linspace(1.2,21.2,100*20)';
for kk = 1:N{1,end}
    disp(kk)
    % Load Kinefly Data
	load(fullfile(PATH, FILES{kk}),'FlyState','AI'); % load pattern x-position
    syncTime = FlyState{1,1};
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get head data
    head.Time = FlyState{:,1} - syncTime;
    head.Pos = rad2deg(FlyState{:,2} - mean(FlyState{:,2}));
    head.Pos = hampel(head.Time,head.Pos,10);
    Head = Fly(head.Pos,head.Time,30,IOFreq,tt); % head object
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= AI{:,1} - syncTime;
    pat.Pos 	= panel2deg( AI{:,2} );  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,tt,false,false); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,[],IOFreq); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between head & pattern
    head.Err = Pat.X(:,1) - Head.X(:,1); % calculate position error between head & pattern [deg]
    Err = Fly(head.Err,Head.Time,0.4*Head.Fs,IOFreq); % error object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate iput-output relationships
    pat2head    = IO_Class(Pat,Head);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    for jj = 1:n.catg
        ALL{kk,jj} = I{kk,jj};
    end
    
	vars = {Pat,Head,Err,pat2head};
    for jj = 1:length(vars)
        ALL{kk,n.catg+jj} = vars{jj};
    end
    
	qq = size(TRIAL{I{kk,1},1},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1},1}{qq+1,ww} = vars{ww};
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