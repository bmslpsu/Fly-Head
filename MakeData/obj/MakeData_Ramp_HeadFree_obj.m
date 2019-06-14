function [] = MakeData_Ramp_HeadFree_obj(rootdir,Amp)
%% MakeData_Ramp_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% rootdir = 'F:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast';
% Amp = 60;
% filename = ['Ramp_HeadFree_' num2str(Amp) '_DATA'];
filename = 'Ramp_HeadFree_SACCD';
rootdir = 'H:\EXPERIMENTS\Experiment_Ramp';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
% root.daq = fullfile(rootdir,num2str(Amp));
root.daq = rootdir;
root.ang = fullfile(root.daq,'\Vid\Angles\');

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = cellstr(FILES)';

PATH.daq = root.daq;

[D,I,N,U,T] = GetFileData(FILES,true);
[D_dir] = GetFileData(FILES,false);
clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
TRIAL = cell(N{1,1},N{1,3});
SACCD.Head = [];
SACCD.Wing = [];
pp = 0;
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
    head.Time = t_v;
    head.Pos = hAngles;
    Head = Fly(head.Pos,head.Time,40,[]); % head object
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
   	Wing = Fly(wing.Pos,t_p,40,[],Head.Time); % wing object
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,Head.Time); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,[],[]); % pattern object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    vars = {Pat,Head,Wing};
	qq = size(TRIAL{I{kk,1},I{kk,3}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1},I{kk,3}}{qq+1,ww} = vars{ww};
    end
    
    [head.SACCD,head.thresh,head.count,head.rate] = GetSaccade(Head,2.5,true);
    [wing.SACCD,wing.thresh,wing.count,wing.rate] = GetSaccade(Wing,1.75,true);
    
    head.match = table(head.SACCD.Direction*sign(D_dir{kk,end}));
    head.match.Properties.VariableNames = {'Match'};
  	wing.match = table(wing.SACCD.Direction*sign(D_dir{kk,end}));
    wing.match.Properties.VariableNames = {'Match'};
    
    head.SACCD = [head.SACCD , head.match];
    wing.SACCD = [wing.SACCD , wing.match];
    
   	Dir = table(sign(D{kk,3}),'VariableNames',{'Dir'});
    I_table = [I(kk,1:2) , D(kk,3) , Dir];
    
    if isnan(head.count)
        head.I_table = I_table;
    else
        head.I_table = repmat(I_table,head.count,1);
    end
    
    if isnan(wing.count)
        wing.I_table = I_table;
    else
        wing.I_table = repmat(I_table,wing.count,1);
    end
    
    head.SACCD = [head.I_table , head.SACCD];
  	wing.SACCD = [wing.I_table , wing.SACCD];
    
    SACCD.Head = [SACCD.Head ; head.SACCD];
    SACCD.Wing = [SACCD.Wing ; wing.SACCD];
    
end

% SACCD.Head( table2array(varfun(@isnan, SACCD.Head)), : ) = [];
% SACCD.Wing( table2array(varfun(@isnan, SACCD.Wing)), : ) = [];

clear jj ii kk pp qq ww n a b  t_v hAngles data head wing pat tt ...
    Head Pat Wing  vars root t_p 
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
    'SACCD','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end