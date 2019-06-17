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

[D,I,N,U,T] = GetFileData(FILES,false);
[D,I,N,U,~] = GetFileData(FILES,false);
% clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
TRIAL = cell(N{1,1},N{1,3});

clear SACD
SACD.Head = [];
SACD.Wing = [];
SACD.Saccade.Head = cell(N{1,3},1);
SACD.Interval.Head = cell(N{1,3},1);
SACD.Stimulus.Saccade.Head = cell(N{1,3},1);
SACD.Stimulus.Interval.Head = cell(N{1,3},1);

Vel = 3.75*U{1,3}{1};
tt = (0:(1/200):9.8)';
Stim = (Vel*tt')';
pp = 0;
for kk = 1:N{1,end}
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.ang, FILES{kk}),'hAngles','t_v'); % load head angles % time arrays
    %-----------------------------------------------------------------------------------------------------------------------------
    % Check WBF
	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        warning('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    else
        pp = pp + 1; % set next index to store data
    end
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get head data
    head.Time = t_v;
    head.Pos = hAngles;
    Head = Fly(head.Pos,head.Time,40,[],tt); % head object
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
   	Wing = Fly(wing.Pos,t_p,40,[],tt); % wing object
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,tt); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,[],[]); % pattern object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    vars = {Pat,Head,Wing};
	qq = size(TRIAL{I{kk,1},I{kk,3}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1},I{kk,3}}{qq+1,ww} = vars{ww};
    end
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get Saccade Stats
    [head.SACCD,head.thresh,head.count,head.rate] = GetSaccade(Head,300,false);
    [wing.SACCD,wing.thresh,wing.count,wing.rate] = GetSaccade(Wing,1.75,false);
    
    head.match = table(head.SACCD.Direction*sign(D{kk,end}));
    head.match.Properties.VariableNames = {'Match'};
  	wing.match = table(wing.SACCD.Direction*sign(D{kk,end}));
    wing.match.Properties.VariableNames = {'Match'};
    
    head.SACCD = [head.SACCD , head.match];
    wing.SACCD = [wing.SACCD , wing.match];
    
   	Dir = table(sign(D{kk,3}),'VariableNames',{'Dir'});
    I_table = [I(kk,1:2) , rowfun(@(x) abs(x), D(kk,3)), Dir];
    
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
    
    SACD.Head = [SACD.Head ; head.SACCD];
    SACD.Wing = [SACD.Wing ; wing.SACCD];
    
    [Saccade,Interval,Stimulus,Error,IntError] = SaccdInter(Head.X(:,1),Head.Time(:,1),head.SACCD,nan,Stim(:,I{kk,3}),false);
    
    var1 = {Saccade.Time, Saccade.Pos,Saccade.Vel, Error.Saccade.Pos, Error.Saccade.Vel,...
                IntError.Saccade.Pos, IntError.Saccade.Vel, Stimulus.Saccade.Pos , Stimulus.Saccade.Vel};
            
	var2 = {Interval.Time, Interval.Pos, Interval.Vel, Error.Interval.Pos, Error.Interval.Vel,...
                IntError.Interval.Pos, IntError.Interval.Vel, Stimulus.Interval.Pos , Stimulus.Interval.Vel};
                
    SACD.Saccade.Head{I{kk,3},1}  = [SACD.Saccade.Head{I{kk,3},1}  ; var1];
    SACD.Interval.Head{I{kk,3},1} = [SACD.Interval.Head{I{kk,3},1} ; var2];
    close all
end

clear jj ii kk pp qq ww n a b  t_v hAngles data head wing pat tt ...
    Head Pat Wing  vars root t_p 
disp('LOADING DONE')

%% Transfrom data to arrays
%---------------------------------------------------------------------------------------------------------------------------------
varnames = {'Time','Position','Velocity','Position_Error','Velocity_Error',...
                'Position_IntError','Velocity_IntError','Stimulus_Position','Stimulus_Velocity'};

clear SACCADE
SACCADE.Head = cell(N{1,3},9);
dR = cell(N{1,3},1);
center = 0;
dim = 1;
for jj = 1:N{1,3}
    [SACCADE.Head{jj,1},~,~,~,dR{jj}] = nancat_center(SACD.Saccade.Head{jj}(:,1), center, dim);
    for ww = 2:size(SACD.Saccade.Head{jj},2)
        for kk = 1:size(SACD.Saccade.Head{jj},1)
            for ii = 1:size(SACD.Interval.Head{jj}{kk,ww},2)
                SACCADE.Head{jj,ww}{kk,1}(:,ii) = cat_pad(SACD.Saccade.Head{jj}{kk,ww}(:,ii), dR{jj}{kk}(:,ii),nan);
            end
        end
    	SACCADE.Head{jj,ww} = cat(2,SACCADE.Head{jj,ww}{:});
    end
end
SACCADE.Head = cell2table(SACCADE.Head,'VariableNames',varnames);

clear INTERVAL
INTERVAL.Head = cell(N{1,3},9);
dR = cell(N{1,3},1);
center = 0;
dim = 1;
for jj = 1:N{1,3}
    [INTERVAL.Head{jj,1},~,~,~,dR{jj}] = nancat_center(SACD.Interval.Head{jj}(:,1), center, dim);
    for ww = 2:size(SACD.Interval.Head{jj},2)
        for kk = 1:size(SACD.Interval.Head{jj},1)
            for ii = 1:size(SACD.Interval.Head{jj}{kk,ww},2)
                INTERVAL.Head{jj,ww}{kk,1}(:,ii) = cat_pad(SACD.Interval.Head{jj}{kk,ww}(:,ii), dR{jj}{kk}(:,ii),nan);
            end
        end
    	INTERVAL.Head{jj,ww} = cat(2,INTERVAL.Head{jj,ww}{:});
    end
end
INTERVAL.Head = cell2table(INTERVAL.Head,'VariableNames',varnames);

INTERVAL.Head



%%
TIME = cell(N{1,3},1);
POS  = cell(N{1,3},1);
dR   = cell(N{1,3},1);
for jj = 1:N{1,3}
    [TIME{jj},~,~,~,dR{jj}] = nancat_center(SACD.Interval.Head{jj}(:,1),0,1);
    for kk = 1:size(SACD.Interval.Head{jj},1)
        for ii = 1:size(SACD.Interval.Head{jj}{kk,2},2)
            POS{jj}{kk,1}(:,ii) = cat_pad(SACD.Interval.Head{jj}{kk,2}(:,ii), dR{jj}{kk}(:,1),nan);
        end
    end
	POS{jj} = cat(2,POS{jj}{:});
end

figure (1) ; clf
for jj = 1:N{1,3}
    subplot(2,3,jj) ; hold on
    plot(TIME{jj},POS{jj})
%     if jj<=3
%         ylim([-100 1000])
%     else
%         ylim([-1000 100])
%     end
%     xlim(0.06*[-1 1])

    xlim([0 2])
end

CC = repmat({'r','g','b'},1,2);
FIG = figure (2) ; clf
FIG.Color = 'k';
ax = gca;
ax.Color = 'k';
set(ax,'YColor','w','XColor','w')
for jj = [1 4 2 5 3 6]
	hold on
    plot(TIME{jj},POS{jj},'Color',CC{jj})
    xlim([0 2])
end

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
    'SACCADE','INTERVAL','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end