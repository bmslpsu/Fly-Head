function [] = MakeData_Ramp_Wave_HeadFree_obj(match,Fc)
%% MakeData_Ramp_Wave_HeadFree_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
Fc = 30;
match = -1;
if match==1
    clss = 'CO';
elseif match==-1
    clss = 'Anti';
elseif isnan(match)
    clss = 'All';
else
    error('Invalid match condition')
end

wave = 30;
filename = ['Ramp_HeadFree_SACCD_' clss '_filt=' num2str(Fc) '_Wave_' num2str(wave)];
rootdir = ['H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;
root.vid = fullfile(root.daq,'\Vid');
root.ang = fullfile(root.daq,'\Vid\tracked');

% Select files
[files, PATH.ang] = uigetfile({'*csv', 'files'}, 'Select files',root.ang, 'MultiSelect', 'on');
FILES = cellstr(files)';

[D,I,N,U,T,~,~,basename] = GetFileData(FILES,false,'fly','trial','vel','wave');

PATH.daq = root.daq;
PATH.vid = root.vid;
clear rootdir

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
   	%-----------------------------------------------------------------------------------------------------------------------------
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
  	%-----------------------------------------------------------------------------------------------------------------------------
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
    
    head2wing = IO_Class(Head,Wing);
    %-----------------------------------------------------------------------------------------------------------------------------
  	% Check WBF & WBA
    if min(wing.f)<150 || mean(wing.f)<200 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    elseif any(wing.Left>10.6) || any(wing.Right>10.6)
        fprintf('WBA out of range: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        % continue
    end
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
  	pat.Pos     = 3.75*round((96/10)*data(:,2)); % panel position [panel index]
    Pat      	= Fly(pat.Pos,pat.Time,[],[],Head.Time); % pattern obje
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get Saccade Stats   
    [head.SACD,head.thresh,head.count,head.rate,head.SACDRmv] = SacdDetect(Head.X(:,1),Head.Time,2.5,false);
    [wing.SACD,wing.thresh,wing.count,wing.rate,wing.SACDRmv] = SacdDetect(Wing.X(:,1),Wing.Time,1.75,false);
    
    HeadRmv = Fly(head.SACDRmv,Head.Time,[],[],tt);
    WingRmv = Fly(wing.SACDRmv,Wing.Time,[],[],tt);
    
    head.match = table(head.SACD.Direction*sign(D.vel(kk)));
    mIdx = 1:head.count;
    if ~isnan(match)
        mIdx = mIdx([head.match{:,1}]==match);
    end
    head.match.Properties.VariableNames = {'Match'};
  	wing.match = table(wing.SACD.Direction*sign(D.vel(kk)));
    wing.match.Properties.VariableNames = {'Match'};
    
    if isnan(head.count)
        head.Rate = table(0);
        wing.Rate = table(0);
    else
        head.Rate = table(nan(head.count,1));
        head.Rate{1,1} = head.rate;
    	wing.Rate = table(nan(wing.count,1));
        wing.Rate{1,1} = wing.rate;
    end
    head.Rate.Properties.VariableNames = {'Rate'};
  	wing.Rate.Properties.VariableNames = {'Rate'};
    
    head.SACD = [head.SACD , head.match];
    wing.SACD = [wing.SACD , wing.match];
     
   	Dir = table(sign(D{kk,3}),'VariableNames',{'Dir'});
    I_table = [I(kk,1:3) , rowfun(@(x) abs(x), D(kk,3)), Dir];
    I_table.Properties.VariableNames{3} = 'velIdx';
    I_table.Properties.VariableNames{4} = 'speed';
    
    [Saccade,Interval,Stimulus,Error,IntError,matchFlag] = SaccdInter(Head.X(:,1),Head.Time,head.SACD, ...
                                                                    match, Stim(:,I{kk,3}), false);
    
    var1 = {Saccade.Time, Saccade.Pos,Saccade.Vel, Error.Saccade.Pos, Error.Saccade.Vel,...
                IntError.Saccade.Pos, IntError.Saccade.Vel, Stimulus.Saccade.Pos , Stimulus.Saccade.Vel};
            
    var2 = {Interval.Time, Interval.Pos, Interval.Vel, Error.Interval.Pos, Error.Interval.Vel,...
                IntError.Interval.Pos, IntError.Interval.Vel, Stimulus.Interval.Pos , Stimulus.Interval.Vel};
    
    if isnan(head.count)
        Err_table = nan(1,5);
        loop = [];
        emptyFlag = true;
    else
        Err_table = nan(head.count,5);
        loop = size(Error.Interval.Pos,2);
       	if matchFlag
            emptyFlag = true;
        else
            emptyFlag = false;
        end
    end
    
    for jj = 1:loop
        pos_err = Error.Interval.Pos(:,jj);
        pos_err = pos_err(~isnan(pos_err));
     	vel_err = Error.Interval.Vel(:,jj);
        vel_err = vel_err(~isnan(vel_err));
        
        pos_int_err = IntError.Interval.Pos(:,jj);
        pos_int_err = pos_int_err(~isnan(pos_int_err));
      	vel_int_err = IntError.Interval.Vel(:,jj);
        vel_int_err = vel_int_err(~isnan(vel_int_err));
        
        stim_pos = Stimulus.Interval.Pos(:,jj);
     	stim_pos = stim_pos(~isnan(stim_pos));
        if ~isempty(pos_err)
            Err_table(mIdx(jj),1) = nanmean(pos_err(end-5:end),1);
            Err_table(mIdx(jj),2) = nanmean(vel_err(end-5:end),1);
            Err_table(mIdx(jj),3) = pos_int_err(end);
            Err_table(mIdx(jj),4) = vel_int_err(end);
            Err_table(mIdx(jj),5) = stim_pos(end);
        end
    end
    Err_table = splitvars(table(Err_table));
    Err_table.Properties.VariableNames = {'Position_Error','Velocity_Error','Position_IntError',...
                                                'Velocity_IntError','Stimulus_Position'};
    if isnan(head.count)
        head.I_table = I_table;
    else
        head.I_table = repmat(I_table,head.count,1);
        head.Rate{1,1} = head.Rate{1,1}*(size(Error.Interval.Pos,2)/head.count);
  	end
    head.SACD = [head.SACD , head.Rate];
    
    if isnan(wing.count)
        wing.I_table = I_table;
    else
        wing.I_table = repmat(I_table,wing.count,1);
    end
    
	head.SACD = [head.I_table , [head.SACD , Err_table]];
  	wing.SACD = [wing.I_table , wing.SACD];
    
    SACD.Head = [SACD.Head ; head.SACD];
    SACD.Wing = [SACD.Wing ; wing.SACD];
    
    if ~emptyFlag
        SACD.Saccade.Head{I{kk,3},1}  = [SACD.Saccade.Head{I{kk,3},1}  ; var1];
        SACD.Interval.Head{I{kk,3},1} = [SACD.Interval.Head{I{kk,3},1} ; var2];
    end
    
 	%-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    vars = {Pat,Head,Wing,HeadRmv,WingRmv,head2wing};
	qq = size(TRIAL{I{kk,1},I{kk,3}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1},I{kk,3}}{qq+1,ww} = vars{ww};
    end

    % pause()
    close all
    % clc
end

clear jj ii kk pp qq ww n a b  t_v hAngles data head wing pat tt I_table Dir loop Saccade Interval Stimulus Error IntError...
    Head Pat Wing  vars root t_p var1 var2 Err_table pos_err vel_err pos_int_err vel_int_err stim_pos matchFlag emptyFlag

%% Normalize Head Saccades
%---------------------------------------------------------------------------------------------------------------------------------
varnames = {'Time','Position','Velocity','Position_Error','Velocity_Error',...
                'Position_IntError','Velocity_IntError','Stimulus_Position','Stimulus_Velocity'};

clear SACCADE
SACCADE.Head = cell(N{1,3},9);
SACCADE.cIdx = cell(N{1,3},1);
dR = cell(N{1,3},1);
center = 0;
dim = 1;
for jj = 1:N{1,3}
    [SACCADE.Head{jj,1},SACCADE.cIdx{jj},~,~,dR{jj}] = nancat_center(SACD.Saccade.Head{jj}(:,1), center, dim, [], []);
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
SACCADE.HeadStats = cell2table(cellfun(@(x) MatStats(x,2), table2cell(SACCADE.Head),...
                            'UniformOutput',false),'VariableNames',varnames);
%%
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
INTERVAL.HeadStats = cell2table(cellfun(@(x) MatStats(x,2), table2cell(INTERVAL.Head),...
                            'UniformOutput',false),'VariableNames',varnames);
                        
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
    'SACD','SACCADE','INTERVAL','Stim','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end