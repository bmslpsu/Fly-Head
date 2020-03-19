function [] = MakeData_Chirp_HeadFree_DFT(rootdir)
%% MakeData_Chirp_HeadFree_DFT:
%   INPUTS:
%       root        : root directory
%   OUTPUTS:
%       -
%
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_ChirpLog_HeadFree';
filename = 'Chirp_HeadFree_DATA';

%% Setup Directories %%
root.daq = rootdir;
root.ang = fullfile(root.daq,'\Vid\Angles\');

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = cellstr(FILES)';

PATH.daq = root.daq;

[D,I,N,U,T] = GetFileData(FILES,true);

clear rootdir
%% Get Data %%
IOFreq = [];
disp('Loading...')
TRIAL = cell(N.fly,N.Amp);
pp = 0;
tintrp = (0:(1/200):20)';
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.ang, FILES{kk}),'hAngles','t_v'); % load head angles % time arrays
    
    % Get head data
    head.Time = t_v;
    head.Pos = hAngles - mean(hAngles);
    Head = Fly(head.Pos,head.Time,40,IOFreq,tintrp); % head object
  	
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
   	Wing            = Fly(wing.Pos,t_p,40,IOFreq,tintrp); % wing object
    
    wing.f          = interp1(wing.Time,wing.f,Head.Time);
   	wing.Left     	= interp1(wing.Time,wing.Left,Head.Time);
   	wing.Right     	= interp1(wing.Time,wing.Right,Head.Time);

    Wing.WBF        = wing.f;
    Wing.WBA        = [wing.Left,wing.Right,wing.Left + wing.Right];
	
    % Check WBF & WBA
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    elseif any(wing.Left>10.6) || any(wing.Right>10.6)
        fprintf('WBA out of range: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        % continue
    else
        pp = pp + 1; % set next index to store data
    end
	
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2)); % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,Head.Time); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,[],IOFreq,tintrp); % pattern object 
	
 	% Calculate error between head & pattern
    head.Err = Pat.X(:,1) - Head.X(:,1); % calculate position error between head & pattern [deg]
    Err = Fly(head.Err,Head.Time,0.4*Head.Fs,IOFreq); % error object
    
    freq_bins = ((0:0.2:100)-0.0001)';
    [binFv,pat_Mag,NN,bins] = DFT_bin(Pat.Fv,Pat.Mag(:,1),freq_bins);    
    [~,head_Mag,~,~] = DFT_bin(Pat.Fv,Head.Mag(:,1),freq_bins);
    [~,wing_Mag,~,~] = DFT_bin(Pat.Fv,Wing.Mag(:,1),freq_bins);
    [~,err_Mag,~,~]  = DFT_bin(Pat.Fv,Err.Mag(:,1),freq_bins);
    
    [~,pat_phase,~,~]  = DFT_bin(Pat.Fv,Pat.Phase(:,1),freq_bins);   
 	[~,head_phase,~,~] = DFT_bin(Pat.Fv,Head.Phase(:,1),freq_bins);
    [~,wing_phase,~,~] = DFT_bin(Pat.Fv,Wing.Phase(:,1),freq_bins);
    [~,err_phase,~,~]  = DFT_bin(Pat.Fv,Err.Phase(:,1),freq_bins);
    
	[~,pat_FREQ,~,~]  = DFT_bin(Pat.Fv,Pat.FREQ(:,1),freq_bins);   
 	[~,head_FREQ,~,~] = DFT_bin(Pat.Fv,Head.FREQ(:,1),freq_bins);
    [~,wing_FREQ,~,~] = DFT_bin(Pat.Fv,Wing.FREQ(:,1),freq_bins);
    [~,err_FREQ,~,~]  = DFT_bin(Pat.Fv,Err.FREQ(:,1),freq_bins);
    
    pat2head_gain  = head_Mag ./ pat_Mag;
    pat2wing_gain  = wing_Mag ./ pat_Mag;
    err2wing_gain  = wing_Mag ./ err_Mag;
    wing2head_gain = wing_Mag ./ head_Mag;
    
 	pat2head_phase  = -(pat_phase  - head_phase);
    pat2wing_phase  = -(pat_phase  - wing_phase);
    err2wing_phase  = -(err_phase  - wing_phase);
    wing2head_phase = -(head_phase - wing_phase);
    
    pat2head  = head_FREQ ./ pat_FREQ;
    pat2wing  = wing_FREQ ./ pat_FREQ;
    err2wing  = wing_FREQ ./ err_FREQ;
    wing2head = wing_FREQ ./ head_FREQ;
    
    % Calculate iput-output relationships
    vars = {binFv,pat_Mag,head_Mag,wing_Mag,err_Mag, ...
            pat2head,pat2wing,err2wing,wing2head,...
            pat2head_gain,err2wing_gain,wing2head_gain,pat2wing_gain, ...
            pat2head_phase,pat2wing_phase,err2wing_phase,wing2head_phase};
        
	qq = size(TRIAL{I{kk,1},I{kk,3}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1},I{kk,3}}{qq+1,ww} = vars{ww};
    end
end

 %% Means
n_vars = length(vars);
FLY = cell(N.fly,N.Amp);
FLY_mean = cell(N.fly,N.Amp);
for jj = 1:N.Amp
    for kk = 1:N.fly
        for ii = 1:n_vars
            FLY{kk,jj}{ii} = cat(2,TRIAL{kk,jj}{:,ii});
            if ii < 14
                FLY_mean{kk,jj}{1,ii} = mean(FLY{kk,jj}{ii},2);
            else
                FLY_mean{kk,jj}{1,ii} = circ_mean(FLY{kk,jj}{ii},[],2);
            end
        end
    end
end

GRAND_all = cell(N.Amp,n_vars);
GRAND = cell(N.Amp,n_vars);
for jj = 1:N.Amp
    for ii = 1:n_vars
        for kk = 1:N.fly
            GRAND_all{jj,ii} = [GRAND_all{jj,ii}  , FLY{kk,jj}{ii}];
            GRAND{jj,ii} = [GRAND{jj,ii}  , FLY_mean{kk,jj}{ii}];  
        end
    end
end
GRAND_mean_mag = cellfun(@(x) mean(x,2), GRAND, 'UniformOutput', false);
GRAND_mean_circ = cellfun(@(x) circ_mean(x,[],2), GRAND, 'UniformOutput', false);

%% HEAD: All Amplitudes
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 N.Amp*2.5 5];
FIG.Name = filename;
movegui(FIG,'center')
CC = hsv(N.Amp+1);
ax = gobjects(3,N.Amp);
clc
p2h = 10;
phs_shift = 4;
pp = 1;
Fv = GRAND_mean_mag{jj,1};
for jj = 1:N.Amp % amplitudes    
    ax(1,jj) = subplot(2,N.Amp,pp);
        hold on
        ylabel('Gain (°/°)')
        ax(1,jj).YTick = unique(sort([ax(1,jj).YTick ax(1,jj).YLim(2)]));
      	ax(1,jj).XTick = unique(sort([min(ax(1,jj).XLim) 2:2:12]));
        
        GAIN_ALL = abs(GRAND_all{jj,6});
        % GAIN_ALL = GRAND_all{jj,p2h};
        
        plot(Fv, GAIN_ALL, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.5)
        
        GAIN = abs(GRAND_mean_mag{jj,6});
        % GAIN = GRAND_mean_mag{jj,p2h};
        
        h.patch = PlotPatch(GAIN, 0*GAIN, Fv, 1, N.fly, CC(jj+1,:), [0.4 0.4 0.6], 0.5, 1);
        vel = round(U.Amp{1}(jj)*2*pi*ax(1,jj).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2,jj) = subplot(2,N.Amp,pp + N.Amp);
        hold on
        title([num2str(U.Amp{1}(jj)) , '°']);
        ylabel('Phase Difference (°)')
        xlabel('Frequency (Hz)');
        ax(2,jj).YTick = -180:60:180;
        
        PHASE_ALL = rad2deg(angle((GRAND_all{jj,6})));
        % PHASE_ALL = rad2deg(GRAND_all{jj,p2h+phs_shift});
        
     	plot(Fv, PHASE_ALL, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.5)
        
        PHASE = rad2deg(angle(GRAND_mean_mag{jj,6}));
        % PHASE = rad2deg(GRAND_mean_circ{jj,p2h+phs_shift});
        
        h.patch = PlotPatch(PHASE, 0*PHASE, Fv, 1, N.fly ,CC(jj+1,:), [0.4 0.4 0.6],0.5,1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
                
	ax(3,jj) = axes;
        ax(3,jj).Position = ax(1,jj).Position + [0 -0.00 0 0];
        ax(3,jj).Color = 'none';
        ax(3,jj).YAxisLocation = 'right';
        ax(3,jj).YAxis.Color = 'none';
        ax(3,jj).XAxisLocation = 'top';
        ax(3,jj).XLim = ax(1,jj).XLim;
        ax(3,jj).XTick = ax(1,jj).XTick;
        ax(3,jj).XTickLabels = velLabel;
      	xlabel('Peak Velocity (°/s)');
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,jj).XLim(1),2:2:12]))
linkaxes(ax,'x')
linkaxes(ax(1,:),'y')
linkaxes(ax(2:3,:),'y')
set(ax(1,:),'YLim',[0 1.2])
set(ax(2:3,:),'YLim',rad2deg(pi*[-1 1]))
%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end