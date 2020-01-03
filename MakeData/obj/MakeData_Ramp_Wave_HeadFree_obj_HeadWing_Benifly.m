function [] = MakeData_Ramp_Wave_HeadFree_obj_HeadWing_Benifly(Fc,wave)
%% MakeData_Ramp_Wave_HeadFree_obj_HeadWing_Benifly: 
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       -
%

Fc = 30;
wave = 22.5;
filename = ['Ramp_HeadFree_SACCD_HeadWing_filt=' num2str(Fc) '_Wave=' num2str(wave)];
rootdir = ['H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

%% Setup Directories %%
root.daq = rootdir;
root.vid = fullfile(root.daq,'\Vid');
root.ang = fullfile(root.daq,'\Vid\tracked');

% Select files
[files, PATH.ang] = uigetfile({'*csv', 'files'}, 'Select files',root.ang, 'MultiSelect', 'on');
FILES = cellstr(files)';

[D,I,N,U,T,~,~,basename] = GetFileData(FILES,false,'fly','trial','vel','wave');

PATH.daq = root.daq;
PATH.vid = root.vid;
PATH.sacd = fullfile(root.ang,'SACD');

%% Get Data %%
disp('Loading...')
clear SACD
badtrial = {};

Vel = U.vel{1};
tt = (0:(1/200):10)'; 
Stim = (Vel*tt')';
bad = 1;
T = 4;
DX = 50;
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
% 	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    disp(basename{kk})
   	
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
  	
    % Get wing data
    wing.Time       = t_v; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = Fc; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = rad2deg(hampel(wing.Time,benifly.LWing)); % left wing
    wing.Right      = rad2deg(hampel(wing.Time,benifly.RWing)); % right wing
 	wing.Left       = filtfilt(b,a,wing.Left); % left wing
    wing.Right      = filtfilt(b,a,wing.Right); % right wing
    
    wing.WBA        = wing.Left - wing.Right; % dWBA (L-R)
    Wing            = Fly(wing.WBA,wing.Time,wing.Fc,[],tt); % head object
	    
    [SACD.Head,~,~,~,~,~,~] = Sacd_Manual(Head.X(:,1),Head.Time,false);
    pause
 	[SACD.Wing,~,~,~,~,~,~] = Sacd_Manual(Wing.X(:,1),Wing.Time,false);
    pause
    close all
    
    % save(fullfile(PATH.sacd, [basename{kk} '.mat']), 'SACD')
end
disp('DONE')

%%

fig = figure (2) ; clf
set(fig,'Color','k','Units','Inches','Position',[2 2 8 3])
movegui(fig,'center')
ax(1) = subplot(1,1,1); hold on
set(ax,'Color','k','XColor','w','YColor','w','FontSize',12,'YLim',22*[-1 1])
h(1) = plot(Head.Time,Head.X(:,1) - mean(Head.X(:,1)),'c','LineWidth',1);
h(2) = plot(Head.Time,Wing.X(:,1) - mean(Wing.X(:,1)),'r','LineWidth',1);
xlabel('Time (s)')
ylabel('Angle (°)')
leg = legend(h,'Head','\Delta WBA');
leg.Box = 'off';
leg.TextColor = 'w';


end