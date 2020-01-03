function [] = HeadWing_Raw(wave)
%% HeadWing_Raw: 
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
wave = 22.5;
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

fig = figure (1) ; clf
set(fig,'Color','w','Units','inches','Position',[2 2 6 6])
movegui(fig,'center')
ax(1) = subplot(4,1,1); ylabel('Head')
ax(2) = subplot(4,1,2); ylabel('Left Wing')
ax(3) = subplot(4,1,3); ylabel('Right Wing')
ax(4) = subplot(4,1,4); ylabel('\delta WBA') ; xlabel('Time')
set(ax,'XLim',[0 10],'FontSize',12)
linkaxes(ax,'x')
linkaxes(ax(2:4),'y')
lw = 1;
for kk = 1:N.file
    % Load HEAD & DAQ data
    disp(basename{kk})
	% load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    
    % Get data
    Time = t_v;
    Head  = rad2deg(benifly.Head);
    LWing = rad2deg(benifly.LWing);
    RWing = rad2deg(benifly.RWing);
    dWBA = LWing - RWing;
    
    subplot(4,1,1) ; hold on
        plot(Time, Head,  'k', 'LineWidth', lw)
    subplot(4,1,2) ; hold on
        plot(Time, LWing, 'b', 'LineWidth', lw)
    subplot(4,1,3) ; hold on
        plot(Time, RWing, 'r', 'LineWidth', lw)
    subplot(4,1,4) ; hold on
        plot(Time, dWBA,  'g', 'LineWidth', lw)
    
%     pause
%     clf
        
end

%% Plot Benifly Output
clear;close all;clc
% root = 'C:\Users\boc5244\Documents\temp\out';
root = 'H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\22.5\Vid\tracked';

[FILE,PATH] = uigetfile({'*.csv'},'Select Benifly output .csv', root, 'MultiSelect','on');
FILE = cellstr(FILE);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 7 7];
movegui(FIG,'center')

T   = 4;
DX  = 50; 
n = length(FILE);
for kk = 1:n
    filedata = ImportBenifly(fullfile(PATH,FILE{kk}));
    tt = linspace(0,21,size(filedata,1))';

    Head(:,kk) = filedata.Head;
    LWing(:,kk) = hampel(tt,filedata.LWing,DX,T);
    RWing(:,kk) = hampel(tt,filedata.RWing,DX,T);
    Abdomen(:,kk) = filedata.Abdomen;
    WBF(:,kk) = filedata.WBF;
end

Head_mean = mean(Head,2);
LWing_mean = mean(LWing,2);
RWing_mean = mean(RWing,2);
Abdomen_mean = mean(Abdomen,2);
WBF_mean = mean(WBF,2);

for kk = 1:n
    filedata = ImportBenifly(fullfile(PATH,FILE{kk}));
    tt = linspace(0,21,size(filedata,1))';
    
    ax(1) = subplot(5,1,1) ; hold on ; title('Head')
    plot(tt, rad2deg(filedata.Head), 'LineWidth', 0.5)

    ax(2) = subplot(5,1,2) ; hold on ; title('Left')
    [filedata.LWing] = hampel(tt,filedata.LWing);
    plot(tt, rad2deg(filedata.LWing), 'LineWidth', 0.5)

    ax(3) = subplot(5,1,3) ; hold on ; title('Right')
    [filedata.RWing] = hampel(tt,filedata.RWing);
    plot(tt, rad2deg(filedata.RWing), 'LineWidth', 0.5)

    ax(4) = subplot(5,1,4) ; hold on ; title('\delta WBA')
    plot(tt, rad2deg(filedata.LWing - filedata.RWing), 'LineWidth', 0.5)

    ax(5) = subplot(5,1,5) ; hold on ; title('Aux')
    plot(tt, filedata.WBF, 'LineWidth', 0.5)
    xlabel('Time')
end
% axes(ax(1)) ; plot(tt,rad2deg(Head_mean),'k','LineWidth',3')
% axes(ax(2)) ; plot(tt,rad2deg(LWing_mean),'k','LineWidth',3')
% axes(ax(3)) ; plot(tt,rad2deg(RWing_mean),'k','LineWidth',3')
% axes(ax(4)) ; plot(tt,rad2deg(Abdomen_mean),'k','LineWidth',3')
% axes(ax(5)) ; plot(tt,WBF_mean,'k','LineWidth',3')

set(ax,'XLim',[0 tt(end)])
set(ax(1),'YLim',20*[-1 1]);
% set(ax(2:3),'YLim',[-10 70]);
set(ax(5),'YLim',[0 1]);
linkaxes(ax,'x')
linkaxes(ax(2:3),'y')
set(ax(1:end-1),'XTickLabels','')

end