function [] = Test_IOFreq()
%% Test_IOFreq:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%

root = 'H:\DATA\Rigid_Data\';

[FILES,~] = uigetfile({'*.mat', 'DAQ-files'},...
    'Select head angle trials', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nAmp = length(FILES);
Amp = nan(nAmp,1);
for ww = 1:nAmp
filedata = textscan(FILES{ww}, '%s', 'delimiter', '_');
    Amp(ww) = str2double(filedata{1}{3});
end

HeadFree = cell(nAmp,1);
for ww = 1:nAmp
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N');
end

%% Complex Gain Calculations
clearvars -except nAmp Amp HeadFree

filename = 'Sine_HeadFree_ComplexGain';

cIdx = 1;
xIdx = 1;

Freq = HeadFree{1}.U{1,3}{1}';
nFreq = HeadFree{1}.N{1,3};

figure (1) ; clf ; hold on
color_freq = prism(nFreq);
Test = cell(nAmp,1);
for ww = 1:nAmp
    for jj = 1:nFreq
        pp = 1;
        for kk = 1:HeadFree{ww}.N{1,1}
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1)
                Fv = HeadFree{ww}.TRIAL{kk,jj}{ii,cIdx}.Fv;
                Mag = HeadFree{ww}.TRIAL{kk,jj}{ii,cIdx}.Mag(:,xIdx);
                Phase = HeadFree{ww}.TRIAL{kk,jj}{ii,cIdx}.Phase(:,xIdx);
                
                IOFreq  = HeadFree{ww}.TRIAL{kk,jj}{ii,cIdx}.IOFreq;
                IOMag   = HeadFree{ww}.TRIAL{kk,jj}{ii,cIdx}.IOMag(:,xIdx);
                IOPhase = HeadFree{ww}.TRIAL{kk,jj}{ii,cIdx}.IOPhase(:,xIdx);
                
                subplot(2,1,1) ; hold on
                cla
                plot(Fv,Mag,'Color',rand(1,3),'LineWidth',1)
                plot(IOFreq,IOMag,'r*','MarkerSize',10)
                xlim([0 13])   
                ylim([0 20])
                
                subplot(2,1,2) ; hold on
                cla
               	plot(Fv,Phase,'Color',rand(1,3),'LineWidth',1)
                plot(IOFreq,IOPhase,'r*','MarkerSize',10)
                xlim([0 13])
                ylim([-2*pi 2*pi])
                
                pause()
                                
                % Test{ww}(pp,jj)     = R;
                
                pp = pp + 1;
            end
        end
    end

end
clear ww jj kk ii pp R I



end