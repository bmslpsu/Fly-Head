%% Select File(s) %%
%---------------------------------------------------------------------------------------------------------------------------------
clear;close all;clc

root = 'H:\EXPERIMENTS\Experiment_ChirpLog_HeadFixed\Vid';
[D,I,N,U,T,FILES,PATH] = GetFileData(root);

WBA_root = fullfile(PATH,'WBA');

%% Load video data & run tracking software %%
%---------------------------------------------------------------------------------------------------------------------------------
prev = 0;
WBA = nan(N.file,3);
for kk = 1:N.file
    load([PATH FILES{kk}],'vidData')
    vid = squeeze(vidData);
    dim = size(vid);
    med = median(vid,3);
    
    if (I.fly(kk) - prev)~=0
        disp('Making Mask...')
        [Mask] = MakeWingMask(med);
        close all
    end
    prev = I.fly(kk);
    
    [Wing] = WingTracker_Edge(med, Mask, 0.35, true);
    
    figure (500)
    WBA_image = getframe(gcf); hold on
    title(FILES{kk},'Interpreter', 'none')
    
    save(fullfile(WBA_root,FILES{kk}),'Wing','WBA_image','Mask','med')
    pause(0.2)
    
%     beep on
%     for jj = 1:5
%         beep; pause (0.5)
%     end
%     pause
    
    WBA(kk,1) = Wing.L;
    WBA(kk,2) = Wing.R;
    
end
close all
WBA(:,3) = WBA(:,1) + WBA(:,2);
