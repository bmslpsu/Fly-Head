function [FIG] = MakeFig_WBA()
%% MakeFig_WBA:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\EXPERIMENTS\Experiment_SOS_v2\WBA';
% root = 'H:\EXPERIMENTS\Experiment_ChirpLog_HeadFree\Vid\WBA';
[~,I.Free,N.Free,~,~,FILES,PATH] = GetFileData(root);

VID.Free = uint8(zeros(630,840,3,N.Free.file));
WBA.Free = nan(N.Free.file,5);
vid.Free = cell(N.Free.fly,1);
pp = 1;
prev = 0;
for kk = 1:N.Free.file
    load([PATH FILES{kk}],'WBA_image','Wing','med')
    VID.Free(:,:,:,kk) = WBA_image.cdata;
    WBA.Free(kk,1) = 90 + Wing.L;
    WBA.Free(kk,2) = 90 + Wing.R;
    
    if (I.Free.fly(kk) - prev)==0
        pp = pp + 1;
    else
        pp = 1;
    end
    prev = I.Free.fly(kk);
    
    vid.Free{I.Free.fly(kk)}(:,:,pp) = med;
end
WBA.Free(:,3) = WBA.Free(:,1) + WBA.Free(:,2);
WBA.Free(:,4) = WBA.Free(:,1) - WBA.Free(:,2);
WBA.Free(:,5) = I.Free.fly;
WBA_Free_med = [1*ones(N.Free.fly,1), accumarray(I.Free.fly, WBA.Free(:, 3), [], @median)];

vid.Free_med = cellfun(@(x) median(x,3), vid.Free, 'Uniformoutput', false);
vid.Free_all = median(cat(3,vid.Free_med{:}),3);

root = 'H:\EXPERIMENTS\Experiment_SOS_v2_HeadFixed\WBA';
% root = 'H:\EXPERIMENTS\Experiment_ChirpLog_HeadFixed\Vid\WBA';
[~,I.Fixed,N.Fixed,~,~,FILES,PATH] = GetFileData(root);

VID.Fixed = uint8(zeros(630,840,3,N.Fixed.file));
WBA.Fixed = nan(N.Fixed.file,5);
vid.Fixed = cell(N.Fixed.fly,1);
pp = 1;
prev = 0;
for kk = 1:N.Fixed.file
    load([PATH FILES{kk}],'WBA_image','Wing','med')
    VID.Fixed(:,:,:,kk) = WBA_image.cdata;
    WBA.Fixed(kk,1) = 90 + Wing.L;
    WBA.Fixed(kk,2) = 90 + Wing.R;
    
    if (I.Fixed.fly(kk) - prev)==0
        pp = pp + 1;
    else
        pp = 1;
    end
    prev = I.Fixed.fly(kk);
    
    vid.Fixed{I.Fixed.fly(kk)}(:,:,pp) = med;
end
WBA.Fixed(:,3) = WBA.Fixed(:,1) + WBA.Fixed(:,2);
WBA.Fixed(:,4) = WBA.Fixed(:,1) - WBA.Fixed(:,2);
WBA.Fixed(:,5) = I.Fixed.fly;
WBA_Fixed_med = [2*ones(N.Fixed.fly,1), accumarray(I.Fixed.fly, WBA.Fixed(:, 3), [], @median)];

vid.Fixed_med = cellfun(@(x) median(x,3), vid.Fixed, 'Uniformoutput', false);
vid.Fixed_all = median(cat(3,vid.Fixed_med{:}),3);

WBA_ALL = cat(1,WBA_Free_med,WBA_Fixed_med);

%% WBA Median Boxplot
FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 3 2];
movegui(FIG,'center')

ax1 = gca;
ax1.FontSize = 8;
bx = boxplot(ax1,WBA_ALL(:,2),WBA_ALL(:,1),'Labels',{'Free','Fixed'});

ax1.YLim = [0 300];
ax1.YTick = 0:50:300;

ylabel(['\Sigma WBA (' char(176) ')'],'FontSize',8)
h = get(bx(5,:),{'XData','YData'});
CC = [ 0.7 0 0 ; 0.35 0 0.65];
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end
set(findobj(gcf,'tag','Median'), 'Color', 'w','LineWidth',1);
set(findobj(gcf,'tag','Box'), 'Color', 'k','LineWidth',1);
set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k','LineWidth',1);
set(findobj(gcf,'tag','Lower Whisker'), 'Color', 'k','LineWidth',1);

ax1.Children = ax1.Children([end 1:end-1]);

%% WBA All flies
FIG = figure (2) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4 4];
movegui(FIG,'center')

ax1 = subplot(2,1,1);
    ax1.FontSize = 8;
    bx = boxplot(ax1,WBA.Free(:,3),I.Free.fly) ; title('Free')

    ylabel(['WBA (' char(176) ')'],'FontSize',8)
    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},[0.7 0 0]);
    end
    set(findobj(gcf,'tag','Median'), 'Color', 'w','LineWidth',1);
    set(findobj(gcf,'tag','Box'), 'Color', [0.5 0.5 0.5]','LineWidth',1);
    set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k','LineWidth',1);
    set(findobj(gcf,'tag','Lower Whisker'), 'Color', 'k','LineWidth',1);

    ax1.Children = ax1.Children([end 1:end-1]);
    ax1.XLim = [0 11];
    ax1.YLim = [0 350];

ax1 = subplot(2,1,2);
    ax1.FontSize = 8;
    bx = boxplot(ax1,WBA.Fixed(:,3),I.Fixed.fly) ; title('Fixed')

    ylabel(['WBA (' char(176) ')'],'FontSize',8)
    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},[0 0 0.5]);
    end
    set(findobj(gcf,'tag','Median'), 'Color', 'w','LineWidth',1);
    set(findobj(gcf,'tag','Box'), 'Color', [0.5 0.5 0.5]','LineWidth',1);
    set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k','LineWidth',1);
    set(findobj(gcf,'tag','Lower Whisker'), 'Color', 'k','LineWidth',1);

    ax1.Children = ax1.Children([end 1:end-1]);
    ax1.XLim = [0 11];
    ax1.YLim = [0 350];
%     xlabel('Fly #')

%% Free vs Fixed Image
FIG = figure (3) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4 4];
movegui(FIG,'center')
colormap(jet)
% colormap(hsv)

subplot(2,1,1)
frame = 10*vid.Free_med{7};
imagesc(frame)
axis off equal
title('Free')

subplot(2,1,2)
frame = 10*vid.Fixed_med{3};
imagesc(frame)
axis off equal
title('Fixed')

end