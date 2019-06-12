[EMD_ALL, Pattern_ALL, Eye_ALL] = EMDsim_v5_7_5deg_Grnd_new(true, 0);

[EMD_ALL_filt] = HR_sim(Eye_ALL);

EMD_ALL_filt = EMD_ALL_filt(:,1:end-1); % resize array 1-63

% make space-time plot
figure; subplot(1,2,1); imagesc(Pattern_ALL)
ylabel('time (s)'); set(gca,'YTick',[0 2.0*160])
set(gca,'YTickLabel',{'0', '2.0'})
set(gca,'XTick',[1 48 96])
set(gca,'XTickLabel',{'-180','0', '180'})
title('Intensity space-time plot');
xlabel('Pattern position'); colormap bone; freezeColors;
ylim([0 2.0*160]);

% EMD
subplot(1,2,2); 
imagesc(EMD_ALL_filt);
set(gca,'YTick',[0 2.0*160]); set(gca,'YTickLabel',{'0', '2.0'})
colormap(jet);
xlabel('ommatidia'); title('EMD'); ylim([0 2.0*160]); xlim([1 8]);
set(gca,'XTick',[1 8])
set(gca,'XTickLabel',{'1','8'})