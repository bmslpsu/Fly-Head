%% Run Make_Data functions %%
%---------------------------------------------------------------------------------------------------------------------------------
clear ; close all ; clc
rootdir = 'H:\EXPERIMENTS\Experiment_Sinusoid\';
Amp = 3.75*(1:5);
for kk = 1:length(Amp)
    [~] = MakeData_Sine_HeadFree(rootdir,Amp(kk));
    beep on
    for jj = 1:3
        beep
        pause(1)
    end
end