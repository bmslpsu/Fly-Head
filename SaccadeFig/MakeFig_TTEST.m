function [] = MakeFig_TTEST()
%% MakeFig_TTEST:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE.Anti,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE.Anti = cellstr(FILE.Anti)';

Anti = load(fullfile(root,FILE.Anti{1}),'SACD','U','N','T');

Vel = Anti.U{1,3}{1};

%% ANOVA %%

Y = Anti.SACD.Head.Position_IntError;
Y = Anti.SACD.Head.Velocity_IntError;

% G1 = Anti.SACD.Head.fly;
G2 = Anti.SACD.Head.speed;
% G3 = Anti.SACD.Head.Dir;
% G4 = Anti.SACD.Head.speed;
% G5 = Anti.SACD.Head.speed;

P = anovan(Y,{G2})




end