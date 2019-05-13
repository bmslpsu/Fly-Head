function [FIG] = Run_MakeFig_Sine_HeadFree_pat2head_BODE_ALL()
%% Run_MakeFig_Sine_HeadFree_pat2head_BODE_ALL:
%   INPUTS:
%       -
%   OUTPUTS:
%    	-
%---------------------------------------------------------------------------------------------------------------------------------
CLR = {'b','r','g','c','k'};
saveFig = false;
for kk = 1:length(CLR)
    [FIG] = MakeFig_Sine_HeadFree_pat2head_BODE_ALL(CLR{kk},saveFig);
end

figdir = 'H:\DATA\Rigid_Data\FIGURE\';
saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
end