function [degData] = panel2deg(voltData)
%% panel2deg: Reads voltage from arena controller & converts to angles
%   INPUTS:
%       voltData    :   raw voltage ouput of arena
%   OUTPUTS:
%       degData     :   angular position of arena pattern  
%---------------------------------------------------------------------------------------------------------------------------------
vRange = 10;    % ouput voltage range
pixels = 96;    % # of pixels around display
spRes  = 3.75;  % deg/pixel
voltData = voltData - mean(voltData); % get rid of offset to center around 0
panelData = round((pixels/vRange)*voltData); % panel position [panel index]
degData = spRes*panelData; % panel angular position [deg]
end