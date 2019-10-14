function [panel] = deg2panel(deg,res)
%% MakePosFunction_Chirp: makes chirp position function
%   INPUTS:
%       deg     :   data in degrees
%   OUTPUTS:
%       panel   :   data in panel #
%---------------------------------------------------------------------------------------------------------------------------------=
if nargin==1
    res = 3.75; % default
end
panel = round(deg/res); % convert to steps [deg]
end