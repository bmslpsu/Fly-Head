function [panel] = deg2panel(deg)
%% MakePosFunction_Chirp: makes chirp position function
%   INPUTS:
%       deg     :   data in degrees
%   OUTPUTS:
%       panel   :   data in panel #
%---------------------------------------------------------------------------------------------------------------------------------=
panel = round(deg/3.75); % convert to steps [deg]
end