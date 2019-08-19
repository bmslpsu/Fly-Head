function [tout,iout,tdiff] = Tindex(tin,tfind)
%% Tindex: find the value and index of a scalr value in a vector
%   INPUTS:
%       tin  	: vector to search
%       tfind  	: scalar to search for
%   OUTPUTS:
%       tout  	: closest value to "tfind"
%       iout 	: index of "tout" in "tin"
%---------------------------------------------------------------------------------------------------------------------------------
[tdiff,iout] = min(abs(tin - tfind));
tout = tin(iout);
end

