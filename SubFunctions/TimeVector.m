function [time_out,tt] = TimeVector(time_in,linFlag)
%% TimeVector: makes chirp position function
%   INPUTS:
%       time_in:      	:   time vector from experimental equipment
%       linFlag        	:   true for linspace distrobution
%   OUTPUTS:
%       time_out        :   evenly spaced output vector
%---------------------------------------------------------------------------------------------------------------------------------
if nargin==1
    linFlag = false;
end

tt.n      = length(time_in);
tt.min    = round(min(time_in),2);
tt.max    = round(max(time_in),2);
tt.T      = range([tt.min tt.max]);
tt.fs     = round(1/mean(diff(time_in)),2);
tt.ts     = 1/tt.fs;

if linFlag
    time_out = linspace(tt.min,tt.max,tt.T*tt.fs)';
else
    time_out = (tt.min:tt.ts:tt.max)';
end

end