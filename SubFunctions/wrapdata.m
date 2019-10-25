function [wrap,n_lim] = wrapdata(raw,lim,debug)
%% wrapdata:  wraps data to specified limit
%   INPUT:
%       raw     : raw data
%       lim     : wrap limit
%       debug   : boolean (showplots)
%   OUTPUT:
%       wrap   	: wrapped data
%       n_lim   : # of time data hits limit
%

if nargin<3
    debug = false;
end

wrap = raw;
if debug
    figure
    hold on
    plot(wrap)
end

n_lim = 0;
go = true;
while any(go)   
    above = wrap > lim;
    go(1) = any(above);
    
    below = wrap < -lim;
    go(2) = any(below);
    
    wrap(above) = wrap(above) - lim;
    wrap(below) = wrap(below) + lim;
    
  	if debug
        if ~any(go)
            plot(wrap,'k','LineWidth',1.5)
        else
            plot(wrap)
        end
    end
    
    n_lim = n_lim + 1;
end
end