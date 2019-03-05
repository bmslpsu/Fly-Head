function [angMean] = angle_mean(angData,dim)
%% angle_mean: calculates mean of circular angle data
%   INPUTS:
%       angData     : circular angle data
%       dim         : dimension to operate on
%   OUTPUTS:
%       angMean     : circular mean
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUTS %
% clear ; clc
% dim = 2;
% angData = [ 10 20 30 40 ;...
%             40 30 20 10  ];
%---------------------------------------------------------------------------------------------------------------------------------
wrapData = wrapTo360(angData);
n = size(angData,dim);
s = (1/n)*sum(sind(wrapData),dim);
c = (1/n)*sum(cosd(wrapData),dim);

if dim==1
    dir = size(angData,2);
    angMean = nan(1,size(angData,2));
elseif dim==2
    dir = size(angData,1);
    angMean = nan(size(angData,1),1);
elseif dim==3
	dir = size(angData,3);
	angMean = nan(size(angData));
else
    error('WUT')
end

for kk = 1:dir
    s_bar = s(kk);
    c_bar = c(kk);
    if (s_bar>0) && (c_bar>0)
        angMean(kk) = atand(s_bar/c_bar);
    elseif (c_bar<0)
        angMean(kk) = atand(s_bar/c_bar) + 180;
    elseif (s_bar<0) && (c_bar>0)
        angMean(kk) = atand(s_bar/c_bar) + 360;
    else
        error('WUT')
    end
end

end
