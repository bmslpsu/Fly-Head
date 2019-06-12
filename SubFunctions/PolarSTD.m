function [rSTD] = PolarSTD(X,Y,sigma)
%% PolarSTD:
%   INPUTS:
%       x       :   x-data
%       y       :   y-data
%       sigma   :   mean (X,Y), default is (0,0)
%   OUTPUTS:
%       mu      : standard deviation
%---------------------------------------------------------------------------------------------------------------------------------
if nargin<3
    sigma = [0,0]; % default mean is (0,0)
end

% Remove missing data
X = rmmissing(X);
Y = rmmissing(Y);

if size(X)~=size(Y)
    error('Error: X & X vectors must be the same size')
else
    X = X(:);
    Y = Y(:);
    nPoint = length(X);
end

% Normalized to mean
nX = X - sigma(1);
nY = Y - sigma(2);

% Find vectors from mean to points
R = sqrt(nX.^(2) + nY.^(2));

% Calculate standard deviation of vectors
rSTD = std(R);
end